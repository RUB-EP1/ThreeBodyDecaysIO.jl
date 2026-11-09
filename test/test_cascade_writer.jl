using Test
using ThreeBodyDecaysIO
import CascadeDecays as CD
import ThreeBodyDecaysIO.ThreeBodyDecays as TBD
using ThreeBodyDecaysIO.HadronicLineshapes

@testset "Cascade writer agrees with three-body writer" begin
    ms = TBD.ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09)
    tbs = TBD.ThreeBodySystem(ms, TBD.ThreeBodySpins(0, 0, 0; two_h0 = 2))
    for k = 1:3
        i, j = TBD.ij_from_k(k)
        tree = ((i, j), k)
        topology = CD.DecayTopology(tree)
        shape = BreitWigner(1.1, 0.1)
        inner = TBD.VertexFunction(TBD.RecouplingLS((2, 0)))
        outer = TBD.VertexFunction(TBD.RecouplingLS((2, 2)))
        original = TBD.DecayChain(;
            k,
            two_j = 2,
            Xlineshape = shape,
            Hij = inner,
            HRk = outer,
            tbs,
        )
        cascade_chain = CD.DecayChain(
            topology,
            CD.SystemSpins(0, 0, 0; two_h0 = 2);
            propagators = ((i, j) => CD.Propagator(2, shape),),
            vertices = (tree => outer, (i, j) => inner),
        )
        legacy = TBD.ThreeBodyDecay(["resonance" => (2.0 + 1.0im, original)])
        cascade = CD.CascadeDecay(
            (cascade_chain,),
            topology;
            couplings = (2.0 + 1.0im,),
            names = ("resonance",),
        )
        expected, expected_functions = serializeToDict(legacy)
        actual, functions = serializeToDict(cascade; masses = Tuple(ms))
        @test actual == expected
        @test add_hs3_fields(actual, functions) ==
              add_hs3_fields(expected, expected_functions)
        @test_throws ArgumentError serializeToDict(cascade; masses = (1.0,))
        @test_throws ArgumentError serializeToDict(cascade; masses = (0.1, 0.2, NaN, 1.0))
    end
end

@testset "Four-body cascade serialization" begin
    tree = ((1, 2), (3, 4))
    topology = CD.DecayTopology(tree)
    vertex = CD.Vertex(TBD.RecouplingLS((0, 0)))
    chain = CD.DecayChain(
        topology,
        CD.SystemSpins(0, 0, 0, 0; two_h0 = 0);
        propagators = (
            (1, 2) => CD.Propagator(0, CD.ConstantLineshape(1.0)),
            (3, 4) => CD.Propagator(0, CD.ConstantLineshape(2.0)),
        ),
        vertices = (tree => vertex, (1, 2) => vertex, (3, 4) => vertex),
    )
    model = CD.CascadeDecay((chain,), topology)
    dict, functions = serializeToDict(model; masses = (0.1, 0.2, 0.3, 0.4, 3.0))
    @test dict["reference_topology"] == [[1, 2], [3, 4]]
    @test length(dict["kinematics"]["final_state"]) == 4
    @test length(dict["chains"][1]["vertices"]) == 3
    @test length(dict["chains"][1]["propagators"]) == 2
    @test length(functions) == 2
    @test sort([first(serializeToDict(f))["coefficients"][1] for f in values(functions)]) == [1.0, 2.0]
end
