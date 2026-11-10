using ThreeBodyDecaysIO
using ThreeBodyDecays
using HadronicLineshapes
using Test

@testset "Symbol to String Key Conversion" begin
    # Create a simple model
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )

    ch1 = DecayChain(;
        k = 1,
        two_j = 2,
        Xlineshape = BreitWigner(1.1, 0.1),
        Hij = RecouplingLS((2, 0)) |> Vertex,
        HRk = RecouplingLS((2, 2)) |> Vertex,
        tbs,
    )

    model = ThreeBodyDecay("test" .=> [(1.0, ch1)]; reference_k = 1)
    dict, appendix = serializeToDict(model)

    # Test that all keys are Strings (not Symbols)
    @test all(isa(k, String) for k in keys(dict))
    @test all(isa(k, String) for k in keys(dict["kinematics"]))
    @test all(isa(k, String) for k in keys(dict["chains"][1]))
    @test all(isa(k, String) for k in keys(dict["chains"][1]["vertices"][1]))
    @test all(isa(k, String) for k in keys(dict["chains"][1]["propagators"][1]))
end
