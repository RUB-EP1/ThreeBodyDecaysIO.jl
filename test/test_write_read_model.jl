using ThreeBodyDecaysIO, ThreeBodyDecays, Test
using ThreeBodyDecaysIO.HadronicLineshapes
using ThreeBodyDecaysIO: JSON, FourVectors

model = let
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )
    chain = DecayChain(
        k = 1,
        two_j = 2,
        Xlineshape = BreitWigner(1.1, 0.1),
        Hij = Vertex(RecouplingLS((2, 0))),
        HRk = Vertex(RecouplingLS((2, 2))),
        tbs = tbs,
    )
    ThreeBodyDecay(
        ["first"=>(1.0, chain), "second"=>(2.0, DecayChain(chain; k = 2))];
        reference_k = 3,
    )
end
@testset "Indexed model roundtrip" begin
    document = JSON.parse(JSON.json(add_hs3_fields(model)))
    parsed = read_model(document)
    @test reference_k(parsed.model) == 3
    @test parsed.model.names == model.names
    @test parsed.model.couplings == model.couplings
    @test !haskey(
        document["distributions"][1]["decay_description"]["kinematics"]["initial_state"],
        "mass",
    )
    sigma = randomPoint(masses(model))
    @test amplitude(parsed.model, sigma) ≈ amplitude(model, sigma)
    chain = first(model.chains)
    shared = ThreeBodyDecay(["one"=>(1.0, chain), "two"=>(2.0, chain)]; reference_k = 2)
    description, functions = serializeToDict(shared)
    @test length(functions) == 1
    @test getindex.(description["chains"], "name") == ["one", "two"]
    @test read_model(add_hs3_fields(shared)).model.names == shared.names
end

@testset "Referenced vertex form-factor roundtrip" begin
    chain = first(model.chains)
    with_factor = DecayChain(
        k = chain.k,
        two_j = chain.two_j,
        Xlineshape = chain.Xlineshape,
        HRk = chain.HRk,
        Hij = Vertex(chain.Hij.h, BlattWeisskopf{1}(2.0)),
        tbs = chain.tbs,
    )
    source = ThreeBodyDecay(["factor" => (1.0, with_factor)]; reference_k = 3)
    document = JSON.parse(JSON.json(add_hs3_fields(source)))
    parsed = read_model(document).model
    sigma = randomPoint(masses(source))
    @test amplitude(parsed, sigma) ≈ amplitude(source, sigma)
    @test any(fn -> fn["type"] == "BlattWeisskopf", document["functions"])
end
