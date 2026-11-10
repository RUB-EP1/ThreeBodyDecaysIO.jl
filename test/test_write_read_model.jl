using ThreeBodyDecaysIO.HadronicLineshapes
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using ThreeBodyDecaysIO
using Test


model = let
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )
    dpp = randomPoint(tbs)
    #
    two_j = 2
    #
    ch1 = DecayChain(;
        k = 1,
        two_j,
        Xlineshape = BreitWigner(1.1, 0.1),
        Hij = RecouplingLS((two_j, 0)) |> VertexFunction,
        HRk = RecouplingLS((two_j, two_j)) |> VertexFunction,
        tbs,
    )
    ch2 = DecayChain(ch1; k = 2)
    ch3 = DecayChain(ch1; k = 3)
    #
    ThreeBodyDecay("K892" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)])
end

decay_description, appendix = serializeToDict(model)
dict = add_hs3_fields(decay_description, appendix, "default-model")

open("test.json", "w") do io
    JSON.print(io, dict, 4)
end

# reading
json_content = open("test.json") do io
    JSON.parse(io)
end
rm("test.json")

@unpack decay_description = json_content["distributions"][1]
@testset "JSON has all key sections" begin
    @test haskey(decay_description, "kinematics")
    @test haskey(decay_description, "reference_topology")
    @test haskey(decay_description, "chains")
    #
    # @test haskey(json_content, "variables")
    # @test haskey(json_content, "validation")
    #
    @unpack chains = decay_description
    chain = chains[1]
    @test haskey(chain, "vertices")
    @test haskey(chain, "propagators")
end

@testset "Specific to the test model" begin
    @test length(decay_description["chains"]) == 3

    @unpack chains = decay_description
    map(chains) do chain
        @test length(chain["vertices"]) == 2
        @test length(chain["propagators"]) == 1
    end
end


input = copy(json_content)

@testset "Parse model" begin
    @unpack decay_description = input["distributions"][1]
    @unpack functions = input
    workspace = Dict{String,Any}()
    for fn in functions
        @unpack name, type = fn
        instance_type = eval(Symbol(type))
        workspace[name] = dict2instance(instance_type, fn)
    end
    @test dict2instance(ThreeBodyDecay, decay_description; workspace) isa ThreeBodyDecay
end



@testset "Parse model matches the original model" begin
    @unpack decay_description = input["distributions"][1]
    @unpack functions = input
    workspace = Dict{String,Any}()
    for fn in functions
        @unpack name, type = fn
        instance_type = eval(Symbol(type))
        workspace[name] = dict2instance(instance_type, fn)
    end
    model_parsed = dict2instance(ThreeBodyDecay, decay_description; workspace)

    σs0 = randomPoint(masses(model_parsed))
    @test unpolarized_intensity(model_parsed, σs0) ≈ unpolarized_intensity(model, σs0)
end

@testset "Chain names are independent of shared propagators" begin
    chain = first(model.chains)
    named_model = ThreeBodyDecay(["first" => (1.0, chain), "second" => (2.0, chain)])
    description, functions = serializeToDict(named_model)
    @test getindex.(description["chains"], "name") == ["first", "second"]
    @test length(functions) == 1
    @test all(
        c["name"] != c["propagators"][1]["parametrization"] for c in description["chains"]
    )
    restored = dict2instance(ThreeBodyDecay, description; workspace = functions)
    @test restored.names == named_model.names
    @test restored.couplings == named_model.couplings
    @test first(serializeToDict(chain))["name"] == "my_decay_chain"
end
