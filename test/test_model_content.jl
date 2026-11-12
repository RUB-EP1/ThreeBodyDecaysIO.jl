using HadronicLineshapes
using ThreeBodyDecaysIO
using ThreeBodyDecays
using Parameters
using DataFrames
using JSON
using Test



test_file_name = joinpath(@__DIR__, "..", "test", "lc2ppi-lhcb-test.json")

let
    @info "⭐ Reading model from $test_file_name ⭐"
    # @testset "Test of $test_file_name" begin
    @testset "Test of $test_file_name" begin
        # get the JSON content
        input = open(test_file_name) do io
            JSON.parse(io)
        end
        model_description = first(input["distributions"])
        @test model_description["type"] == "HadronicUnpolarizedIntensity"

        @test haskey(model_description, "decay_description")
        @unpack decay_description = model_description

        @test haskey(decay_description, "reference_topology")
        @unpack reference_topology = decay_description

        @test flatten_topology(reference_topology) |> sort == [1, 2, 3]

        @info "🐰 Parsing kinematics 🐰"
        @unpack kinematics = decay_description
        tbs = dict2instance(ThreeBodySystem, kinematics)

        @test haskey(input, "functions")
    end
end



# More serious stuff
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# extra lineshapes for testing
@with_kw struct BreitWignerWidthExpLikeBugg <: HadronicLineshapes.AbstractFlexFunc
    m::Float64
    Γ::Float64
    γ::Float64
end
function (bw::BreitWignerWidthExpLikeBugg)(σ)
    mK = 0.493677
    mπ = 0.13957018
    σA = mK^2 - mπ^2 / 2
    @unpack m, Γ, γ = bw
    Γt = (σ - σA) / (m^2 - σA) * Γ * exp(-γ * σ)
    1 / (m^2 - σ - 1im * m * Γt)
end
function ThreeBodyDecaysIO.dict2instance(::Type{BreitWignerWidthExpLikeBugg}, dict)
    @unpack mass, width, slope, x = dict
    bw = BreitWignerWidthExpLikeBugg(mass, width, slope)
    variables = [x]
    return NamedArgFunc(bw, variables)
end

let
    @info "⭐ Reading model from $test_file_name ⭐"
    # @testset "Test of $test_file_name" begin
    @testset "Test of $test_file_name" begin
        # get the JSON content
        input = open(joinpath(@__DIR__, "..", "models", test_file_name)) do io
            JSON.parse(io)
        end

        @info "🔥 Building lineshapes functions 🔥"
        workspace = Dict{String,Any}()
        @unpack functions = input
        for fn in functions
            @unpack name, type = fn
            @info "   🎈$name"
            instance_type = eval(Symbol(type))
            workspace[name] = dict2instance(instance_type, fn)
        end
        # @test dict2instance(ThreeBodyDecay, decay_description; workspace) isa ThreeBodyDecay

        @info "🦊 Building distributions functions 🦊"
        # build distributions from JSON array,
        @unpack distributions = input
        for dist in distributions
            @unpack name, type = dist
            @info "   🎈$name"
            instance_type = eval(Symbol(type))
            workspace[name] = dict2instance(instance_type, distributions[1]; workspace)
        end

        @info "👡 Performing validation 👡"
        @unpack misc, parameter_points = input
        @unpack amplitude_model_checksums = misc

        @testset "Validation" begin
            map(amplitude_model_checksums) do check_point_info
                @unpack name, value, distribution = check_point_info
                #
                # pull distribution
                dist = workspace[distribution]

                # pull correct parameter point
                parameter_points_dict = array2dict(parameter_points; key = "name")
                # find the point in the list of points
                parameter_point = parameter_points_dict[name]
                # compute, compare
                _parameters = array2dict(
                    parameter_point["parameters"];
                    key = "name",
                    apply = v -> v["value"],
                )
                @test ThreeBodyDecays.reference_k(dist.model) == 1
                mass_angles = ThreeBodyDecaysIO.mass_angles_replace(
                    _parameters,
                    dist.mass_angles_cascade_names,
                )
                _, σs = angles_invariants(mass_angles, masses(dist.model); k = 1)
                @test σs.σ1 ≈ 0.7980703453578917
                @test σs.σ2 ≈ 3.6486261122281745
                @test value ≈ dist(_parameters)
                return "🟢"
            end
        end
    end
end
