using ThreeBodyDecaysIO: parse_into_function
using ThreeBodyDecaysIO
using Parameters
using Test

e(x) = 1 / (0.77^2 - x - 1im * 0.77 * 0.15 * sin(x))
#
expression = "1/(0.77^2 - x - i*0.77*0.15*sin(x))"

@test expression_argument(Meta.parse(expression)) == :x
@test expression_argument(Meta.parse(replace(expression, "x" => "m12sq"))) == :m12sq

parse_into_function(expression)

@testset "Parse into a function works" begin
    f, _ = parse_into_function(expression)
    @test f(1.1) == e(1.1)
    #
    f, _ = parse_into_function(expression * "*8")
    @test f(1.1) == e(1.1) * 8
    #
    @test_throws ArgumentError parse_into_function(expression * "*y")
end

@testset "BreitWignerWidthExpLikeBugg as a string" begin

    K700_BuggBW_lineshape = """
    1/(0.824^2 - σ - i * 0.824 *
        (σ - 0.23397706275638377) / (0.824^2 - 0.23397706275638377) * 0.478 * exp(-0.941060 * σ))
    """

    X_K700_BuggBW, _ = parse_into_function(K700_BuggBW_lineshape)
    @test X_K700_BuggBW(1.1) ≈ -1.6748633697157216 + 1.083006788846748im


    K1430_BuggBW_lineshape = """
    1/(1.375^2 - σ - i * 1.375 *
        (σ - 0.23397706275638377) / (1.375^2 - 0.23397706275638377) * 0.190 * exp(-0.020981 * σ))
    """
    X_K1430_BuggBW, _ = parse_into_function(K1430_BuggBW_lineshape)
    @test X_K1430_BuggBW(1.1) ≈ 1.2297831004988296 + 0.20758230110948936im
end

@testset "KatchaevSigma as a string" begin
    KatchaevSigma_lineshape = """
    1.0 / (
        0.1131 / (σ + 0.0073999) + 0.0337 -
        0.3185 * (σ / 0.982657846681 - 1.0) -
        0.0942 * (σ / 0.982657846681 - 1.0)^2 -
        0.5927 * (σ / 0.982657846681 - 1.0)^3 -
            i * sqrt(1 - 4*0.13956755^2 / σ)
        )
    """
    KatchaevSigma, _ = parse_into_function(KatchaevSigma_lineshape)
    @test KatchaevSigma(1.1) ≈ 0.10172436035046228 + 1.0273440332286132im
end


let
    dict = Dict(
        "type" => "generic_function",
        "expression" => "1/(0.77^2 - x - i*0.77*0.15*sin(x))",
    )
    @unpack type = dict
    instance_type = eval(Symbol(type))
    custom_function = dict2instance(instance_type, dict)
    @testset "generic_function works" begin
        @test custom_function(Dict("x" => 1.1)) == e(1.1)
    end
end


@testset "generic with arguments" begin
    d = Dict(
        "type" => "generic_function",
        "expression" => "(m_12 + a) / (m_12 + b)",
        "a" => 1.0,
        "b" => -3.0,
    )
    f = dict2instance(generic_function, d)
    @test f.f(1.1) ≈ (1.1 + 1) / (1.1 + -3)
    #
    d = Dict("type" => "generic_function", "expression" => "(m_12 + 1) / (m_12 - 3)")
    g = dict2instance(generic_function, d)
    @test g.f(2.2) ≈ -4
    #
    d = Dict(
        "type" => "generic_function",
        "name" => "X(1001)",
        "expression" => "x+l*log(x)-2*o+g",
        "l" => 2.0,
        "o" => -3.0,
        "g" => 1.5,
    )
    g = dict2instance(generic_function, d)
    @test g(Dict("x" => 3)) ≈ 3 + 2 * log(3) - 2 * (-3) + 1.5
end
