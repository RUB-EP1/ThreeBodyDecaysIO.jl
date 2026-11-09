using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.HadronicLineshapes
using ThreeBodyDecaysIO.OrderedCollections
using Test


@testset "BreitWigner from plane Dict" begin
    d = Dict(
        "type" => "BreitWigner",
        "name" => "L1520_BW",
        "mass" => 1.0,
        "width" => 0.1,
        "ma" => 0.0,
        "mb" => 0.0,
        "l" => 0,
        "d" => 1.5,
        "x" => "m23sq",
    )
    bw1 = dict2instance(BreitWigner, d)
    @test bw1 isa NamedArgFunc{<:HadronicLineshapes.AbstractFlexFunc}
    @test bw1(Dict("m23sq" => 1.1)) ≈ -5 + 5im
    @test bw1(LittleDict("m23sq" => 1.1)) ≈ -5 + 5im
    @test bw1(OrderedDict("m23sq" => 1.1)) ≈ -5 + 5im
    @test_throws KeyError bw1(LittleDict("msq" => 1.1))
end

@testset "MomentumPower from plane Dict" begin
    d = Dict("type" => "MomentumPower", "l" => 5)
    ff1 = dict2instance(MomentumPower, d)
    @test ff1 isa HadronicLineshapes.MomentumPower
    @test ff1(2.0) == 32
end

@testset "BlattWeisskopf from plane Dict" begin
    d = Dict("type" => "BlattWeisskopf", "l" => 0, "radius" => 1.5)
    bw1 = dict2instance(BlattWeisskopf, d)
    @test bw1 isa HadronicLineshapes.BlattWeisskopf
end

@testset "MultichannelBreitWigner from a nested Dict" begin
    d = LittleDict(
        "type" => "MultichannelBreitWigner",
        "name" => "L1520_BW",
        "mass" => 1.0,
        "channels" =>
            [LittleDict("gsq" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)],
        "x" => "msq",
    )
    bw1 = dict2instance(MultichannelBreitWigner, d)
    @show typeof(bw1)
    @test bw1 isa NamedArgFunc{<:HadronicLineshapes.AbstractFlexFunc}
end

@testset "Deserialize MultichannelBreitWigner" begin
    d = LittleDict(
        "type" => "MultichannelBreitWigner",
        "name" => "L1520_BW",
        "mass" => 1.0,
        "x" => "msq",
        "channels" => [
            LittleDict("gsq" => 1.1, "ma" => 0.1, "mb" => 0.2, "l" => 1, "d" => 1.5),
            LittleDict("gsq" => 2.1, "ma" => 1.0, "mb" => 2.0, "l" => 3, "d" => 1.5),
        ],
    )
    bw1 = dict2instance(MultichannelBreitWigner, d)
    dict, _ = serializeToDict(bw1)
    @test all(values(dict["channels"][1]) .== values(d["channels"][1]))
    @test all(values(dict["channels"][2]) .== values(d["channels"][2]))
end

@testset "Polynomial from plane Dict" begin
    d = Dict("coefficients" => [1.0, 2.0, 3.0], "x" => "m23sq")
    bw1 = dict2instance(Polynomial, d)
    @test bw1.f(1) == 6
    fsq = bw1.f * WrapFlexFunction(x -> 3x)
    @test fsq(1) == 18
    @test bw1(Dict("m23sq" => 1)) == 6
    @test bw1 isa NamedArgFunc{<:HadronicLineshapes.AbstractFlexFunc}
end

@testset "Named function round trips" begin
    for (T, data) in (
        (
            BreitWigner,
            Dict(
                "mass" => 1.0,
                "width" => 0.1,
                "ma" => 0.0,
                "mb" => 0.0,
                "l" => 0,
                "d" => 1.5,
                "x" => "s",
            ),
        ),
        (
            MultichannelBreitWigner,
            Dict(
                "mass" => 1.0,
                "x" => "s",
                "channels" =>
                    [Dict("gsq" => 0.1, "ma" => 0.0, "mb" => 0.0, "l" => 0, "d" => 1.5)],
            ),
        ),
        (Polynomial, Dict("coefficients" => [1.0, 2.0, 3.0], "x" => "s")),
    )
        wrapped = dict2instance(T, data)
        serialized, _ = serializeToDict(wrapped)
        restored = dict2instance(T, serialized)
        @test restored.variable_names == ["s"]
        for value in (0.5, 1.1, 2.0)
            @test restored(Dict("s" => value, "unused" => 42)) ≈ wrapped.f(value)
        end
        @test_throws KeyError restored(Dict("wrong" => 1.1))
    end
    @test_throws ArgumentError NamedArgFunc(identity, String[])
    @test_throws ArgumentError NamedArgFunc(identity, ["s", "t"])
end
