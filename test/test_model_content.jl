@testset "Pinned PR111 reference models" begin
    for name in (
        "lb2pkg-lhcb-2765817.json",
        "lc2ppik-lhcb-2683025.json",
        "x2pipipi-compass-1391643.json",
    )
        results = validate(joinpath(@__DIR__, "reference_models", name))
        @test !isempty(results)
        @test all(r -> r.error === nothing, results)
        @test all(r -> r.passed, results)
    end
end
