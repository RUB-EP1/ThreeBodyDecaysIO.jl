using ThreeBodyDecaysIO.Random
Random.seed!(3441)
rotate(p, a, b, g) = p |> FourVectors.Rz(g) |> FourVectors.Ry(b) |> FourVectors.Rz(a)
@testset "Independent geometry reconstruction, poles, canonical boost" begin
    ms=ThreeBodyMasses(1, 2, 3; m0 = 7.9)
    for n = 1:20, k = 1:3
        s=randomPoint(ms)
        aligned=map(
            t->FourVectors.FourVector(t[1], t[2], t[3]; E = t[4]),
            aligned_four_vectors(s, ms; k),
        )
        for (a, b, g) in (
            (0.78, 0.62, -0.51),
            (-2.4, 2.7, 2.9),
            (1.2, 0.0, -0.9),
            (1.2, Float64(pi), -0.9),
            (0.0, 0.0, 0.0),
        )
            event=map(p->rotate(p, a, b, g), aligned)
            angles, inv=angles_invariants(event, ms; k)
            rebuilt=map(p->rotate(p, angles.α, acos(angles.cosβ), angles.γ), aligned)
            @test all(
                isapprox(collect(p), collect(q); atol = 1e-10) for
                (p, q) in zip(event, rebuilt)
            )
            @test collect(inv) ≈ collect(s) atol=1e-10
            moving=map(event) do p
                p |>
                FourVectors.Rz(-0.39) |>
                FourVectors.Ry(-0.85) |>
                FourVectors.Bz(1.7) |>
                FourVectors.Ry(0.85) |>
                FourVectors.Rz(0.39)
            end
            restored=map(p->ThreeBodyDecaysIO.canonical_rest_vector(p, sum(moving)), moving)
            @test all(
                isapprox(collect(p), collect(q); atol = 1e-10) for
                (p, q) in zip(event, restored)
            )
            ang2, inv2=angles_invariants(moving, ms; k, frame = :canonical)
            rebuilt2=map(p->rotate(p, ang2.α, acos(ang2.cosβ), ang2.γ), aligned)
            @test all(
                isapprox(collect(p), collect(q); atol = 1e-7) for
                (p, q) in zip(event, rebuilt2)
            )
        end
    end
end

@testset "Physical input and reference rejection" begin
    ms = masses(model)
    vectors = map(FourVectors.FourVector, aligned_four_vectors(randomPoint(ms), ms; k = 3))
    for k in (true, 1.0, 0, 4)
        @test_throws ArgumentError angles_invariants(vectors, ms; k)
    end
    @test_throws MethodError angles_invariants(vectors, model; k = 1)
    @test_throws MethodError angles_invariants(vectors, model; reference_k = 1)
    @test angles_invariants(vectors, model) == angles_invariants(vectors, ms; k = 3)
    domain = PhaseSpaceDomain("default", ms)
    point = four_vector_point(vectors, domain)
    for component in ("E", "px")
        bad = deepcopy(point)
        bad["final_state"][1]["p"][component] += 1
        @test_throws ArgumentError four_vectors(bad, domain)
    end
    for energy in (0.0, -1.0, Inf, NaN)
        bad = deepcopy(point)
        bad["final_state"][1]["p"]["E"] = energy
        @test_throws ArgumentError four_vectors(bad, domain)
    end
    moving = map(p -> p |> FourVectors.Bz(1.3), vectors)
    @test_throws ArgumentError angles_invariants(moving, ms; k = 3)
    for index in (0, 4, true, 1.0)
        bad=deepcopy(point)
        bad["final_state"][1]["index"]=index
        @test_throws ArgumentError four_vectors(bad, domain)
    end
    bad=deepcopy(point)
    pop!(bad["final_state"])
    @test_throws ArgumentError four_vectors(bad, domain)
    description=first(serializeToDict(domain))
    for mass in (-1.0, Inf, NaN, 0.0)
        bad=deepcopy(description)
        bad["initial_state"]["mass"]=mass
        @test_throws ArgumentError dict2instance(PhaseSpaceDomain, bad)
    end
    @test_throws ArgumentError HadronicUnpolarizedIntensity(
        model,
        PhaseSpaceDomain("other", ThreeBodyMasses(0.1, 0.2, 0.3; m0 = 4.0)),
    )
    document=add_hs3_fields(model)
    push!(document["distributions"], deepcopy(first(document["distributions"])))
    @test_throws ArgumentError read_model(document)
end

@testset "Polarized orientation roundtrip" begin
    ms=masses(model)
    for k = 1:3
        referenced=ThreeBodyDecay(
            model.names .=> zip(model.couplings, model.chains);
            reference_k = k,
        )
        sigma=randomPoint(ms)
        aligned=map(FourVectors.FourVector, aligned_four_vectors(sigma, ms; k))
        expected=(α = 0.34, cosβ = 0.42, γ = -0.73)
        event=map(p->rotate(p, expected.α, acos(expected.cosβ), expected.γ), aligned)
        recovered, invariants=angles_invariants(event, ms; k)
        @test amplitude(referenced, recovered, invariants) ≈
              amplitude(referenced, expected, sigma)
        @test unpolarized_intensity(referenced, recovered, invariants) ≈
              sum(abs2, amplitude(referenced, expected, sigma))
    end
end
