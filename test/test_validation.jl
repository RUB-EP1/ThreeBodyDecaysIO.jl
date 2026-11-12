using Test
using ThreeBodyDecaysIO
using ThreeBodyDecays

@testset "Order in angles_invariants" begin
    ms = ThreeBodyMasses(1.0, 2.0, 1.3; m0 = 8.0)
    σs = x2σs([0.1, 0.2], ms; k = 2)

    _cosθ12 = cosθ12(σs, ms^2)
    mass_angles_cascade = [(ms.m0, 0.0, 1.0), (sqrt(σs.σ3), 0, _cosθ12)]
    angles, _σs = angles_invariants(mass_angles_cascade, ms; k = 3)
    @test all(collect(_σs) .≈ collect(σs))

    _cosθ23 = cosθ23(σs, ms^2)
    mass_angles_cascade = [(ms.m0, 0.0, 1.0), (sqrt(σs.σ1), 0, _cosθ23)]
    angles, _σs = angles_invariants(mass_angles_cascade, ms; k = 1)
    @test all(collect(_σs) .≈ collect(σs))

    _cosθ31 = cosθ31(σs, ms^2)
    mass_angles_cascade = [(ms.m0, 0.0, 1.0), (sqrt(σs.σ2), 0, _cosθ31)]
    angles, _σs = angles_invariants(mass_angles_cascade, ms; k = 2)
    @test all(collect(_σs) .≈ collect(σs))
end




model = let
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )
    dpp = randomPoint(tbs)
    #
    two_j = 2
    ch1 = DecayChain(;
        k = 1,
        two_j,
        Xlineshape = BreitWigner(1.1, 0.1),
        Hij = RecouplingLS((two_j, 0)) |> Vertex,
        HRk = RecouplingLS((two_j, two_j)) |> Vertex,
        tbs,
    )
    ch2 = DecayChain(ch1; k = 2)
    ch3 = DecayChain(ch1; k = 3)
    #
    ThreeBodyDecay("K892" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)]; reference_k = 2)
end

@testset "Order and values in validation_fields" begin
    ms = masses(model)
    σs = x2σs([0.2, 0.3], ms; k = 1)
    #
    model_name, point_name = "my_model", "my_point"
    _fields1 = validation_fields(model, σs; k = 2, point_name, model_name)
    _fields2 = validation_fields(
        model,
        (m = sqrt(σs[2]), cosθ = cosθ31(σs, ms^2));
        k = 2,
        point_name,
        model_name,
    )
    #
    @test _fields1[2]["value"] ≈ _fields2[2]["value"]
    @test _fields1[2]["point"] == _fields2[2]["point"] == point_name
    @test _fields1[2]["distribution"] == _fields2[2]["distribution"] == model_name
    d1 = array2dict(_fields1[1]["parameters"]; key = "name", apply = x -> x["value"])
    d2 = array2dict(_fields2[1]["parameters"]; key = "name", apply = x -> x["value"])
    #
    @test collect(keys(d1)) ==
          collect(keys(d2)) ==
          ["m_31_2", "phi_31_2", "cos_theta_31_2", "m_31", "phi_31", "cos_theta_31"]
    @test all(values(d1) .≈ values(d2))

    section = validation_section(
        model,
        [σs, σs];
        k = 1,
        point_names = point_name .* string(1:2),
        model_name,
    )
    @test array2dict(section[:parameter_points][1]["parameters"]; key = "name") |>
          keys |>
          collect ==
          ["m_23_1", "phi_23_1", "cos_theta_23_1", "m_23", "phi_23", "cos_theta_23"]
end
