@testset "Domain selection and checksum failure isolation" begin
    @test_throws ArgumentError read_model(
        JSON.parsefile(joinpath(@__DIR__, "lc2ppi-lhcb-test.json")),
    )
    document = add_hs3_fields(model)
    ms = masses(model)
    vectors = map(
        FourVectors.FourVector,
        aligned_four_vectors(randomPoint(ms), ms; k = reference_k(model)),
    )
    merge!(
        document,
        validation_section(
            model,
            [vectors];
            point_names = ["point"],
            model_name = "my_amplitude_model",
        ),
    )
    @test only(validate(document)).passed
    parsed = read_model(document)
    point = only(document["parameter_points"])
    @test parsed(point) > 0
    @test_throws MethodError parsed(point; k = 1)
    permuted = deepcopy(document)
    reverse!(permuted["domains"][1]["final_state"])
    reverse!(permuted["distributions"][1]["decay_description"]["kinematics"]["final_state"])
    reverse!(permuted["parameter_points"][1]["final_state"])
    @test only(validate(permuted)).passed
    broken = deepcopy(document)
    broken["parameter_points"][1]["final_state"][1]["index"] = 2
    @test !only(validate(broken)).passed
    badchecks = deepcopy(document)
    badchecks["misc"] = Dict{String,Any}(badchecks["misc"])
    badchecks["misc"]["amplitude_model_checksums"] =
        Any[nothing, only(document["misc"]["amplitude_model_checksums"])]
    @test getproperty.(validate(badchecks), :passed) == [false, true]
    @test_throws ArgumentError validate(add_hs3_fields(model))
    two_domains = deepcopy(document)
    extra = deepcopy(first(two_domains["domains"]))
    extra["name"] = "other"
    push!(two_domains["domains"], extra)
    @test_throws ArgumentError read_model(two_domains)
    @test read_model(two_domains; domain = "default")(point) ≈ parsed(point)
    @test only(validate(two_domains)).passed
    @test_throws ArgumentError read_model(two_domains; domain = "missing")
    wrong = deepcopy(point)
    wrong["domain"] = "other"
    @test_throws ArgumentError parsed(wrong)
    legacy = deepcopy(document)
    legacy["distributions"][1]["decay_description"]["kinematics"]["initial_state"]["mass"] =
        ms.m0
    @test_throws ArgumentError read_model(legacy)
    unsupported = deepcopy(document)
    unsupported["functions"][1]["type"] = "exit"
    @test_throws ArgumentError read_model(unsupported)
end

@testset "Numeric and analytic input restrictions" begin
    @test string2complex("1.5-2.3i") == 1.5-2.3im
    for value in ("exit()", "1; error(1)", "NaN", "Inf")
        @test_throws ArgumentError string2complex(value)
    end
    for expression in
        ("run(x)", "Base.sin(x)", "x[1]", "begin x; exit() end", "x=1", "@time x")
        @test_throws ArgumentError ThreeBodyDecaysIO.parse_into_function(expression)
    end
end
