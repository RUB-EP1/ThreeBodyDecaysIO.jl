"""
    validate(document; rtol=1e-10, atol=1e-12, custom_functions=Dict())

Evaluate each distribution checksum at its indexed four-vector point and selected
domain. Errors are recorded per checksum as `passed=false` and an `error` string;
one bad checksum does not stop subsequent checks. Scalar lineshape checksums are
not supported by the phase-space format.
"""
function validate(
    document::AbstractDict;
    rtol = 1e-10,
    atol = 1e-12,
    custom_functions = Dict(),
)
    checks = get(get(document, "misc", Dict()), "amplitude_model_checksums", [])
    isempty(checks) && throw(ArgumentError("document has no amplitude model checksums"))
    return map(checks) do check
        point_name = distribution = expected = nothing
        try
            check isa AbstractDict || throw(ArgumentError("checksum must be an object"))
            point_name = get(check, "point", nothing)
            distribution = get(check, "distribution", nothing)
            expected = get(check, "value", nothing)
            expected_value = finite_real(expected, "checksum value")
            point = named_entry(document["parameter_points"], point_name, "point")
            dist = read_model(
                document;
                distribution,
                domain = point["domain"],
                custom_functions,
            )
            value = dist(point)
            isfinite(value) ||
                throw(ArgumentError("distribution returned a nonfinite value"))
            return (;
                point = point_name,
                distribution,
                value,
                expected = expected_value,
                passed = isapprox(value, expected_value; rtol, atol),
                error = nothing,
            )
        catch err
            err isa InterruptException && rethrow()
            return (;
                point = point_name,
                distribution,
                value = nothing,
                expected,
                passed = false,
                error = sprint(showerror, err),
            )
        end
    end
end

validate(path::AbstractString; kw...) = validate(JSON.parsefile(path); kw...)

function validation_fields(
    model::ThreeBodyDecay,
    vectors;
    point_name,
    model_name,
    domain_name = "default",
)
    domain = PhaseSpaceDomain(domain_name, masses(model))
    point = four_vector_point(vectors, domain; name = point_name)
    angles, sigmas = angles_invariants(vectors, model)
    value = unpolarized_intensity(model, angles, sigmas)
    amplitude_model_checksums =
        LittleDict("distribution" => model_name, "point" => point_name, "value" => value)
    return (; point, amplitude_model_checksums)
end

function validation_section(model, events; point_names, model_name, domain_name = "default")
    length(events) == length(point_names) ||
        throw(ArgumentError("provide one name per validation event"))
    length(unique(point_names)) == length(point_names) ||
        throw(ArgumentError("validation point names must be unique"))
    fields = [
        validation_fields(model, event; point_name, model_name, domain_name) for
        (point_name, event) in zip(point_names, events)
    ]
    return LittleDict(
        "misc" => Dict(
            "amplitude_model_checksums" =>
                getfield.(fields, :amplitude_model_checksums),
        ),
        "parameter_points" => getfield.(fields, :point),
    )
end
