"""A distribution bound to a model and its indexed phase-space domain."""
struct HadronicUnpolarizedIntensity{M,D}
    model::M
    domain::D
    function HadronicUnpolarizedIntensity(model::ThreeBodyDecay, domain::PhaseSpaceDomain)
        masses(model) == domain.masses ||
            throw(ArgumentError("model and domain masses differ"))
        new{typeof(model),typeof(domain)}(model, domain)
    end
end

function (dist::HadronicUnpolarizedIntensity)(point::AbstractDict)
    vectors = four_vectors(point, dist.domain)
    angles, sigmas = angles_invariants(vectors, dist.model)
    return unpolarized_intensity(dist.model, angles, sigmas)
end

function dict2instance(
    ::Type{HadronicUnpolarizedIntensity},
    dict;
    workspace,
    domain::PhaseSpaceDomain,
)
    dict["type"] == "HadronicUnpolarizedIntensity" ||
        throw(ArgumentError("unsupported distribution type $(dict["type"])"))
    haskey(dict, "variables") && throw(
        ArgumentError(
            "angle variables are obsolete; supply a four-vector phase-space point",
        ),
    )
    model = dict2instance(ThreeBodyDecay, dict["decay_description"]; workspace, domain)
    return HadronicUnpolarizedIntensity(model, domain)
end

function named_entry(entries, name, label)
    matches = filter(entry -> entry["name"] == name, entries)
    length(matches) == 1 || throw(ArgumentError("expected one $label named $name"))
    return only(matches)
end

const function_types = Dict(
    "BreitWigner" => BreitWigner,
    "MultichannelBreitWigner" => MultichannelBreitWigner,
    "Polynomial" => Polynomial,
    "BlattWeisskopf" => BlattWeisskopf,
    "MomentumPower" => MomentumPower,
    "generic_function" => generic_function,
)

function function_workspace(document; custom_functions = Dict())
    workspace = Dict{String,Any}()
    for fn in document["functions"]
        name = fn["name"]
        haskey(workspace, name) && throw(ArgumentError("duplicate function name $name"))
        if haskey(custom_functions, name)
            workspace[name] = custom_functions[name]
        else
            type = fn["type"]
            haskey(function_types, type) || throw(
                ArgumentError(
                    "unsupported function type $type; provide custom_functions[$name]",
                ),
            )
            workspace[name] = dict2instance(function_types[type], fn)
        end
    end
    return workspace
end

"""
    read_model(document; distribution=nothing, domain=nothing, custom_functions=Dict())

Read one distribution using masses from its phase-space domain. Select a name
explicitly when a document contains multiple distributions or domains. A checksum
uses the domain named by its point instead. No first-entry fallback is applied.
"""
function read_model(
    document::AbstractDict;
    distribution = nothing,
    domain = nothing,
    custom_functions = Dict(),
)
    distributions = document["distributions"]
    selected = if distribution === nothing
        length(distributions) == 1 ||
            throw(ArgumentError("select a distribution explicitly"))
        only(distributions)
    else
        named_entry(distributions, distribution, "distribution")
    end
    selected_domain = phase_space_domain(document; domain)
    return dict2instance(
        HadronicUnpolarizedIntensity,
        selected;
        workspace = function_workspace(document; custom_functions),
        domain = selected_domain,
    )
end

read_model(path::AbstractString; kw...) = read_model(JSON.parsefile(path); kw...)
