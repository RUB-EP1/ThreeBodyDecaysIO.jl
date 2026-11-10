"""
    PhaseSpaceDomain(name, masses)

Three-body phase space with final-state masses indexed 1, 2, 3 and parent index 0.
The JSON domain is the sole source of particle masses.
"""
struct PhaseSpaceDomain{M}
    name::String
    masses::M
    function PhaseSpaceDomain(name::AbstractString, masses::ThreeBodyDecays.MassTuple)
        all(m -> isfinite(m) && m >= 0, masses) ||
            throw(ArgumentError("phase-space masses must be finite and nonnegative"))
        masses.m0 > sum(masses[i] for i = 1:3) ||
            throw(ArgumentError("parent mass must exceed the sum of final-state masses"))
        new{typeof(masses)}(String(name), masses)
    end
end

function indexed_states(states, expected; label)
    result = Dict{Int,Any}()
    for state in states
        index = state["index"]
        index isa Integer && !(index isa Bool) ||
            throw(ArgumentError("$label indices must be integers"))
        haskey(result, index) && throw(ArgumentError("duplicate $label index $index"))
        result[index] = state
    end
    Set(keys(result)) == Set(expected) ||
        throw(ArgumentError("$label indices must be $(collect(expected))"))
    return result
end

function finite_real(value, label)
    value isa Real && !(value isa Bool) && isfinite(value) ||
        throw(ArgumentError("$label must be a finite real number"))
    return Float64(value)
end

function dict2instance(::Type{PhaseSpaceDomain}, dict::AbstractDict)
    dict["type"] == "phase_space" || throw(ArgumentError("expected a phase_space domain"))
    initial = indexed_states([dict["initial_state"]], [0]; label = "initial-state domain")
    finals = indexed_states(dict["final_state"], 1:3; label = "final-state domain")
    m0 = finite_real(initial[0]["mass"], "parent mass")
    ms = ntuple(i -> finite_real(finals[i]["mass"], "mass $i"), 3)
    all(>=(0), ms) && m0 >= 0 || throw(ArgumentError("masses must be nonnegative"))
    m0 > sum(ms) || throw(ArgumentError("parent mass must exceed the final-state masses"))
    return PhaseSpaceDomain(dict["name"], ThreeBodyMasses(ms...; m0))
end

function phase_space_domain(document::AbstractDict; domain = nothing)
    candidates = [d for d in document["domains"] if d["type"] == "phase_space"]
    names = [d["name"] for d in candidates]
    length(unique(names)) == length(names) ||
        throw(ArgumentError("duplicate phase-space domain names"))
    if domain === nothing
        length(candidates) == 1 || throw(
            ArgumentError(
                "select a domain explicitly when there is not exactly one phase_space domain",
            ),
        )
        return dict2instance(PhaseSpaceDomain, only(candidates))
    end
    matches = filter(d -> d["name"] == domain, candidates)
    length(matches) == 1 || throw(ArgumentError("unknown phase_space domain $domain"))
    return dict2instance(PhaseSpaceDomain, only(matches))
end

function serializeToDict(domain::PhaseSpaceDomain)
    return LittleDict{String,Any}(
        "name" => domain.name,
        "type" => "phase_space",
        "initial_state" => LittleDict{String,Any}("index" => 0, "mass" => domain.masses.m0),
        "final_state" =>
            [LittleDict{String,Any}("index" => i, "mass" => domain.masses[i]) for i = 1:3],
    ),
    Dict()
end

function check_four_vectors(
    vectors,
    ms;
    atol_mass2 = 1e-6,
    atol_momentum = 1e-8,
    frame = :rest,
)
    length(vectors) == 3 ||
        throw(ArgumentError("expected exactly three final-state four-vectors"))
    frame in (:rest, :canonical) ||
        throw(ArgumentError("frame must be :rest or :canonical"))
    for (i, p) in enumerate(vectors)
        all(isfinite, p) && p.E > 0 ||
            throw(ArgumentError("four-vector $i must be finite with positive energy"))
        isapprox(FourVectors.mass2(p), ms[i]^2; atol = atol_mass2, rtol = 0) ||
            throw(ArgumentError("four-vector $i is off shell for its domain mass"))
    end
    parent = sum(vectors)
    isapprox(FourVectors.mass2(parent), ms.m0^2; atol = atol_mass2, rtol = 0) ||
        throw(ArgumentError("total invariant mass does not match the domain parent mass"))
    if frame == :rest
        maximum(abs, (parent.px, parent.py, parent.pz, parent.E - ms.m0)) <=
        atol_momentum ||
            throw(ArgumentError("phase-space point must be in the parent rest frame"))
    end
    return parent
end

"""
    four_vectors(point, domain; atol_mass2=1e-6, atol_momentum=1e-8)

Read PR111 `final_state = [{index, p = {E, px, py, pz}}, ...]` points.
Match by index, validate masses and rest-frame momentum conservation, and return
four-vectors in index order. Scalar and legacy `four_vectors` points are rejected.
"""
function four_vectors(
    point::AbstractDict,
    domain::PhaseSpaceDomain;
    atol_mass2 = 1e-6,
    atol_momentum = 1e-8,
)
    point["domain"] == domain.name ||
        throw(ArgumentError("point and selected domain do not match"))
    states = indexed_states(point["final_state"], 1:3; label = "momentum")
    vectors = ntuple(3) do i
        p = states[i]["p"]
        FourVectors.FourVector(
            finite_real(p["px"], "px"),
            finite_real(p["py"], "py"),
            finite_real(p["pz"], "pz");
            E = finite_real(p["E"], "E"),
        )
    end
    check_four_vectors(vectors, domain.masses; atol_mass2, atol_momentum)
    return vectors
end

function canonical_rest_vector(p, parent)
    spatial = FourVectors.spatial_magnitude(parent)
    spatial <= eps(Float64) * parent.E && return p
    phi = atan(parent.py, parent.px)
    theta = acos(clamp(parent.pz / spatial, -1.0, 1.0))
    gamma = parent.E / sqrt(FourVectors.mass2(parent))
    return p |>
           FourVectors.Rz(-phi) |>
           FourVectors.Ry(-theta) |>
           FourVectors.Bz(-gamma) |>
           FourVectors.Ry(theta) |>
           FourVectors.Rz(phi)
end

"""
    angles_invariants(vectors, masses; k, frame=:rest)

Recover `(α, cosβ, γ)` and the three pair invariants for the reference topology
`[[i,j],k]`. The isobar points along +z and daughter i defines the positive-x
half-plane before the active rotation Rz(α) Ry(β) Rz(γ).

The default requires a parent-rest-frame event. `frame=:canonical` first applies
a rotation-free boost to the parent rest frame; resulting polarization axes are
those of that canonical frame. The parent helicity frame is not inferred.
"""
function angles_invariants(
    vectors::Union{
        Tuple{Vararg{FourVectors.FourVector}},
        AbstractVector{<:FourVectors.FourVector},
    },
    ms::ThreeBodyDecays.MassTuple;
    k,
    frame = :rest,
    kw...,
)
    k isa Integer && !(k isa Bool) && k in 1:3 ||
        throw(ArgumentError("reference topology k must be 1, 2, or 3"))
    parent = check_four_vectors(vectors, ms; frame, kw...)
    rest =
        frame == :canonical ? map(p -> canonical_rest_vector(p, parent), vectors) : vectors
    i, j = ij_from_k(k)
    isobar = rest[i] + rest[j]
    radius = FourVectors.spatial_magnitude(isobar)
    radius > 1e-12 * ms.m0 ||
        throw(ArgumentError("reference isobar is at rest: orientation is undefined"))
    transverse = hypot(isobar.px, isobar.py)
    alpha = transverse <= 1e-12 * radius ? 0.0 : atan(isobar.py, isobar.px)
    cosbeta = clamp(isobar.pz / radius, -1.0, 1.0)
    daughter = rest[i] |> FourVectors.Rz(-alpha) |> FourVectors.Ry(-acos(cosbeta))
    gamma =
        hypot(daughter.px, daughter.py) <= 1e-12 * ms.m0 ? 0.0 :
        atan(daughter.py, daughter.px)
    sigmas = ntuple(3) do index
        a, b = ij_from_k(index)
        FourVectors.mass2(rest[a] + rest[b])
    end
    return (α = alpha, cosβ = cosbeta, γ = gamma), MandelstamTuple{Float64}(sigmas)
end

function four_vector_point(vectors, domain::PhaseSpaceDomain; name = "validation_point")
    check_four_vectors(vectors, domain.masses)
    return LittleDict{String,Any}(
        "name" => name,
        "domain" => domain.name,
        "final_state" => [
            LittleDict(
                "index" => i,
                "p" => LittleDict("E" => p.E, "px" => p.px, "py" => p.py, "pz" => p.pz),
            ) for (i, p) in enumerate(vectors)
        ],
    )
end

"""Recover orientation and invariants using the model’s stored reference topology."""
angles_invariants(vectors, model::ThreeBodyDecay) =
    angles_invariants(vectors, masses(model); k = reference_k(model))
