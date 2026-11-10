function dict2instance(
    ::Type{ThreeBodySystem},
    dict::AbstractDict;
    domain::PhaseSpaceDomain,
)
    initial =
        indexed_states([dict["initial_state"]], [0]; label = "initial-state kinematics")
    finals = indexed_states(dict["final_state"], 1:3; label = "final-state kinematics")
    any(haskey(state, "mass") for state in (initial[0], values(finals)...)) && throw(
        ArgumentError(
            "kinematics masses are obsolete; masses belong only in the phase_space domain",
        ),
    )
    two_js = ThreeBodySpins(
        (x2(finals[i]["spin"]) for i = 1:3)...;
        two_h0 = x2(initial[0]["spin"]),
    )
    return ThreeBodySystem(domain.masses, two_js)
end

function dict2instance(::Type{DecayChain}, dict; tbs, workspace = Dict())
    coupling = string2complex(dict["weight"])
    name = dict["name"]
    #
    @unpack vertices, propagators = dict
    @unpack topology = dict
    #
    k = topology2k(topology)
    i, j = ij_from_k(k)
    # spin of subchannel resonance
    resonance = first(propagators)
    @assert resonance["node"] == [i, j]
    #
    spin = resonance["spin"]
    two_j = spin |> x2
    two_js = tbs.two_js
    # build two vertices manually
    ind_ij = findfirst(v -> v["node"] == [i, j], vertices)
    vertex_ij = vertices[ind_ij]
    Hij = dict2instance(
        ThreeBodyDecays.Vertex,
        vertex_ij,
        (; two_j_fin = [two_js[i], two_js[j]], two_j_ini = two_j);
        workspace,
    )
    #
    ind_Rk = findfirst(v -> v["node"] == [[i, j], k], vertices)
    vertex_Rk = vertices[ind_Rk]
    HRk = dict2instance(
        ThreeBodyDecays.Vertex,
        vertex_Rk,
        (; two_j_fin = [two_j, two_js[k]], two_j_ini = two_js[4]);
        workspace,
    )

    # build lineshape
    scattering = resonance["parametrization"]
    scattering isa Dict && error(
        "Deserialization of lineshape directly from the chains is not implemented yet. Use `functions` field instead.",
    )
    !(scattering isa String) && error("The scattering should be a string. Got: $scattering")
    X = workspace[scattering].f
    chain = DecayChain(; k, two_j, Xlineshape = X, Hij, HRk, tbs)
    (; coupling, chain, name)
end


function dict2instance(
    ::Type{ThreeBodyDecay},
    decay_description;
    workspace,
    domain::PhaseSpaceDomain,
)
    @unpack kinematics, chains = decay_description
    tbs = dict2instance(ThreeBodySystem, kinematics; domain)
    df = dict2instance.(DecayChain, chains; tbs, workspace) |> DataFrame
    model = ThreeBodyDecay(
        Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(
            df.name .=> zip(df.coupling, df.chain),
        );
        reference_k = topology2k(decay_description["reference_topology"]),
    )
    return model
end

function dict2instance(::Type{ThreeBodyDecays.Vertex}, dict, properties; workspace = Dict())

    formfactor_name = dict["formfactor"]
    # Handle formfactor: if empty string, use NoFormFactor(), otherwise get from workspace
    if formfactor_name == ""
        formfactor = NoFormFactor()
    else
        formfactor = workspace[formfactor_name]
    end
    helicity = dict2instance(ThreeBodyDecays.Recoupling, dict, properties)

    return Vertex(helicity, formfactor)
end

function dict2instance(::Type{ThreeBodyDecays.Recoupling}, dict, properties)
    if dict["type"] == "ls"
        @unpack l, s = dict
        @unpack two_j_ini = properties
        return RecouplingLS(; two_ls = (l, s) .|> x2)
    end
    if dict["type"] == "helicity"
        @unpack helicities = dict
        two_λa = helicities[1] |> x2
        two_λb = helicities[2] |> x2
        return NoRecoupling(two_λa, two_λb)
    end
    if dict["type"] == "parity"
        @unpack helicities, parity_factor = dict
        two_λa = helicities[1] |> x2
        two_λb = helicities[2] |> x2
        return ParityRecoupling(two_λa, two_λb, parity_factor == "+")
    end
    error("Unknown type: $(dict["type"])")
end
