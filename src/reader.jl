function dict2instance(::Type{ThreeBodySystem}, dict::AbstractDict)
    initial_state = dict["initial_state"]
    final_states = dict["final_state"]

    # Parse initial state
    initial_index = initial_state["index"]
    initial_name = initial_state["name"]
    initial_spin = initial_state["spin"]
    initial_mass = initial_state["mass"]

    # Parse final states
    final_indices = [state["index"] for state in final_states]
    final_names = [state["name"] for state in final_states]
    final_spins = [state["spin"] for state in final_states]
    final_masses = [state["mass"] for state in final_states]

    # Create DataFrame
    df = DataFrame(
        indices = [initial_index; final_indices],
        names = [initial_name; final_names],
        spins = [initial_spin; final_spins],
        masses = [initial_mass; final_masses],
    )

    sort!(df, :indices)
    N = size(df, 1) - 1
    first(df.indices) != 0 && error("particle 0 should be indicating the overall system")
    !all(df.indices .== 0:N) &&
        error("Expected indices 0,1,...,$N for the reaction with N final particles")
    N != 3 && error("ThreeBodyDecays.jl only works with three-body systems only")

    ms = ThreeBodyMasses(df.masses[2:end]...; m0 = df.masses[1])
    two_spins = df.spins .|> x2
    two_js = ThreeBodySpins(two_spins[2:end]...; two_h0 = two_spins[1])

    tbs = ThreeBodySystem(ms, two_js)
    return tbs
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

function build_or_fetch(key, workspace)
    key isa Dict && return dict2object(key)
    return workspace[key]
end

function dict2instance(::Type{ThreeBodyDecay}, decay_description; workspace)
    @unpack kinematics, chains, reference_topology = decay_description
    reference_k = topology2k(reference_topology)
    tbs = dict2instance(ThreeBodySystem, kinematics)
    df = dict2instance.(DecayChain, chains; tbs, workspace) |> DataFrame
    model = ThreeBodyDecay(
        Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(
            df.name .=> zip(df.coupling, df.chain),
        );
        reference_k,
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
