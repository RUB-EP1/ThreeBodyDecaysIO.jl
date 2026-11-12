function serializeToDict(H::Vertex)
    type = "ls"
    @unpack h, ff = H
    ff_name, ff_appendix = serializeToDict(ff)
    h_dict, _ = serializeToDict(h) # helicity coupling does not produce an appendix
    H_dict = LittleDict{String,Any}("type" => type, "formfactor" => ff_name)
    merge!(H_dict, h_dict)
    (H_dict, ff_appendix)
end

function serializeToDict(H::RecouplingLS)
    type = "ls"
    l, s = H.two_ls .|> d2
    H_dict = LittleDict{String,Any}("type" => type, "l" => l, "s" => s)
    (H_dict, Dict())
end

function serializeToDict(H::ParityRecoupling)
    type = "parity"
    helicities = [H.two_λa, H.two_λb] .|> d2
    parity_factor = H.ηηηphaseisplus ? '+' : '-'
    H_dict = LittleDict{String,Any}(
        "type" => type,
        "helicities" => helicities,
        "parity_factor" => parity_factor,
    )
    (H_dict, Dict())
end

function serializeToDict(H::NoRecoupling)
    type = "helicity"
    helicities = [H.two_λa, H.two_λb] .|> d2
    H_dict = LittleDict{String,Any}("type" => type, "helicities" => helicities)
    (H_dict, Dict())
end


"""
    serializeToDict(chain::AbstractDecayChain;
        name::AbstractString="my_decay_chain")

Writes a `DecayChain` model to a dictionary including `propagators`, and `vertices`, and `topology`.
"""
function serializeToDict(chain::AbstractDecayChain; name::AbstractString = "my_decay_chain")
    @unpack k = chain
    i, j = ij_from_k(k)
    #
    appendix = Dict()

    # propagator
    propagator_name = "propagator_$(k)_$(hash(chain.Xlineshape))"
    named_func = NamedArgFunc(chain.Xlineshape, ["s"])
    merge!(appendix, Dict(propagator_name => named_func))
    propagator = LittleDict(
        "spin" => chain.two_j |> d2,
        "parametrization" => propagator_name,
        "node" => [i, j],
    )
    propagators = [propagator]
    #
    H1, a = serializeToDict(chain.HRk)
    push!(H1, "node" => [[i, j], k])
    merge!(appendix, a)
    #
    H2, a = serializeToDict(chain.Hij)
    push!(H2, "node" => [i, j])
    merge!(appendix, a)
    #
    vertices = [H1, H2]
    #
    topology = [[i, j], k]
    chain_dict = LittleDict{String,Any}(
        "vertices" => vertices,
        "propagators" => propagators,
        "topology" => topology,
        "name" => name,
    )
    (chain_dict, appendix)
end

function serializeToDict(
    tbs::ThreeBodySystem;
    particle_labels::NTuple{4,String} = ("A", "B", "C", "X"),
)
    #
    @unpack ms, two_js = tbs
    system_dict = LittleDict(
        "initial_state" => LittleDict(
            "name" => particle_labels[4],
            "mass" => ms[4],
            "index" => 0,
            "spin" => d2(two_js[4]),
        ),
        "final_state" => [
            LittleDict(
                "name" => particle_labels[1],
                "mass" => ms[1],
                "index" => 1,
                "spin" => d2(two_js[1]),
            ),
            LittleDict(
                "name" => particle_labels[2],
                "mass" => ms[2],
                "index" => 2,
                "spin" => d2(two_js[2]),
            ),
            LittleDict(
                "name" => particle_labels[3],
                "mass" => ms[3],
                "index" => 3,
                "spin" => d2(two_js[3]),
            ),
        ],
    )
    system_dict, Dict()
end


"""
serializeToDict::ThreeBodyDecay;
        particle_labels::NTuple{4,String}=("A", "B", "C", "X"))

Writes a `ThreeBodyDecay` model to a dictionary.
The argument `particle_labels` is passed to the kinematics serialization function.

## Arguments

- `model::ThreeBodyDecay`: The model to serialize.
- `particle_labels`: a tuple with labels of the particles in the decay.
```
"""
function serializeToDict(
    model::ThreeBodyDecay;
    particle_labels::NTuple{4,String} = ("A", "B", "C", "X"),
)
    #
    @unpack chains, names, couplings = model

    appendix = Dict()
    _kinematics, a = serializeToDict(model.chains[1].tbs; particle_labels)
    merge!(appendix, a)
    #
    _chains = map(zip(chains, names, couplings)) do (chain, name, coupling)
        dict, a = serializeToDict(chain; name)
        merge!(appendix, a)
        push!(dict, "weight" => replace(string(coupling), "im" => "i"))
    end
    #
    i, j, k = ijk(ThreeBodyDecays.reference_k(model))
    #
    model_dict = LittleDict{String,Any}(
        "kinematics" => _kinematics,
        "reference_topology" => [[i, j], k],
        "chains" => _chains,
    )
    model_dict, appendix
end


function add_hs3_fields(decay_description, appendix, model_name = "my_amplitude_model")

    k = topology2k(decay_description["reference_topology"])
    variable_groups = variablesToDict(k)
    variable_names = vcat(map(variable_groups) do v
        v["mass_phi_costheta"]
    end...)

    dict = OrderedDict(
        "distributions" => [
            OrderedDict(
                "type" => "HadronicUnpolarizedIntensity",
                "name" => model_name,
                "decay_description" => decay_description,
                "variables" => variable_groups,
            ),
        ],
        "functions" => [
            begin
                if v isa NamedArgFunc
                    v_dict, _ = serializeToDict(v)
                    v_dict["name"] = k
                    v_dict
                else
                    v["name"] = k
                    v
                end
            end for (k, v) in appendix
        ],
        "domains" => [
            OrderedDict(
                "name" => "default",
                "type" => "product_domain",
                "axes" => [
                    OrderedDict("name" => name, "min" => -1, "max" => 1) for
                    name in variable_names
                ],
            ),
        ],
        "misc" => Dict(),
        "parameter_points" => Dict(),
    )
    return dict
end


function variablesToDict(k::Int)
    i, j = ij_from_k(k)
    ij_str = "$(i)$(j)"
    ij_k_str = "$(i)$(j)_$(k)"
    return [
        Dict(
            "node" => [i, j],
            "mass_phi_costheta" => ["m_" * ij_str, "phi_" * ij_str, "cos_theta_" * ij_str],
        ),
        Dict(
            "node" => [[i, j], k],
            "mass_phi_costheta" =>
                ["m_" * ij_k_str, "phi_" * ij_k_str, "cos_theta_" * ij_k_str],
        ),
    ]
end
