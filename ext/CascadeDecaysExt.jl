module CascadeDecaysExt

import CascadeDecays as CD
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.OrderedCollections: LittleDict
using ThreeBodyDecaysIO.ThreeBodyDecays: d2
import ThreeBodyDecaysIO: serializeToDict

json_topology(x::Integer) = x
json_topology(x::Tuple) = Any[json_topology(child) for child in x]
topology_dict(topology) = json_topology(
    CD.vertex_address(topology, CD.consumed_by(topology, CD.root_line_ind(topology))),
)

function serializeToDict(chain::CD.DecayChain; name::AbstractString = "my_decay_chain")
    topology = chain.topology
    appendix = Dict{String,Any}()
    propagators = map(
        collect(zip(CD.propagating_line_inds(chain), chain.propagators)),
    ) do (line, shape)
        vertex = CD.consumed_by(topology, line)
        node = json_topology(CD.vertex_address(topology, vertex))
        # Use the same key as the three-body writer for equivalent isobar chains.
        prefix =
            CD.nfinal(chain) == 3 ?
            string(only(setdiff(1:3, CD.final_descendants(chain, line)))) :
            join(CD.final_descendants(chain, line), "_")
        key = "propagator_$(prefix)_$(hash(shape))"
        appendix[key] = NamedArgFunc(shape, ["s"])
        LittleDict{String,Any}(
            "spin" => d2(CD.line_two_js(chain)[line]),
            "parametrization" => key,
            "node" => node,
        )
    end
    vertices = map(collect(enumerate(chain.vertices))) do (index, vertex)
        dict, extra = serializeToDict(vertex)
        dict["node"] = json_topology(CD.vertex_address(topology, index))
        merge!(appendix, extra)
        dict
    end
    return LittleDict{String,Any}(
        "vertices" => vertices,
        "propagators" => propagators,
        "topology" => topology_dict(topology),
        "name" => name,
    ),
    appendix
end

"""
    serializeToDict(model::CascadeDecays.CascadeDecay; masses, particle_labels)

Serialize a cascade using external masses ordered as final states followed by
its parent. CascadeDecay stores spins and topology, so masses must be supplied.
The returned decay description and function appendix follow the three-body
writer's format. Cascade deserialization is not implemented.
"""
function serializeToDict(model::CD.CascadeDecay; masses, particle_labels = nothing)
    chain = first(model.chains)
    n = CD.nfinal(chain)
    length(masses) == n + 1 ||
        throw(ArgumentError("supply $n final masses followed by the parent mass"))
    all(m -> m isa Real && isfinite(m) && m >= 0, masses) ||
        throw(ArgumentError("masses must be finite nonnegative real numbers"))
    labels =
        isnothing(particle_labels) ?
        (n == 3 ? ("A", "B", "C", "X") : (ntuple(i -> "particle_$i", n)..., "X")) :
        particle_labels
    length(labels) == n + 1 || throw(ArgumentError("particle_labels must match masses"))
    spins = CD.line_two_js(chain)
    initial = LittleDict(
        "name" => labels[end],
        "mass" => masses[end],
        "index" => 0,
        "spin" => d2(spins[CD.root_line_ind(chain)]),
    )
    finals = [
        LittleDict(
            "name" => labels[i],
            "mass" => masses[i],
            "index" => i,
            "spin" => d2(spins[line]),
        ) for (i, line) in enumerate(CD.final_line_inds(chain))
    ]
    appendix = Dict{String,Any}()
    chains = map(
        collect(zip(model.chains, model.names, model.couplings)),
    ) do (chain, name, coupling)
        dict, extra = serializeToDict(chain; name)
        merge!(appendix, extra)
        dict["weight"] = replace(string(coupling), "im" => "i")
        dict
    end
    return LittleDict{String,Any}(
        "kinematics" => LittleDict("initial_state" => initial, "final_state" => finals),
        "reference_topology" => topology_dict(CD.reference_topology(model)),
        "chains" => chains,
    ),
    appendix
end

function serializeToDict(obj::NamedArgFunc{<:CD.ConstantLineshape})
    value = obj.f.value
    dict = LittleDict{String,Any}(
        "type" => "Polynomial",
        "x" => first(obj.variable_names),
        "coefficients" => [value],
    )
    value isa Real ||
        throw(ArgumentError("constant lineshape serialization requires a real value"))
    return dict, Dict()
end

end
