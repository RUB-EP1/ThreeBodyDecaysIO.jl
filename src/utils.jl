


function update2values(x, ref)
    if x isa String && haskey(ref, x)
        return ref[x]
    end
    if x isa Dict
        _d = Dict{String,Any}()
        for (k, v) in x
            _d[k] = update2values(v, ref)
        end
        return _d
    end
    if x isa AbstractArray
        return update2values.(x, Ref(ref))
    end
    return x
end

function string2complex(value::AbstractString)
    # Base's numeric parser accepts complex literals without evaluating Julia code.
    result = tryparse(ComplexF64, replace(value, r"i$" => "im"))
    result === nothing && throw(ArgumentError("invalid complex numeric literal: $value"))
    isfinite(real(result)) && isfinite(imag(result)) ||
        throw(ArgumentError("complex value must be finite"))
    return result
end

function topology2k(topology::AbstractArray)
    indices = flatten_topology(topology)
    length(indices) != 3 &&
        error("Topology with more that three particles is not implemented")
    !(topology ∈ ([[1, 2], 3], [[2, 3], 1], [[3, 1], 2])) && error(
        "Only regular topologies [[1, 2], 3], [[2, 3], 1], [[3, 1], 2] are implemented",
    )
    #
    return topology[2]
end

function reorder(last_index)
    k = last_index
    i, j = ij_from_k(k)
    t -> invpermute!(collect(t), [i, j, k]) |> Tuple
end

"""
    array2dict(a::AbstractArray; key, apply = identity)

Convert an array of objects to a dictionary by extracting a key from each object.

This function transforms structured data from JSON arrays (which preserve order)
into Julia dictionaries (which enable fast key-based lookups). For each object in the array,
it extracts the specified key to use as the dictionary key, and applies a transformation
function to the remaining data.

# Arguments
- `a`: Array of objects (typically from JSON parsing)
- `key`: The field name to extract and use as the dictionary key
- `apply`: Function to transform the remaining data (default: identity)

# Example
```julia
data = [Dict("node" => [1,2], "vars" => ["a","b"]),
        Dict("node" => [3,4], "vars" => ["c","d"])]
result = array2dict(data; key="node", apply=x->x["vars"])
# Returns: LittleDict([1,2] => ["a","b"], [3,4] => ["c","d"])
```
"""
function array2dict(a::AbstractArray; key, apply = identity)
    map(a) do p
        _key = p[key]
        # Extract remaining data without the key
        remaining_data = Dict(k => v for (k, v) in p if k != key)
        _key => apply(remaining_data)
    end |> LittleDict
end

"""
    flatten_topology(topology)

Converts the topology structure into a flat array of indices.

## Examples:
```juliadoc
julia> flatten_topology([[1, 2], 3])
[1, 2, 3]
```

Can be used for checking the validity of the topology.
```julia
    @assert flatten_topology(reference_topology) |> sort == [1, 2, 3] "Error: allowed indices are only 1,2,3"
```
"""
flatten_topology(topology) =
    topology isa Array ? vcat(flatten_topology.(topology)...) : topology


function label_diff(diff; levels = [1e-2, 1e-12])
    _diff = abs(diff)
    _diff < levels[2] && return '🟢'
    _diff < levels[1] && return '🟡'
    return '🔴'
end
