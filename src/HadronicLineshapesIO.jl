# deserialization

"""
    NamedArgFunc(f, variable_names)

Wrap a scalar function with its serialized argument name. Call the wrapper with
an `AbstractDict` containing that name; other dictionary entries are ignored.
Exactly one argument name is required.
"""
struct NamedArgFunc{M,D}
    f::M
    variable_names::D
    function NamedArgFunc(f::M, variable_names::D) where {M,D}
        length(variable_names) == 1 ||
            throw(ArgumentError("NamedArgFunc requires exactly one argument name"))
        new{M,D}(f, variable_names)
    end
end

NamedArgFunc(; f, variable_names) = NamedArgFunc(f, variable_names)

function (obj::NamedArgFunc)(dict::AbstractDict)
    @unpack variable_names = obj
    variable_value = dict[first(variable_names)]
    obj.f(variable_value)
end

# ## Functions of s

function dict2instance(::Type{BreitWigner}, dict::AbstractDict)
    @unpack mass, width, ma, mb, l, d, x = dict
    bw = BreitWigner(mass, width, ma, mb, l, d)
    variables = [x]
    NamedArgFunc(bw, variables)
end

function dict2instance(::Type{MultichannelBreitWigner}, dict::AbstractDict)
    @unpack mass, channels, x = dict
    # convert dict channels into NamedTuple
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        (; gsq, ma, mb, l, d)
    end
    bw = MultichannelBreitWigner(mass, _channels)
    variables = [x]
    return NamedArgFunc(bw, variables)
end

function dict2instance(::Type{Polynomial}, dict::AbstractDict)
    @unpack coefficients, x = dict
    P = Polynomial(coefficients, Symbol(x))
    arithmetic_P = WrapFlexFunction(P)
    variables = [x]
    NamedArgFunc(arithmetic_P, variables)
end

# ## Functions of p

function dict2instance(::Type{BlattWeisskopf}, dict)
    @unpack radius, l = dict
    return BlattWeisskopf{l}(radius)
end

function dict2instance(::Type{MomentumPower}, dict)
    @unpack l = dict
    return MomentumPower{l}()
end



# serialization

function serializeToDict(obj::NamedArgFunc{BreitWigner})
    type = "BreitWigner"
    @unpack f, variable_names = obj
    @unpack ma, mb, l, d = f
    x = first(variable_names)
    dict = LittleDict{String,Any}(
        "type" => type,
        "mass" => f.m,
        "width" => f.Γ,
        "ma" => ma,
        "mb" => mb,
        "l" => l,
        "d" => d,
        "x" => x,
    )
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(obj::NamedArgFunc{<:MultichannelBreitWigner})
    type = "MultichannelBreitWigner"
    @unpack f, variable_names = obj
    @unpack m, channels = f
    x = first(variable_names)
    _channels = map(channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        LittleDict{String,Any}("gsq" => gsq, "ma" => ma, "mb" => mb, "l" => l, "d" => d)
    end
    appendix = Dict()
    dict = LittleDict{String,Any}(
        "type" => type,
        "mass" => m,
        "channels" => _channels,
        "x" => x,
    )
    return (dict, appendix)
end


function serializeToDict(x::BlattWeisskopf)
    l = orbital_momentum(x)
    radius = x.d
    type = "BlattWeisskopf"
    dict = LittleDict{String,Any}("type" => type, "l" => l, "radius" => radius)
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(x::MomentumPower)
    l = orbital_momentum(x)
    type = "MomentumPower"
    dict = LittleDict{String,Any}("type" => type, "l" => l)
    appendix = Dict()
    return (dict, appendix)
end
function serializeToDict(H::NoFormFactor)
    appendix = Dict()
    ("", appendix)
end

function serializeToDict(obj::NamedArgFunc{<:WrapFlexFunction{<:Polynomial}})
    polynomial = obj.f.F
    dict = LittleDict{String,Any}(
        "type" => "Polynomial",
        "coefficients" => collect(coeffs(polynomial)),
        "x" => first(obj.variable_names),
    )
    return dict, Dict()
end
