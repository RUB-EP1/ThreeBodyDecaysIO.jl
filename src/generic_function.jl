
make_julia_complex(expr) = replace(string(expr), r"(?<![a-zA-Z])i(?![a-zA-Z])" => "im")

const operator_symbols = Set([
    :+,
    :-,
    :*,
    :/,
    :%,
    :^,
    :im,
    :i,
    :sin,
    :log,
    :cos,
    :tan,
    :asin,
    :acos,
    :atan,
    :sinh,
    :cosh,
    :tanh,
    :asinh,
    :acosh,
    :atanh,
    :exp,
    :sqrt,
    :abs,
    :real,
    :imag,
    :conj,
    :pow,
    :log10,
    :angle,
    :abs2,
    :sign,
    :round,
    :trunc,
    :floor,
    :ceil,
])

# Step 2: Function to traverse the syntax tree and collect symbols using MacroTools
function collect_nonregistered_symbols(expr)
    symbols = Set{Symbol}()
    MacroTools.prewalk(expr) do node
        if isa(node, Symbol) && !(node in operator_symbols)
            push!(symbols, node)
        end
        node
    end
    return symbols
end

expression_arguments(expr) = collect_nonregistered_symbols(expr)
function expression_argument(expr)
    args = expression_arguments(expr)
    length(args) != 1 && error("The expression has multiple arguments: $args")
    return first(args)
end

function checked_analytic_expression(node)
    node isa Number && isfinite(node) && return node
    node isa Symbol && return node
    node isa Expr && node.head == :call ||
        throw(ArgumentError("Only numeric analytic expressions are supported"))
    op = first(node.args)
    op isa Symbol && op in operator_symbols && !(op in (:i, :im)) ||
        throw(ArgumentError("Unsupported analytic operation: $op"))
    return Expr(
        :call,
        op == :pow ? :^ : op,
        map(checked_analytic_expression, node.args[2:end])...,
    )
end

function parse_into_function(body, pars = Dict("i" => 1im))
    expression = checked_analytic_expression(Meta.parse(body))
    parameters = Dict{Symbol,Any}(:i => 1im, :im => 1im)
    for (key, value) in pars
        symbol = Symbol(key)
        Base.isidentifier(string(key)) || throw(ArgumentError("Invalid parameter name"))
        symbol in operator_symbols &&
            !(symbol in (:i, :im)) &&
            throw(ArgumentError("Parameter shadows an analytic operation"))
        number = value isa AbstractString ? string2complex(value) : value
        number isa Number && isfinite(number) ||
            throw(ArgumentError("Analytic parameters must be finite numbers"))
        parameters[symbol] = number
    end
    variables = setdiff(expression_arguments(expression), keys(parameters))
    length(variables) == 1 ||
        throw(ArgumentError("Expected one independent variable: $variables"))
    arg = only(variables)
    substituted = MacroTools.prewalk(expression) do node
        node isa Symbol && haskey(parameters, node) ? parameters[node] : node
    end
    f = @RuntimeGeneratedFunction Expr(:->, arg, substituted)
    return f, arg
end

struct generic_function end

function dict2instance(::Type{generic_function}, dict)
    @unpack expression = dict
    parameters_dict = Dict{String,Any}(dict)
    pop!(parameters_dict, "expression")
    haskey(parameters_dict, "name") && pop!(parameters_dict, "name")
    haskey(parameters_dict, "type") && pop!(parameters_dict, "type")
    parameters_dict["i"] = 1im
    f, arg = parse_into_function(expression, parameters_dict)
    # arg = expression_argument(parsed_expression)
    return NamedArgFunc(WrapFlexFunction(f), [string(arg)])
end
