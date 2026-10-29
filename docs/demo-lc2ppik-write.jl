## This notebook take models from Lc2ppiKSemileptonicModelLHCb.jl
# and serialize them into a dictionary that is saved to JSON file

using Lc2ppiKSemileptonicModelLHCb
using ThreeBodyDecaysIO
using ThreeBodyDecays
using HadronicLineshapes
using Parameters
using JSON

model = published_model("Default amplitude model")

function lineshape_parser(Xlineshape::BreitWignerMinL)
    appendix = Dict()
    Xlineshape_str = Xlineshape.name * "_BW"
    scattering, a = Xlineshape_str,
    Dict(
        Xlineshape_str => Dict(
            :type => "BreitWigner",
            :mass => Xlineshape.pars.m,
            :width => Xlineshape.pars.Γ,
            :l => Xlineshape.l,
            :ma => Xlineshape.m1,
            :mb => Xlineshape.m2,
        ),
    )
    merge!(appendix, a)
    FF_production = "BlattWeisskopf(resonance)"
    FF_decay = "BlattWeisskopf(b-decay)"
    # hard coded
    dR, dΛc = 1.5, 5.0 # /GeV
    a = Dict(
        "BlattWeisskopf(resonance)" =>
            Dict(:type => "BlattWeisskopf", :radius => dR, :l => Xlineshape.l),
        "BlattWeisskopf(b-decay)" =>
            Dict(:type => "BlattWeisskopf", :radius => dΛc, :l => 0),
    )
    merge!(appendix, a)
    (; scattering, FF_production, FF_decay), appendix
end

function lineshape_parser(Xlineshape::Flatte1405)
    appendix = Dict()
    Xlineshape_str = Xlineshape.name * "_Flatte"
    scattering, a = Xlineshape_str,
    Dict(
        Xlineshape_str => Dict(
            :type => typeof(Xlineshape),
            :mass => Xlineshape.pars.m,
            :width => Xlineshape.pars.Γ,
        ),
    )
    merge!(appendix, a)
    FF_production = ""
    FF_decay = ""
    (; scattering, FF_production, FF_decay), appendix
end

function lineshape_parser(Xlineshape::BuggBreitWignerMinL)
    appendix = Dict()
    Xlineshape_str = Xlineshape.name * "_BuggBW"
    scattering, a = Xlineshape_str,
    Dict(
        Xlineshape_str => Dict(
            :type => typeof(Xlineshape),
            :mass => Xlineshape.pars.m,
            :width => Xlineshape.pars.Γ,
            :gamma => Xlineshape.pars.γ,
        ),
    )
    merge!(appendix, a)
    FF_production = ""
    FF_decay = ""
    (; scattering, FF_production, FF_decay), appendix
end

lineshape_parser(model.chains[6].Xlineshape) isa Tuple
lineshape_parser(model.chains[1].Xlineshape) isa Tuple
lineshape_parser(model.chains[end].Xlineshape) isa Tuple

dict, appendix = serializeToDict(model; lineshape_parser)
dict[:kinematics][:names] = ["p", "pi", "K", "Lc"]
dict[:functions] = [(v[:name] = k; v) for (k, v) in appendix]
#
test_file_name = "lc2ppi-test.json"
open(test_file_name, "w") do io
    JSON.print(io, dict, 4)
end

dict[:kinematics]

## Reading and building the model

json_content = open(test_file_name) do io
    JSON.parse(io)
end
input = copy(json_content)

workspace = let
    _workspace = Dict{String,Any}()
    @unpack functions = input
    for fn in functions
        @show fn
        @unpack name, type = fn
        @info "   🎈$name"
        instance_type = eval(Symbol(type))
        _workspace[name] = dict2instance(instance_type, fn)
    end
    _workspace
end

# updated_input = update2values(input, json_content["appendix"])
tbs = dict2instance(ThreeBodySystem, input["kinematics"])
cdn = dict2instance(DecayChain, input["chains"][1]; tbs, workspace)

map(updated_input["chains"]) do chain_dict
    chain_dict["propagators"][1]["parametrization"]["type"]
end
rm(test_file_name)
