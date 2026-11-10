module ThreeBodyDecaysIO

using ThreeBodyDecays
import FourVectors
using ThreeBodyDecays: RecouplingLS, ParityRecoupling, NoRecoupling
using ThreeBodyDecays: AbstractDecayChain
using JSON
using OrderedCollections
using Parameters
using DataFrames
using HadronicLineshapes
using MacroTools

import Polynomials: Polynomial, coeffs
export Polynomial

using Random # for name generation
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)


export serializeToDict
export add_hs3_fields

export NamedArgFunc
include("HadronicLineshapesIO.jl")

export generic_function
export expression_argument, expression_arguments
include("generic_function.jl")

export PhaseSpaceDomain, phase_space_domain, four_vectors, four_vector_point
include("phase_space.jl")
include("writer.jl")

export read_model, validate
export HadronicUnpolarizedIntensity
include("HadronicUnpolarizedIntensity.jl")

export dict2instance
include("reader.jl")

export reorder
export array2dict
export string2complex
export update2values
export topology2k, flatten_topology
export label_diff
include("utils.jl")

export angles_invariants
export validation_fields
export validation_section
include("validation.jl")

end # module ThreeBodyDecaysIO
