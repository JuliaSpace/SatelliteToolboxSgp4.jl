module SatelliteToolboxSgp4

using Dates
using StaticArrays
using SatelliteToolboxTle
using SnoopPrecompile

############################################################################################
#                                         Types
############################################################################################

include("types.jl")

############################################################################################
#                                         Includes
############################################################################################

include("gmst.jl")
include("sgp4_model.jl")
include("helpers.jl")

include("precompile.jl")

end
