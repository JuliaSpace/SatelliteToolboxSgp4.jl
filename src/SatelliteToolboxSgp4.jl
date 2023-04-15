module SatelliteToolboxSgp4

using Dates
using Reexport
using StaticArrays
using SatelliteToolboxBase
using SnoopPrecompile

@reexport using SatelliteToolboxTle

############################################################################################
#                                         Types
############################################################################################

include("types.jl")

############################################################################################
#                                        Constants
############################################################################################

# Julian Day related to the epoch 1900-01-01T12:00:00.000.
const _JD_1900 = DateTime(1900, 1, 1, 12, 0, 0) |> datetime2julian

############################################################################################
#                                         Includes
############################################################################################

include("sgp4_model.jl")
include("helpers.jl")

include("precompile.jl")

end
