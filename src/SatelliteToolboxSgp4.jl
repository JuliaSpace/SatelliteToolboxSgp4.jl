module SatelliteToolboxSgp4

using Crayons
using Dates
using LinearAlgebra
using Reexport
using Printf
using StaticArrays
using SatelliteToolboxBase

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

# Escape sequences related to the crayons.
const _D = Crayon(reset = true)
const _B = crayon"bold"
const _Y = crayon"yellow bold"

############################################################################################
#                                         Includes
############################################################################################

include("sgp4_model.jl")
include("tle.jl")

include("precompile.jl")

end
