module SatelliteToolboxSgp4

using Dates
using LinearAlgebra
using Printf
using StyledStrings

using ForwardDiff
using Reexport
using StaticArrays

@reexport using SatelliteToolboxBase
@reexport using SatelliteToolboxTle
@reexport using SatelliteToolboxOrbitDataMessages

import Base: copy

############################################################################################
#                                          Types                                           #
############################################################################################

include("types.jl")

############################################################################################
#                                        Constants                                         #
############################################################################################

# Julian Day related to the epoch 1900-01-01T12:00:00.000.
const _JD_1900 = DateTime(1900, 1, 1, 12, 0, 0) |> datetime2julian

# The decorated strings printed by the TLE fitting algorithm are rendered once here with
# **StyledStrings.jl**, in the plain and in the colored versions, so that the algorithm only
# prints plain strings. Otherwise, printing an annotated string inside the algorithm adds
# many allocation sites to it, even though they are only reachable when `verbose` is `true`.
const _FIT_ACTION_TAG = (
    "ACTION:",
    sprint(
        print, styled"{(foreground=yellow,weight=bold):ACTION:}"; context = :color => true
    ),
)

const _FIT_PROGRESS_TAG = (
    "PROGRESS:", sprint(print, styled"{bold:PROGRESS:}"; context = :color => true)
)

const _FIT_HEADER = let
    header = @sprintf(
        "%10s %20s %20s %20s %20s",
        "Iteration",
        "Position RMSE",
        "Velocity RMSE",
        "Total RMSE",
        "RMSE Variation"
    )
    (
        header,
        sprint(
            print,
            styled"{(foreground=yellow,weight=bold):$header}";
            context = :color => true,
        ),
    )
end

const _FIT_UNITS = let
    units = @sprintf("%10s %20s %20s %20s", "", "[km]", "[km / s]", "[ ]")
    (units, sprint(print, styled"{bold:$units}"; context = :color => true))
end

############################################################################################
#                                         Includes                                         #
############################################################################################

include("copy.jl")
include("sgp4_model.jl")
include("tle.jl")
include("omm.jl")

include("precompile.jl")

end
