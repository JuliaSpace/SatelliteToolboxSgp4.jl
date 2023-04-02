# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Helpers for SGP4 algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sgp4

"""
    sgp4(Δt, args...; kwargs...)

Initialize the SGP4 structure and propagate the orbit until the time Δt [s].

# Returns

- `SVector{3, T}`: The position vector [km].
- `SVector{3, T}`: The velocity vector [km/s].
- [`Sgp4Propagator`](@ref): The SGP4 orbit propagator structure.
"""
function sgp4(Δt, args...; kwargs...)
    sgp4d = sgp4_init(args...; kwargs...)
    r_teme, v_teme = sgp4!(sgp4d, Δt)
    return r_teme, v_teme, sgp4d
end
