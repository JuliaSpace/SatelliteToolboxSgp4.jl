## Description #############################################################################
#
# Functions to initialize the SGP4 orbit propagator using Orbit Mean-Elements Messages
# (OMMs).
#
## References ##############################################################################
#
# [1] CCSDS 502.0-B-3 (2023). Orbit Data Messages. Recommended Standard, Blue Book,
#     Consultative Committee for Space Data Systems.
#
############################################################################################

"""
    sgp4_init(omm::OrbitMeanElementsMessage; kwargs...) -> Sgp4Propagator{Float64, T}

Create and initialize the data structure of the SGP4 orbit propagator using the Orbit
Mean-Elements Message (OMM) `omm` as defined in **[1]**.

The mean element theory of `omm` must be `"SGP4"`. The mean motion is obtained from the
field `MEAN_MOTION` [rev/day] or, if it is absent, from the fields `SEMI_MAJOR_AXIS` [km]
and `GM` [km³/s²]. The epoch of the mean elements is converted to the Julian Day in the
time system of the message, which is UTC for the messages distributed by Space-Track and
Celestrak. If the drag term `BSTAR` [1 / er] is absent, it is set to 0. The function can
fail if the message does not contain the required information.

See also: [`sgp4_init!`](@ref), [`sgp4`](@ref).

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default**: `SGP4C_WGS84`)

# Returns

- `Sgp4Propagator{Float64, T}`: The structure with the initialized parameters.

# References

- **[1]** CCSDS 502.0-B-3 (2023). Orbit Data Messages. Recommended Standard, Blue Book,
    Consultative Committee for Space Data Systems.

# Extended help

## Throws

- `ArgumentError`: If the mean element theory of `omm` is not `"SGP4"`.
- `ArgumentError`: If `omm` provides neither the mean motion nor the semi-major axis
    together with the gravitational coefficient.

## Examples

```julia-repl
julia> omm = read_omm("amazonia_1.xml");

julia> sgp4d = sgp4_init(omm);

julia> r_teme, v_teme = sgp4!(sgp4d, 10)
```
"""
function sgp4_init(
    omm::OrbitMeanElementsMessage; sgp4c::Sgp4Constants{T} = SGP4C_WGS84
) where {T <: Number}
    sgp4d = Sgp4Propagator{Float64}(sgp4c)
    sgp4_init!(sgp4d, omm)
    return sgp4d
end

"""
    sgp4_init!(
        sgp4d::Sgp4Propagator{Tepoch, T},
        omm::OrbitMeanElementsMessage,
    ) where {Tepoch <: Number, T <: Number} -> Nothing

Initialize the SGP4 data structure `sgp4d` with the mean elements in the Orbit
Mean-Elements Message (OMM) `omm` as defined in **[1]**.

The mean element theory of `omm` must be `"SGP4"`. The mean motion is obtained from the
field `MEAN_MOTION` [rev/day] or, if it is absent, from the fields `SEMI_MAJOR_AXIS` [km]
and `GM` [km³/s²]. The epoch of the mean elements is converted to the Julian Day in the
time system of the message, which is UTC for the messages distributed by Space-Track and
Celestrak. If the drag term `BSTAR` [1 / er] is absent, it is set to 0. The function can
fail if the message does not contain the required information.

!!! warning

    The propagation constants `sgp4c::Sgp4Constants` in `sgp4d` will not be changed.
    Hence, they must be initialized.

See also: [`sgp4_init`](@ref).

# References

- **[1]** CCSDS 502.0-B-3 (2023). Orbit Data Messages. Recommended Standard, Blue Book,
    Consultative Committee for Space Data Systems.

# Extended help

## Throws

- `ArgumentError`: If the mean element theory of `omm` is not `"SGP4"`.
- `ArgumentError`: If `omm` provides neither the mean motion nor the semi-major axis
    together with the gravitational coefficient.
"""
function sgp4_init!(
    sgp4d::Sgp4Propagator{Tepoch, T}, omm::OrbitMeanElementsMessage
) where {Tepoch <: Number, T <: Number}
    epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar = _omm_sgp4_elements(omm)

    d2r = T(π / 180)
    sgp4_init!(
        sgp4d,
        epoch,
        n₀ * T(2π / (24 * 60)),
        e₀,
        i₀ * d2r,
        Ω₀ * d2r,
        ω₀ * d2r,
        M₀ * d2r,
        bstar,
    )

    return nothing
end

"""
    sgp4(
        Δt::Number,
        omm::OrbitMeanElementsMessage;
        kwargs...,
    ) -> SVector{3, T}, SVector{3, T}, Sgp4Propagator{Float64, T}

Initialize the SGP4 structure using the Orbit Mean-Elements Message (OMM) `omm` as defined
in **[1]** and propagate the orbit until the time `Δt` [min] from the message epoch.

For more information about the required fields in `omm`, see [`sgp4_init`](@ref). The
function can fail if the message does not contain the required information.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default**: `SGP4C_WGS84`)

# Returns

- `SVector{3, T}`: The position vector represented in the TEME reference frame [km].
- `SVector{3, T}`: The velocity vector represented in the TEME reference frame [km / s].
- `Sgp4Propagator{Float64, T}`: The SGP4 orbit propagator structure.

# References

- **[1]** CCSDS 502.0-B-3 (2023). Orbit Data Messages. Recommended Standard, Blue Book,
    Consultative Committee for Space Data Systems.

# Extended help

## Throws

- `ArgumentError`: If the mean element theory of `omm` is not `"SGP4"`.
- `ArgumentError`: If `omm` provides neither the mean motion nor the semi-major axis
    together with the gravitational coefficient.
"""
function sgp4(
    Δt::Number, omm::OrbitMeanElementsMessage; sgp4c::Sgp4Constants{T} = SGP4C_WGS84
) where {T <: Number}
    epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar = _omm_sgp4_elements(omm)

    d2r = T(π / 180)
    return sgp4(
        Δt,
        epoch,
        n₀ * T(2π / (24 * 60)),
        e₀,
        i₀ * d2r,
        Ω₀ * d2r,
        ω₀ * d2r,
        M₀ * d2r,
        bstar;
        sgp4c = sgp4c,
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _omm_sgp4_elements(omm::OrbitMeanElementsMessage) -> NTuple{8, Float64}

Obtain from the Orbit Mean-Elements Message `omm` the elements required to initialize the
SGP4 orbit propagator, using the same units as the TLE. The function can fail if the
message does not contain the required information.

# Returns

- `Float64`: Epoch of the mean elements [Julian Day].
- `Float64`: Mean motion [rev/day].
- `Float64`: Eccentricity [-].
- `Float64`: Inclination [°].
- `Float64`: Right ascension of the ascending node [°].
- `Float64`: Argument of perigee [°].
- `Float64`: Mean anomaly [°].
- `Float64`: Drag term B* [1 / er].

# Extended help

## Throws

- `ArgumentError`: If the mean element theory of `omm` is not `"SGP4"`.
- `ArgumentError`: If `omm` provides neither the mean motion nor the semi-major axis
    together with the gravitational coefficient.
"""
function _omm_sgp4_elements(omm::OrbitMeanElementsMessage)
    theory = ODM.mean_element_theory(omm)

    theory != "SGP4" && throw(
        ArgumentError(
            "The mean element theory of the OMM must be \"SGP4\", but it is \"$theory\".",
        ),
    )

    # Obtain the mean motion [rev/day], which can be provided directly or computed from the
    # semi-major axis and the gravitational coefficient.
    mean_motion = ODM.mean_motion(omm)

    if isnothing(mean_motion)
        a  = ODM.semi_major_axis(omm)
        GM = ODM.GM(omm)

        (isnothing(a) || isnothing(GM)) && throw(
            ArgumentError(
                "The OMM must provide either the mean motion or the semi-major axis " *
                "together with the gravitational coefficient GM.",
            ),
        )

        mean_motion = √(GM / a^3) * 86400 / 2π
    end

    # The drag term is optional in the OMM. If it is absent, we assume a drag-free
    # propagation.
    bstar = something(ODM.bstar(omm), 0.0)

    return (
        _omm_epoch_to_julian_day(omm),
        mean_motion,
        ODM.eccentricity(omm),
        ODM.inclination(omm),
        ODM.raan(omm),
        ODM.arg_of_pericenter(omm),
        ODM.mean_anomaly(omm),
        bstar,
    )
end

"""
    _omm_epoch_to_julian_day(omm::OrbitMeanElementsMessage) -> Float64

Convert the epoch of the mean elements in the Orbit Mean-Elements Message `omm` to the
Julian Day in the time system of the message, keeping the sub-millisecond information.
"""
function _omm_epoch_to_julian_day(omm::OrbitMeanElementsMessage)
    epoch = ODM.epoch(omm)

    # The `DateTime` conversion truncates the epoch to milliseconds. Hence, we must add the
    # remaining microseconds and nanoseconds provided by the `NanoDate`.
    jd = datetime2julian(DateTime(epoch))
    jd +=
        (1000 * Dates.value(Microsecond(epoch)) + Dates.value(Nanosecond(epoch))) / 86_400e9

    return jd
end
