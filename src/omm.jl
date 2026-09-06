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

# == Fitting Support =======================================================================

"""
    _mean_elements_epoch(omm::OrbitMeanElementsMessage) -> Float64

Return the epoch of the mean elements in `omm` [Julian Day] in the time system of the
message.
"""
_mean_elements_epoch(omm::OrbitMeanElementsMessage) = _omm_epoch_to_julian_day(omm)

"""
    _mean_state_vector(omm::OrbitMeanElementsMessage; sgp4c::Sgp4Constants{T} = SGP4C_WGS84) where {T <: Number} -> SVector{7, T}

Convert the Orbit Mean-Elements Message `omm` to the SGP4 mean state vector using the
constants `sgp4c`. The function can fail if the message does not contain the required
information (see [`_omm_sgp4_elements`](@ref)). The state vector has the following
structure:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘
"""
function _mean_state_vector(
    omm::OrbitMeanElementsMessage; sgp4c::Sgp4Constants{T} = SGP4C_WGS84
) where {T}
    ~, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar = _omm_sgp4_elements(omm)
    return _elements_to_mean_state_vector(n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar, sgp4c)
end

"""
    _build_mean_elements(::Type{OrbitMeanElementsMessage}, sv::SVector{7}, epoch::Number; kwargs...) -> OrbitMeanElementsMessage

Create an Orbit Mean-Elements Message (OMM) for the `epoch` [Julian Day, UTC] given the
SGP4 mean state vector `sv`, which must have the following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

The creation date of the message is set to the current time [UTC]. The first and second
time derivatives of the mean motion are set to 0, since they are not estimated.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 propagator constants.
    (**Default**: `SGP4C_WGS84`)
- `template::Union{Nothing, OrbitMeanElementsMessage, NamedTuple}`: Message from which
    the header, the metadata, the spacecraft parameters, and the TLE-related parameters
    are copied, or a `NamedTuple` whose entries are passed as keywords to the constructor
    of `OrbitMeanElementsMessage` on top of the default metadata. The mean element theory
    and the reference frame are always set to `"SGP4"` and `"TEME"`. If it is `nothing`,
    the metadata is filled with default values.
    (**Default**: `nothing`)
- `covariance::Union{Nothing, SMatrix{6, 6}}`: Covariance matrix of the mean position [km]
    and velocity [km / s] represented in the TEME reference frame, stored in the covariance
    matrix section of the message. If it is `nothing`, the section is omitted.
    (**Default**: `nothing`)

# Extended help

## Throws

- `ArgumentError`: If a `NamedTuple` template contains a field set by the fit.
"""
function _build_mean_elements(
    ::Type{OrbitMeanElementsMessage},
    sv::SVector{7},
    epoch::Number;
    sgp4c::Sgp4Constants = SGP4C_WGS84,
    template::Union{Nothing, OrbitMeanElementsMessage, NamedTuple} = nothing,
    covariance::Union{Nothing, SMatrix{6, 6}} = nothing,
)
    n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar = _mean_state_vector_to_elements(sv, sgp4c)

    # Assemble the covariance matrix section, if requested.
    P = covariance

    covariance_matrix =
        isnothing(P) ? nothing :
        OmmCovarianceMatrix(;
            cov_ref_frame = "TEME",
            cx_x = P[1, 1],
            cy_x = P[2, 1],
            cy_y = P[2, 2],
            cz_x = P[3, 1],
            cz_y = P[3, 2],
            cz_z = P[3, 3],
            cx_dot_x = P[4, 1],
            cx_dot_y = P[4, 2],
            cx_dot_z = P[4, 3],
            cx_dot_x_dot = P[4, 4],
            cy_dot_x = P[5, 1],
            cy_dot_y = P[5, 2],
            cy_dot_z = P[5, 3],
            cy_dot_x_dot = P[5, 4],
            cy_dot_y_dot = P[5, 5],
            cz_dot_x = P[6, 1],
            cz_dot_y = P[6, 2],
            cz_dot_z = P[6, 3],
            cz_dot_x_dot = P[6, 4],
            cz_dot_y_dot = P[6, 5],
            cz_dot_z_dot = P[6, 6],
        )

    # Keywords with the fitted values, which are the same regardless of the template.
    fitted = (;
        creation_date     = NanoDate(now(UTC)),
        ref_frame         = "TEME",
        epoch             = NanoDate(julian2datetime(epoch)),
        semi_major_axis   = nothing,
        mean_motion       = n₀,
        eccentricity      = e₀,
        inclination       = i₀,
        raan              = Ω₀,
        arg_of_pericenter = ω₀,
        mean_anomaly      = M₀,
        bstar             = bstar,
        bterm             = nothing,
        mean_motion_dot   = 0.0,
        mean_motion_ddot  = 0.0,
        agom              = nothing,
        covariance_matrix = covariance_matrix,
    )

    template isa OrbitMeanElementsMessage &&
        return OrbitMeanElementsMessage(template; mean_element_theory = "SGP4", fitted...)

    # The default metadata mirrors the default TLE fields so that the message can be
    # converted to a TLE.
    metadata = (;
        originator          = "SatelliteToolboxSgp4.jl",
        object_name         = "UNDEFINED",
        object_id           = "UNDEFINED",
        center_name         = "EARTH",
        time_system         = "UTC",
        mean_element_theory = "SGP4",
        ephemeris_type      = 0,
        classification_type = 'U',
        norad_cat_id        = 9999,
        element_set_number  = 0,
        rev_at_epoch        = 0,
    )

    if template isa NamedTuple
        _check_template(template, (keys(fitted)..., :mean_element_theory))
        metadata = merge(metadata, template)
    end

    return OrbitMeanElementsMessage(; metadata..., fitted...)
end
