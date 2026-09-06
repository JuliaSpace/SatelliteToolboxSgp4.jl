## Description #############################################################################
#
# Functions to convert between TLEs and the SGP4 mean state vector used by the fitting
# algorithm.
#
############################################################################################

"""
    _mean_elements_epoch(tle::TLE) -> Float64

Return the epoch of the mean elements in `tle` [Julian Day, UTC].
"""
_mean_elements_epoch(tle::TLE) = tle_epoch(tle)

"""
    _mean_state_vector(tle::TLE; sgp4c::Sgp4Constants{T} = SGP4C_WGS84) where {T <: Number} -> SVector{7, T}

Convert the `tle` to the SGP4 mean state vector using the constants `sgp4c`. The state
vector has the following structure:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘
"""
function _mean_state_vector(tle::TLE; sgp4c::Sgp4Constants{T} = SGP4C_WGS84) where {T}
    return _elements_to_mean_state_vector(
        tle.mean_motion,
        tle.eccentricity,
        tle.inclination,
        tle.raan,
        tle.argument_of_perigee,
        tle.mean_anomaly,
        tle.bstar,
        sgp4c,
    )
end

"""
    _build_mean_elements(::Type{TLE}, sv::SVector{7}, epoch::Number; kwargs...) -> TLE

Create a TLE for the `epoch` [Julian Day, UTC] given the SGP4 mean state vector `sv`, which
must have the following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 propagator constants.
    (**Default**: `SGP4C_WGS84`)
- `template::Union{Nothing, TLE}`: TLE from which the satellite name, number,
    classification, international designator, element set number, and revolution number
    are copied. If it is `nothing`, default values are used.
    (**Default**: `nothing`)
- `covariance::Union{Nothing, SMatrix{6, 6}}`: Not used, since a TLE cannot store the
    covariance matrix.
    (**Default**: `nothing`)
"""
function _build_mean_elements(
    ::Type{TLE},
    sv::SVector{7},
    epoch::Number;
    sgp4c::Sgp4Constants = SGP4C_WGS84,
    template::Union{Nothing, TLE} = nothing,
    covariance::Union{Nothing, SMatrix{6, 6}} = nothing,
)
    n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar = _mean_state_vector_to_elements(sv, sgp4c)

    # Compute the epoch as required by the TLE format.
    dt  = julian2datetime(epoch)
    dt₀ = DateTime(Year(dt))

    dt_year    = year(dt)
    epoch_year = dt_year < 2000 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt - dt₀).value / 1000 / 86400 + 1

    return TLE(;
        name                     = isnothing(template) ? "UNDEFINED" : template.name,
        satellite_number         = isnothing(template) ? 9999 : template.satellite_number,
        classification           = isnothing(template) ? 'U' : template.classification,
        international_designator = isnothing(template) ? "999999" : template.international_designator,
        epoch_year               = epoch_year,
        epoch_day                = epoch_day,
        dn_o2                    = 0,
        ddn_o6                   = 0,
        bstar                    = bstar,
        element_set_number       = isnothing(template) ? 0 : template.element_set_number,
        inclination              = i₀,
        raan                     = Ω₀,
        eccentricity             = e₀,
        argument_of_perigee      = ω₀,
        mean_anomaly             = M₀,
        mean_motion              = n₀,
        revolution_number        = isnothing(template) ? 0 : template.revolution_number,
    )
end
