## Description #############################################################################
#
# Functions to print the structures related to the SGP4 propagator.
#
############################################################################################

function Base.show(io::IO, sgp4d::Sgp4Propagator)
    # The field `algorithm` is only assigned by `sgp4_init!`. Hence, it indicates whether
    # the structure has been initialized.
    if !isdefined(sgp4d, :algorithm)
        print(io, _sgp4_propagator_name(sgp4d), " (not initialized)")
        return nothing
    end

    SatelliteToolboxBase.print_compact(io, _sgp4_propagator_name(sgp4d), sgp4d.epoch)

    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", sgp4d::Sgp4Propagator)
    if !isdefined(sgp4d, :algorithm)
        println(io, _sgp4_propagator_name(sgp4d), ":")
        SatelliteToolboxBase.print_field(io, "  Status : ", "not initialized")
        return nothing
    end

    labels = (
        "Mean motion",
        "Eccentricity",
        "Inclination",
        "RAAN",
        "Arg. of perigee",
        "Mean anomaly",
        "B*",
        "Last propagation",
    )

    values = (
        _sgp4_show_number(720 * sgp4d.n₀ / π),
        _sgp4_show_number(sgp4d.e₀),
        _sgp4_show_number(rad2deg(sgp4d.i₀)),
        _sgp4_show_number(rad2deg(sgp4d.Ω₀)),
        _sgp4_show_number(rad2deg(sgp4d.ω₀)),
        _sgp4_show_number(rad2deg(sgp4d.M₀)),
        SatelliteToolboxBase.compact_string(io, sgp4d.bstar),
        SatelliteToolboxBase.compact_string(io, sgp4d.Δt),
    )

    units = ("rev / day", "", "°", "°", "°", "°", "1 / er", "min")

    SatelliteToolboxBase.print_elements(
        io, _sgp4_propagator_name(sgp4d), sgp4d.epoch, labels, values, units
    )

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _sgp4_algorithm_name(algorithm::Symbol) -> String

Return the name of the propagation `algorithm` selected in `sgp4_init!` to be printed.
"""
function _sgp4_algorithm_name(algorithm::Symbol)
    algorithm === :sgp4 && return "SGP4"
    algorithm === :sgp4_lowper && return "SGP4 (low perigee)"
    algorithm === :sdp4 && return "SDP4"
    return string(algorithm)
end

"""
    _sgp4_propagator_name(sgp4d::Sgp4Propagator{Tepoch, T}) where {Tepoch, T} -> String

Return the name of the propagator structure `sgp4d` with its type parameters, followed by
the selected algorithm in parentheses if `sgp4d` has been initialized, as used in the
headers of the printed representations.
"""
function _sgp4_propagator_name(sgp4d::Sgp4Propagator{Tepoch, T}) where {Tepoch, T}
    name = string("Sgp4Propagator{", Tepoch, ", ", T, "}")
    isdefined(sgp4d, :algorithm) || return name
    return string(name, " (", _sgp4_algorithm_name(sgp4d.algorithm), ")")
end

"""
    _sgp4_show_number(x::Number) -> String

Format the number `x` to be printed with 8 decimal digits if it is a floating-point number.
Otherwise, e.g. for dual numbers, it is printed as is.
"""
_sgp4_show_number(x::AbstractFloat) = @sprintf("%.8f", x)
_sgp4_show_number(x::Number) = string(x)
