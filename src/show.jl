## Description #############################################################################
#
# Functions to print the structures related to the SGP4 propagator.
#
# The representations follow the layout of SatelliteToolboxBase.jl: the compact form prints
# the type with its parameters and the epoch, whereas the rich form is a tree with the
# sections holding the mean elements and their epoch, the gravitational constants, and the
# last propagation instant.
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
    SatelliteToolboxBase.print_tree(io, _sgp4_propagator_name(sgp4d), sgp4d)
    return nothing
end

# The body of the rich representation is overloaded so that the wrappers of the propagator
# can print it under their own header.
function SatelliteToolboxBase.print_tree_body(io::IO, sgp4d::Sgp4Propagator)
    if !isdefined(sgp4d, :algorithm)
        fields   = SatelliteToolboxBase.PrintedField[("Status", "not initialized", "")]
        sections = SatelliteToolboxBase.PrintedSection[]
        SatelliteToolboxBase.print_tree_body(io, fields, sections)
        return nothing
    end

    sgp4c = sgp4d.sgp4c

    format_value = SatelliteToolboxBase.format_value

    # The semi-major axis is recovered from the mean motion as in the SGP4 theory.
    semi_major_axis = (sgp4c.XKE / sgp4d.n₀)^(2 // 3) * sgp4c.R0

    epoch_str = SatelliteToolboxBase.epoch_string(sgp4d.epoch)

    sections = SatelliteToolboxBase.PrintedSection[
        "Mean Elements" => [
            ("Epoch",             epoch_str,                        ""),
            ("Semi-Major Axis",   format_value(semi_major_axis),    "km"),
            ("Mean Motion",       format_value(720 * sgp4d.n₀ / π), "rev/day"),
            ("Eccentricity",      format_value(sgp4d.e₀),           ""),
            ("Inclination",       format_value(rad2deg(sgp4d.i₀)),  "°"),
            ("RA of Asc. Node",   format_value(rad2deg(sgp4d.Ω₀)),  "°"),
            ("Arg. of Periapsis", format_value(rad2deg(sgp4d.ω₀)),  "°"),
            ("Mean Anomaly",      format_value(rad2deg(sgp4d.M₀)),  "°"),
            ("B*",                format_value(sgp4d.bstar),        "1/er"),
        ],
        "Constants" => [
            ("R₀",  format_value(sgp4c.R0),  "km"),
            ("XKE", format_value(sgp4c.XKE), "er^(3/2)/min"),
            ("J₂",  format_value(sgp4c.J2),  ""),
            ("J₃",  format_value(sgp4c.J3),  ""),
            ("J₄",  format_value(sgp4c.J4),  ""),
        ],
        "Propagation" => [
            ("Last Instant", format_value(sgp4d.Δt), "min"),
        ],
    ]

    SatelliteToolboxBase.print_tree_body(io, SatelliteToolboxBase.PrintedField[], sections)

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
