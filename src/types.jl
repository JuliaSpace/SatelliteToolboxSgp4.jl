# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Types and structures of SGP4 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export Sgp4Constants, Sgp4Propagator

"""
    struct Sgp4Constants{T}

Gravitational constants for SGP4.

# Fields

- `R0::T`: Earth equatorial radius [km].
- `XKE::T`: 60 ⋅ √(GM / R0^3) [er/min]^(3/2).
- `J2::T`: The second gravitational zonal harmonic of the Earth.
- `J3::T`: The thrid gravitational zonal harmonic of the Earth.
- `J4::T`: The fourth gravitational zonal harmonic of the Earth.
"""
struct Sgp4Constants{T}
    R0::T
    XKE::T
    J2::T
    J3::T
    J4::T
end

"""
    struct Sgp4DeepSpace{T}

Store the internal SGP4 variables to account for deep space perturbations.
"""
mutable struct Sgp4DeepSpace{T}
    atime::T
    xli::T
    xni::T
    xnq::T
    xfact::T
    ssl::T
    ssg::T
    ssh::T
    sse::T
    ssi::T
    xlamo::T
    omegaq::T
    omgdt::T
    gmst::T
    del1::T
    del2::T
    del3::T
    fasx2::T
    fasx4::T
    fasx6::T
    d2201::T
    d2211::T
    d3210::T
    d3222::T
    d4410::T
    d4422::T
    d5220::T
    d5232::T
    d5421::T
    d5433::T
    xnddt::T
    xndot::T
    xldot::T
    zmos::T
    se2::T
    se3::T
    si2::T
    si3::T
    sl2::T
    sl3::T
    sl4::T
    sgh2::T
    sgh3::T
    sgh4::T
    sh2::T
    sh3::T
    zmol::T
    ee2::T
    e3::T
    xi2::T
    xi3::T
    xl2::T
    xl3::T
    xl4::T
    xgh2::T
    xgh3::T
    xgh4::T
    xh2::T
    xh3::T
    pe::T
    pinc::T
    pgh::T
    ph::T
    pl::T
    pgh0::T
    ph0::T
    pe0::T
    pinc0::T
    pl0::T

    isynfl::Bool
    iresfl::Bool
    ilsz::Bool

    # Constructors
    # ======================================================================================

    Sgp4DeepSpace{T}(args...) where T<:Number = new(args...)
    Sgp4DeepSpace{T}() where T<:Number = new()
end

"""
    Sgp4Propagator{Tepoch, T}

Low-level SGP4 propagator structure.
"""
mutable struct Sgp4Propagator{Tepoch<:Number, T<:Number}
    # TLE parameters.
    epoch::Tepoch
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    M_0::T
    bstar::T
    # Propagation time from epoch.
    Δt::T
    # Current mean orbital elements.
    a_k::T
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    M_k::T
    n_k::T
    # Parameters related with the orbit.
    all_0::T
    nll_0::T
    # Useful constants to decrease the computational burden.
    AE::T
    QOMS2T::T
    β_0::T
    ξ::T
    η::T
    sin_i_0::T
    θ::T
    θ²::T
    A_30::T
    k_2::T
    k_4::T
    C1::T
    C3::T
    C4::T
    C5::T
    D2::T
    D3::T
    D4::T
    dotM::T
    dotω::T
    dotΩ1::T
    dotΩ::T
    # Selected algorithm.
    algorithm::Symbol
    # SGP4 gravitational constants.
    sgp4c::Sgp4Constants{T}
    # SGP4 deep space structure.
    sgp4ds::Sgp4DeepSpace{T}

    # Constructors
    # ======================================================================================

    Sgp4Propagator{Tepoch, T}(args...) where {Tepoch<:Number, T<:Number} = new(args...)
    Sgp4Propagator{Tepoch, T}() where {Tepoch<:Number, T<:Number} = new()
end
