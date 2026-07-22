## Description #############################################################################
#
# SGP4 orbit propagator model.
#
# This is a independent implementation of the algorithm presented in [1].
#
## References ##############################################################################
#
# [1] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD Elements Set.
#     Spacetrack Report No. 3.
#
# [2] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting Spacetrack
#     Report #3: Rev1. AIAA.
#
# [3] SGP4 Source code of STRF: https://github.com/cbassa/strf
#     The SGP4 C code available on STRF was converted by Paul. S. Crawford and Andrew R.
#     Brooks.
#
############################################################################################

export sgp4c_wgs72, sgp4c_wgs84
export sgp4c_wgs72_f32, sgp4c_wgs84_f32
export sgp4_init, sgp4_init!, sgp4, sgp4!


############################################################################################
#                                        Constants                                         #
############################################################################################

# WGS-84 / EGM-08 gravitational constants.
const sgp4c_wgs84 = Sgp4Constants{Float64}(
    6378.137,
    60.0 / sqrt(6378.137^3 / 398600.5),
     0.00108262998905,
    -0.00000253215306,
    -0.00000161098761
)

const sgp4c_wgs84_f32 = Sgp4Constants{Float32}(
    6378.137,
    60.0 / sqrt(6378.137^3 / 398600.5),
     0.00108262998905,
    -0.00000253215306,
    -0.00000161098761
)

# WGS-72 gravitational constants.
const sgp4c_wgs72 = Sgp4Constants{Float64}(
    6378.135,
    60.0 / sqrt(6378.135^3 / 398600.8),
     0.001082616,
    -0.00000253881,
    -0.00000165597
)

const sgp4c_wgs72_f32 = Sgp4Constants{Float32}(
    6378.135,
    60.0 / sqrt(6378.135^3 / 398600.8),
     0.001082616,
    -0.00000253881,
    -0.00000165597
)

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    sgp4_init(epoch::Tepoch, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) where {Tepoch<:Number, T<:Number}
    sgp4_init(tle::TLE; kwargs...) where T

Create and initialize the data structure of SGP4 orbit propagator.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/min].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: TLE to initialize the SGP4 (see `TLE`).

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)

# Returns

- [`Sgp4Propagator`](@ref): The structure with the initialized parameters.
"""
function sgp4_init(tle::TLE; sgp4c::Sgp4Constants{T} = sgp4c_wgs84) where T<:Number
    # We must initialize the SGP4 propagator structure together with any mutable fields.
    Tepoch = typeof(tle_epoch(tle))
    sgp4d = Sgp4Propagator{Tepoch, T}()
    sgp4d.sgp4c = sgp4c
    sgp4d.sgp4ds = Sgp4DeepSpace{T}()

    sgp4_init!(sgp4d, tle)
    return sgp4d
end

function sgp4_init(
    epoch::Tepoch,
    n₀::N,
    e₀::E,
    i₀::I,
    Ω₀::O,
    ω₀::W,
    M₀::M,
    bstar::B;
    sgp4c::Sgp4Constants{T} = sgp4c_wgs84
) where {
    Tepoch<:Number,
    N<:AbstractFloat, E<:AbstractFloat, I<:AbstractFloat,
    O<:AbstractFloat, W<:AbstractFloat, M<:AbstractFloat, B<:AbstractFloat,
    T<:Number
}
    sgp4d = Sgp4Propagator{Tepoch, T}()
    sgp4d.sgp4c = Sgp4Constants{T}(sgp4c.R0, sgp4c.XKE, sgp4c.J2, sgp4c.J3, sgp4c.J4)
    sgp4d.sgp4ds = Sgp4DeepSpace{T}()

    sgp4_init!(sgp4d, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)
    return sgp4d
end

function sgp4_init(
    epoch::Tepoch,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants{T} = sgp4c_wgs84
) where {Tepoch<:Number, T<:Number}
    Tprom = promote_type(
        T, typeof(n₀), typeof(e₀), typeof(i₀),
        typeof(Ω₀), typeof(ω₀), typeof(M₀), typeof(bstar)
    )

    sgp4d = Sgp4Propagator{Tepoch, Tprom}()
    sgp4d.sgp4c = Sgp4Constants{Tprom}(sgp4c.R0, sgp4c.XKE, sgp4c.J2, sgp4c.J3, sgp4c.J4)
    sgp4d.sgp4ds = Sgp4DeepSpace{Tprom}()

    sgp4_init!(sgp4d, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)
    return sgp4d
end

"""
    sgp4_init!(sgp4d::Sgp4Propagator{Tepoch, T}, epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number) where {Tepoch, T} -> Nothing
    sgp4_init!(sgp4d::Sgp4Propagator{Tepoch, T}, tle::TLE) where {Tepoch, T} -> Nothing

Initialize the SGP4 data structure `sgp4d` with the initial orbit specified by the
arguments.

!!! warning

    The propagation constants `sgp4c::Sgp4PropagatorConstants` in `sgp4d` will not be
    changed. Hence, they must be initialized.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/min].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: TLE to initialize the SPG4 (see `TLE`).
"""
function sgp4_init!(
    sgp4d::Sgp4Propagator{Tepoch, T},
    tle::TLE
) where {Tepoch<:Number, T<:Number}
    d2r = T(π / 180)
    sgp4_init!(
        sgp4d,
        tle_epoch(tle),
        tle.mean_motion * T(2π / (24 * 60)),
        tle.eccentricity,
        tle.inclination * d2r,
        tle.raan * d2r,
        tle.argument_of_perigee * d2r,
        tle.mean_anomaly * d2r,
        tle.bstar
    )

    return nothing
end

function sgp4_init!(
    sgp4d::Sgp4Propagator{Tepoch, ST},
    epoch::EpT,
    n₀::NT,
    e₀::ET,
    i₀::IT,
    Ω₀::OT,
    ω₀::WT,
    M₀::MT,
    bstar::BT
) where {Tepoch<:Number, EpT<:Number, NT<:Number, ET<:Number, IT<:Number, OT<:Number, WT<:Number, MT<:Number, BT<:Number, ST<:Number}

    T = ST

    # Unpack the gravitational constants to improve code readability.
    sgp4c = sgp4d.sgp4c
    R0    = sgp4c.R0
    XKE   = sgp4c.XKE
    J2    = sgp4c.J2
    J3    = sgp4c.J3
    J4    = sgp4c.J4

    # == Constants =========================================================================
    #
    # Note: [er] = Earth radii.

    # Distance units / Earth radii.
    AE = T(1)

    k₂  = +(1 // 2) * J2 * AE * AE
    k₂² = k₂ * k₂
    k₄  = -(3 // 8) * J4 * AE * AE * AE * AE
    A₃₀ = -J3 * AE * AE * AE

    # Kilometers / Earth radii.
    XKMPER = R0

    # Parameters for the SGP4 density function.
    s  =  78 / XKMPER + 1
    q₀ = 120 / XKMPER + 1

    # (q₀ - s)^4 [er]^4
    QOMS2T = (q₀ - s) * (q₀ - s) * (q₀ - s) * (q₀ - s)

    # == Auxiliary Variables to Improve the Performance ====================================

    e₀² = T(e₀)^2

    sin_i₀, θ = sincos(T(i₀))
    θ²        = θ  * θ
    θ³        = θ² * θ
    θ⁴        = θ² * θ²

    # ======================================================================================

    # Recover the original mean motion (nll₀) and semi-major axis (all₀) from the input
    # elements.
    aux = (3θ² - 1) / √((1 - e₀²)^3)
    a₁  = (XKE / T(n₀))^(T(2 / 3))
    δ₁  = (3 // 2) * k₂ / (a₁ * a₁)* aux
    a₀  = a₁ * @evalpoly(δ₁, 1, -(1 // 3), -1, -(134 // 81))
    δ₀  = (3 // 2) * k₂ / (T(a₀) * T(a₀)) * aux

    nll₀ = T(n₀) / (1 + δ₀)

    # Vallado's implementation of SGP4 [2] compute the semi-major axis considering the new
    # angular velocity, which is called `no_unkozai`. In the original SGP4 technical report
    # [1], the semi-major axis was computed considering:
    #
    #   all₀ = a₀ / (1 - δ₀)
    #

    all₀  = (XKE / nll₀)^(2 // 3)
    all₀² = all₀  * all₀
    all₀⁴ = all₀² * all₀²

    # == Initialization ====================================================================

    # Compute the orbit perigee [ER].
    perigee = (all₀ * (1 - T(e₀)) - AE) * XKMPER

    # For perigee below 156 km, the values of S and QOMS2T are altered.
    if perigee < 156
        if perigee < 98
            s = 20 / XKMPER + AE
        # Perigee between 98km and 156km.
        else
            s = all₀ * (1 - T(e₀)) - s + AE
        end

        QOMS2T = (q₀ - s) * (q₀ - s) * (q₀ - s) * (q₀ - s)
    end

    # Compute SGP4 constants.
    ξ  = 1 / (all₀ - s)
    ξ² = ξ  * ξ
    ξ³ = ξ² * ξ
    ξ⁴ = ξ² * ξ²
    ξ⁵ = ξ⁴ * ξ

    β₀  = √(1 - e₀²)
    β₀² = β₀  * β₀
    β₀³ = β₀² * β₀
    β₀⁴ = β₀² * β₀²
    β₀⁷ = β₀⁴ * β₀³
    β₀⁸ = β₀⁴ * β₀⁴

    η  = all₀ * T(e₀) * ξ
    η² = η  * η
    η³ = η² * η
    η⁴ = η² * η²

    # Vallado's implementation of SGP4 [2] considers the absolute value of (1-η^2) here and
    # in the C2 and C4 computation. Notice that, if (1-η^2) < 0, then aux1 cannot be
    # computed. The original SGP4 technical report [1] does not mention anything about this.

    aux0 = abs(1 - η²)
    aux1 = 1 / (√aux0^7)           # ......................................... aux0^(-7 / 2)
    aux2 = ξ⁴ * all₀ * β₀² * aux1

    C2 = QOMS2T * ξ⁴ * nll₀ * aux1 * (
        all₀ * (1 + (3 // 2) * η² + 4T(e₀) * η + T(e₀) * η³) +
        (3 // 2) * (k₂ * ξ) / aux0 * (-(1 // 2) + (3 // 2) * θ²) * (8 + 24η² + 3η⁴)
    )

    C1  = T(bstar) * C2
    C1² = C1  * C1
    C1³ = C1² * C1
    C1⁴ = C1² * C1²

    C3 = T(e₀) > 1e-4 ? QOMS2T * ξ⁵ * A₃₀ * nll₀ * AE * sin_i₀ / (k₂ * T(e₀)) : T(0)

    C4 = 2nll₀ * QOMS2T * aux2 * (
        2η * (1 + T(e₀) * η) + (1 // 2) * (T(e₀) + η³) -
        2k₂ * ξ / (all₀ * aux0) * (
            3 * (1 - 3θ²) * (1 + (3 // 2) * η² - 2T(e₀) * η - (1 // 2) * T(e₀) * η³) +
            (3 // 4) * (1 - θ²) * (2η² - T(e₀) * η - T(e₀) * η³) * cos(2T(ω₀))
        )
    )

    C5 = 2QOMS2T * aux2 * (1 + (11 // 4) * η * (η + T(e₀)) + T(e₀) * η³)

    D2 = 4all₀ * ξ * C1²

    D3 = T(4 / 3) * all₀ * ξ² * (17all₀ + s) * C1³

    # Vallado's implementation of SGP4 [2] uses all₀^2, instead of only all₀ that is seen
    # in the original SGP4 Technical Report [1].
    D4 = T(2 / 3) * all₀² * ξ³ * (221all₀ + 31s) * C1⁴

    # Compute the time-derivative of some orbital elements.
    ∂M = (
        1 + 3k₂  * (-1 +  3θ²        ) / ( 2all₀² * β₀³) +
            3k₂² * (13 - 78θ² + 137θ⁴) / (16all₀⁴ * β₀⁷)
    ) * nll₀

    ∂ω = (
        -3k₂  * (1 -   5θ²        ) / ( 2all₀² * β₀⁴) +
         3k₂² * (7 - 114θ² + 395θ⁴) / (16all₀⁴ * β₀⁸) +
         5k₄  * (3 -  36θ² +  49θ⁴) / ( 4all₀⁴ * β₀⁸)
    ) * nll₀

    ∂Ω1 = -3k₂ * θ / (all₀² * β₀⁴) * nll₀

    ∂Ω  = ∂Ω1 + (
        3k₂² * (4θ - 19θ³) / (2all₀⁴ * β₀⁸) +
        5k₄  * (3θ -  7θ³) / (2all₀⁴ * β₀⁸)
    ) * nll₀

    # If the orbit period is higher than 225 min., then we must consider the deep space
    # perturbations. This is indicated by selecting the algorithm `:sdp4`.
    if 2π / T(n₀) >= 225.0
        algorithm = :sdp4

        # Initialize the values for the SDP4 (deep space) algorithm.
        _dsinit!(
            sgp4d.sgp4ds,
            Tepoch(epoch),
            nll₀,
            all₀,
            T(e₀),
            T(i₀),
            T(Ω₀),
            T(ω₀),
            T(M₀),
            ∂M,
            ∂ω,
            ∂Ω
        )
    else
        # For perigee lower than 220 km, the equations are truncated to a linear variation
        # in `sqrt(a)` and quadratic variation in mean anomaly. Also, the C5 term, the δω
        # term, and the δM term are dropped. This is indicated by selecting the algorithm
        # `:sgp4_lowper`. Otherwise, if perigee is higher or equal 220 km and the orbit
        # period is lower than 225 min., then we use the normal SGP4 algorithm by selecting
        # `:sgp4`.
        algorithm = (perigee / AE >= (220 + (AE - 1) * XKMPER)) ? :sgp4 : :sgp4_lowper
    end

    # Initialize the structure with the data.
    sgp4d.epoch     = epoch
    sgp4d.n₀        = n₀
    sgp4d.e₀        = e₀
    sgp4d.i₀        = i₀
    sgp4d.Ω₀        = Ω₀
    sgp4d.ω₀        = ω₀
    sgp4d.M₀        = M₀
    sgp4d.bstar     = bstar
    sgp4d.Δt        = 0
    sgp4d.a_k       = all₀
    sgp4d.e_k       = e₀
    sgp4d.i_k       = i₀
    sgp4d.Ω_k       = Ω₀
    sgp4d.ω_k       = ω₀
    sgp4d.M_k       = M₀
    sgp4d.n_k       = nll₀
    sgp4d.all₀      = all₀
    sgp4d.nll₀      = nll₀
    sgp4d.AE        = AE
    sgp4d.QOMS2T    = QOMS2T
    sgp4d.β₀        = β₀
    sgp4d.ξ         = ξ
    sgp4d.η         = η
    sgp4d.sin_i₀    = sin_i₀
    sgp4d.θ         = θ
    sgp4d.θ²        = θ²
    sgp4d.A₃₀       = A₃₀
    sgp4d.k₂        = k₂
    sgp4d.k₄        = k₄
    sgp4d.C1        = C1
    sgp4d.C3        = C3
    sgp4d.C4        = C4
    sgp4d.C5        = C5
    sgp4d.D2        = D2
    sgp4d.D3        = D3
    sgp4d.D4        = D4
    sgp4d.∂M        = ∂M
    sgp4d.∂ω        = ∂ω
    sgp4d.∂Ω        = ∂Ω
    sgp4d.algorithm = algorithm

    return nothing
end

"""
    sgp4(Δt::Number, tle::TLE; kwargs...)
    sgp4(epoch::Tepoch, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) where {Tepoch<:Number, T<:Number}

Initialize the SGP4 structure and propagate the orbit until the time Δt [min].

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/min].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: TLE to initialize the SGP4 (see `TLE`).

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)

# Returns

- `SVector{3, T}`: The position vector [km].
- `SVector{3, T}`: The velocity vector [km/s].
- [`Sgp4Propagator`](@ref): The SGP4 orbit propagator structure.
"""
function sgp4(Δt::Number, tle::TLE; sgp4c::Sgp4Constants{T} = sgp4c_wgs84) where T<:Number
    d2r = T(π / 180)
    return sgp4(
        Δt,
        tle_epoch(tle),
        tle.mean_motion * T(2π / (24 * 60)),
        tle.eccentricity,
        tle.inclination * d2r,
        tle.raan * d2r,
        tle.argument_of_perigee * d2r,
        tle.mean_anomaly * d2r,
        tle.bstar;
        sgp4c = sgp4c
    )
end

function sgp4(
    Δt::D,
    epoch::Tepoch,
    n₀::N,
    e₀::E,
    i₀::I,
    Ω₀::O,
    ω₀::W,
    M₀::M,
    bstar::B;
    sgp4c::Sgp4Constants{T} = sgp4c_wgs84
) where {
    Tepoch<:Number,
    D<:AbstractFloat, N<:AbstractFloat, E<:AbstractFloat, I<:AbstractFloat,
    O<:AbstractFloat, W<:AbstractFloat, M<:AbstractFloat, B<:AbstractFloat,
    T<:Number
}
    sgp4d = sgp4_init(epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar; sgp4c)
    r_teme, v_teme = sgp4!(sgp4d, Δt)
    return r_teme, v_teme, sgp4d
end

function sgp4(
    Δt::Number,
    epoch::Tepoch,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants{T} = sgp4c_wgs84
) where {Tepoch<:Number, T<:Number}
    Tprom = promote_type(
        T, typeof(Δt), typeof(n₀), typeof(e₀), typeof(i₀),
        typeof(Ω₀), typeof(ω₀), typeof(M₀), typeof(bstar)
    )
    sgp4c_p = Sgp4Constants{Tprom}(sgp4c.R0, sgp4c.XKE, sgp4c.J2, sgp4c.J3, sgp4c.J4)

    sgp4d = sgp4_init(epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar; sgp4c = sgp4c_p)
    r_teme, v_teme = sgp4!(sgp4d, Δt)
    return r_teme, v_teme, sgp4d
end

"""
    sgp4!(sgp4d::Sgp4Propagator{Tepoch, T}, t::Number) where T

Propagate the orbit defined in `sgp4d` (see [`Sgp4Propagator`](@ref)) until the time `t`
[min].

!!! note

    The internal values in `sgp4d` will be modified.

# Returns

- `SVector{3, T}`: The position vector represented in TEME frame at time `t` [km].
- `SVector{3, T}`: The velocity vector represented in TEME frame at time `t` [km/s].
"""
function sgp4!(sgp4d::Sgp4Propagator{Tepoch, T}, t::Number) where {Tepoch, T}
    # Unpack variables.
    e₀        = sgp4d.e₀
    i₀        = sgp4d.i₀
    Ω₀        = sgp4d.Ω₀
    ω₀        = sgp4d.ω₀
    M₀        = sgp4d.M₀
    bstar     = sgp4d.bstar
    a_k       = sgp4d.a_k
    e_k       = sgp4d.e_k
    i_k       = sgp4d.i_k
    Ω_k       = sgp4d.Ω_k
    ω_k       = sgp4d.ω_k
    M_k       = sgp4d.M_k
    n_k       = sgp4d.n_k
    all₀      = sgp4d.all₀
    nll₀      = sgp4d.nll₀
    AE        = sgp4d.AE
    QOMS2T    = sgp4d.QOMS2T
    β₀        = sgp4d.β₀
    ξ         = sgp4d.ξ
    η         = sgp4d.η
    sin_i₀    = sgp4d.sin_i₀
    θ         = sgp4d.θ
    θ²        = sgp4d.θ²
    A₃₀       = sgp4d.A₃₀
    k₂        = sgp4d.k₂
    C1        = sgp4d.C1
    C3        = sgp4d.C3
    C4        = sgp4d.C4
    C5        = sgp4d.C5
    D2        = sgp4d.D2
    D3        = sgp4d.D3
    D4        = sgp4d.D4
    ∂M        = sgp4d.∂M
    ∂ω        = sgp4d.∂ω
    ∂Ω        = sgp4d.∂Ω
    algorithm = sgp4d.algorithm
    sgp4c     = sgp4d.sgp4c
    sgp4ds    = sgp4d.sgp4ds

    R0  = sgp4c.R0
    XKE = sgp4c.XKE

    # After unpacking sgp4d, we have two sets of orbit elements:
    #
    #   (n₀, e₀, i₀, Ω₀, ω₀, M₀),
    #
    # and
    #
    #   (n_k, e_k, i_k, Ω_k, ω_k, M_k).
    #
    # The first are those initial elements from the orbit defined in `sgp4_init` function.
    # The second are the current elements. During this functions, the second set is updated
    # by adding the many effects considered in SGP4.

    # Time elapsed since epoch.
    #
    # We convert to `T` to avoid numerical problems with very big numbers as pointed out in:
    #
    #   https://github.com/JuliaLang/julia/issues/27355
    Δt = T(t)

    # Initialization of the current elements with the values of the epoch.
    n_k = nll₀
    a_k = all₀
    e_k = e₀
    i_k = i₀
    Ω_k = Ω₀
    ω_k = ω₀
    M_k = M₀

    # Auxiliary variables to improve code performance.
    sin_i_k = sin_i₀

    # == Secular Effects of Atmospheric Drag and Gravitation ===============================

    M_k = M₀ + ∂M * Δt
    Ω_k = Ω₀ + ∂Ω * Δt - (21 // 2) * (nll₀ * k₂ * θ) / (all₀^2 * β₀^2) * C1 * Δt^2
    ω_k = ω₀ + ∂ω * Δt

    # Check if we need to use SDP4 (deep space) algorithm.
    if algorithm == :sdp4
        # Compute the elements perturbed by the secular effects.
        n_k, e_k, i_k, Ω_k, ω_k, M_k = _dssec!(
            sgp4ds,
            nll₀,
            e₀,
            i₀,
            ω₀,
            Ω_k,
            ω_k,
            M_k,
            ∂ω,
            Δt
        )

        a_k  = (XKE / n_k)^(2 // 3) * (1 - C1 * Δt)^2
        e_k += -bstar * C4 * Δt
        M_k += (3 // 2) * nll₀ * C1 * Δt^2

    # Check if perigee is above 220 km.
    elseif algorithm === :sgp4

        sin_M₀, cos_M₀ = sincos(M₀)
        δω  = bstar * C3 * cos(ω₀) * Δt

        # TODO: sin(M_k) and cos(M_k) can be computed faster here.

        δM  = (e₀ > 1e-4) ?
            -(2 // 3) * QOMS2T * bstar * ξ^4 * AE / (e₀ * η) * (
                (1 + η * cos(M_k))^3 - (1 + η * cos_M₀)^3
            ) : T(0)
        M_k += +δω + δM
        ω_k += -δω - δM
        e_k  = e₀ - bstar * C4 * Δt - bstar * C5 * (sin(M_k) - sin_M₀)
        a_k  = all₀ * (@evalpoly(Δt, 1, -C1, -D2, -D3, -D4))^2
        IL   = M_k + ω_k + Ω_k + nll₀ * @evalpoly(
            Δt,
            0,
            0,
            (3 // 2) * C1,
            +(D2 + 2C1^2),
            +(3D3 + 12C1 * D2 + 10C1^3) / 4,
            +(3D4 + 12C1 * D3 + 6D2^2 + 30C1^2 * D2 + 15C1^4) / 5
        )

    elseif algorithm === :sgp4_lowper
        # If so, then
        #     1. Drop all terms after C1 in `a` and `IL`.
        #     2. Drop all terms involving C5.
        #     3. Drop δω.
        #     4. Drop δM.
        e_k = e₀ - bstar * C4 * Δt
        a_k = all₀ * (1 - C1 * Δt)^2
        IL  = M_k + ω_k + Ω_k + (3 // 2) * nll₀ * C1 * Δt^2
    else
        error("Unknown algorithm :$algorithm. Possible values are :sgp4, :sgp4_lowper, :sdp4.")
    end

    # TODO: Vallado's implementation [2] apply this normalization to the mean anomaly. It is
    # necessary to verify the reason for that.
    M_k_aux = M_k + ω_k + Ω_k
    Ω_k     = rem2pi(Ω_k, RoundToZero)
    ω_k     = rem2pi(ω_k, RoundToZero)
    M_k_aux = rem2pi(M_k_aux, RoundToZero)
    M_k     = rem2pi(M_k_aux - ω_k - Ω_k, RoundToZero)

    # == Lunar-Solar Periodics for Deep Space Orbits =======================================

    # This is only necessary if we are using SDP4 algorithm.
    if algorithm === :sdp4
        # Compute the elements perturbed by the Lunar-Solar periodics.
        e_k, i_k, Ω_k, ω_k, M_k = _dsper!(sgp4ds, e_k, i_k, Ω_k, ω_k, M_k, Δt)

        IL = M_k + ω_k + Ω_k

        # Make sure that the inclination is always positive.
        if i_k < 0
            i_k = -i_k
            Ω_k += T(π)
            ω_k -= T(π)
        end

        # The inclination was changed, hence some auxiliary variables must be
        # recomputed.
        sin_i_k, θ = sincos(i_k)
        θ²         = θ^2
    end

    # Vallado's code does not let the eccentricity to be smaller than 1e-6.
    #
    # TODO: Verify why this is necessary. I did not find any reason for that.
    e_k = max(e_k, T(1e-6))

    β = √(1 - e_k^2)

    # Compute the angular velocity [rad/min].
    n_k = XKE / √(a_k^3)

    # == Long-Period Periodic Term =========================================================

    sin_ω_k, cos_ω_k = sincos(ω_k)

    a_xN = e_k * cos_ω_k

    # TODO: Vallado's implementation of SGP4 uses another equation here.  However, both
    # produces the same result. Verify which one is better.
    a_yNL = A₃₀ * sin_i_k / (4k₂ * a_k * β^2)
    a_yN  = e_k * sin_ω_k + a_yNL
    IL_L  = (1 // 2) * a_yNL * a_xN * (3 + 5θ) / (1 + θ)
    IL_T  = IL + IL_L

    # == Solve Kepler's Equation for (E + ω) ===============================================

    U = mod(IL_T - Ω_k, T(2π))

    E_ω = U

    # Define the following variables that will be modified inside the loop so that we can
    # use them after the loop.
    sin_E_ω = T(0)
    cos_E_ω = T(0)

    for _ in 1:10
        sin_E_ω, cos_E_ω = sincos(E_ω)

        ΔE_ω = (U - a_yN * cos_E_ω + a_xN * sin_E_ω - E_ω) /
               (1 - a_yN * sin_E_ω - a_xN * cos_E_ω)

        # Vallado proposes to limit the maximum increment.
        abs(ΔE_ω) >= 0.95 && (ΔE_ω = sign(ΔE_ω) * T(0.95))

        E_ω += ΔE_ω

        # If the increment is less than a threshold, break the loop.
        #
        # Vallado proposes a threshold of 10^-12 instead of 10^-6.
        abs(ΔE_ω) < 1e-12 && break
    end

    # == Short-Term Periodic Terms =========================================================

    # Auxiliary variables.
    # NOTE: the sine and cosine of E + ω was already computed in the previous loop.
    e_cos_E = a_xN * cos_E_ω + a_yN * sin_E_ω
    e_sin_E = a_xN * sin_E_ω - a_yN * cos_E_ω
    e_L²    = a_xN^2 + a_yN^2
    p_L     = a_k * (1 - e_L²)
    p_L²    = p_L^2
    r       = a_k * (1 - e_cos_E)
    ṙ       = XKE * √a_k * e_sin_E / r
    rḟ      = XKE * √p_L / r
    auxsp   = e_sin_E / (1 + √(1 - e_L²))
    cos_u   = a_k / r * (cos_E_ω - a_xN + a_yN * auxsp)
    sin_u   = a_k / r * (sin_E_ω - a_yN - a_xN * auxsp)
    cos_2u  = 1 - 2sin_u^2
    sin_2u  = 2cos_u * sin_u
    u       = atan(sin_u, cos_u)

    # Short-term periodic terms.
    Δr  = +k₂ / (2p_L) * (1 - θ²) * cos_2u
    Δu  = -k₂ / (4p_L²) * (7θ² - 1) * sin_2u
    ΔΩ  = +3k₂ * θ / (2p_L²) * sin_2u
    Δi  = +3k₂ * θ / (2p_L²) * sin_i_k * cos_2u
    Δṙ  = -k₂ * n_k / p_L * (1 - θ²) * sin_2u
    Δrḟ = +k₂ * n_k / p_L * ((1 - θ²) * cos_2u - (3 // 2) * (1 - 3θ²))

    # The short-term periodics are added to give the osculating quantities.
    r_k  = r * (1 - (3 // 2) * k₂ * √(1 - e_L²) / p_L² * (3θ² - 1)) + Δr
    u_k  = u   + Δu
    Ω_k  = Ω_k + ΔΩ
    i_k  = i_k + Δi
    ṙ_k  = ṙ   + Δṙ
    rḟ_k = rḟ  + Δrḟ

    # Orientation vectors.
    sin_Ω_k, cos_Ω_k = sincos(Ω_k)
    sin_i_k, cos_i_k = sincos(i_k)
    sin_u_k, cos_u_k = sincos(u_k)

    M = SVector{3, T}(-sin_Ω_k * cos_i_k, +cos_Ω_k * cos_i_k, sin_i_k)
    N = SVector{3, T}(+cos_Ω_k,           +sin_Ω_k,           T(0))

    Uv = M * sin_u_k + N * cos_u_k
    Vv = M * cos_u_k - N * sin_u_k

    r_teme = r_k * Uv * R0
    v_teme = (ṙ_k * Uv + rḟ_k * Vv) * R0 / 60

    # Update the variables.
    sgp4d.Δt  = Δt
    sgp4d.a_k = a_k
    sgp4d.e_k = e_k
    sgp4d.i_k = i_k
    sgp4d.Ω_k = Ω_k
    sgp4d.ω_k = ω_k
    sgp4d.M_k = M_k
    sgp4d.n_k = n_k

    return r_teme, v_teme
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# == Deep Space Functions ==================================================================

"""
    _dsinit!(sgp4ds::Sgp4DeepSpace{T}, args...) where {Tepoch<:Number, T<:Number} -> Nothing

Initialize the deep space structure `sgp4ds` using the parameters in `args...`.

This function computes several parameters in `sgp4ds` that will be used when calling the
functions `_dsper!` and `_dssec!`.

# Arguments

- `sgp4ds::Sgp4DeepSpace`: Structure that will be initialized.
- `epoch::Tepoch`: Epoch of the initial orbit [Julian Day].
- `nll₀::T`: Initial mean motion [rad/min].
- `all₀::T`: Initial semi-major axis [ER].
- `e₀::T`: Initial eccentricity.
- `i₀::T`: Initial inclination [rad].
- `Ω₀::T`: Initial right ascension of the ascending node [rad].
- `ω₀::T`: Initial argument of perigee [rad].
- `M₀::T`: Initial mean anomaly [rad].
- `∂M::T`: Time-derivative of the mean motion [rad/min].
- `∂ω::T`: Time-derivative of the argument of perigee [rad/min].
- `∂Ω::T`: Time-derivative of the RAAN [rad/min].
"""
function _dsinit!(
    sgp4ds::Sgp4DeepSpace{ST},
    epoch::Tepoch,
    nll₀::NT,
    all₀::AT,
    e₀::ET,
    i₀::IT,
    Ω₀::OT,
    ω₀::WT,
    M₀::MT,
    ∂M::MT,
    ∂ω::WT,
    ∂Ω::OT
) where {Tepoch<:Number, NT<:Number, AT<:Number, ET<:Number, IT<:Number, OT<:Number, WT<:Number, MT<:Number, ST<:Number}

    T = ST

    # Unpack variables.
    atime  = sgp4ds.atime
    xli    = sgp4ds.xli
    xni    = sgp4ds.xni
    xnq    = sgp4ds.xnq
    xfact  = sgp4ds.xfact
    ssl    = sgp4ds.ssl
    ssg    = sgp4ds.ssg
    ssh    = sgp4ds.ssh
    sse    = sgp4ds.sse
    ssi    = sgp4ds.ssi
    xlamo  = sgp4ds.xlamo
    omegaq = sgp4ds.omegaq
    omgdt  = sgp4ds.omgdt
    gmst   = sgp4ds.gmst
    del1   = sgp4ds.del1
    del2   = sgp4ds.del2
    del3   = sgp4ds.del3
    fasx2  = sgp4ds.fasx2
    fasx4  = sgp4ds.fasx4
    fasx6  = sgp4ds.fasx6
    d2201  = sgp4ds.d2201
    d2211  = sgp4ds.d2211
    d3210  = sgp4ds.d3210
    d3222  = sgp4ds.d3222
    d4410  = sgp4ds.d4410
    d4422  = sgp4ds.d4422
    d5220  = sgp4ds.d5220
    d5232  = sgp4ds.d5232
    d5421  = sgp4ds.d5421
    d5433  = sgp4ds.d5433
    xnddt  = sgp4ds.xnddt
    xndot  = sgp4ds.xndot
    xldot  = sgp4ds.xldot
    zmos   = sgp4ds.zmos
    se2    = sgp4ds.se2
    se3    = sgp4ds.se3
    si2    = sgp4ds.si2
    si3    = sgp4ds.si3
    sl2    = sgp4ds.sl2
    sl3    = sgp4ds.sl3
    sl4    = sgp4ds.sl4
    sgh2   = sgp4ds.sgh2
    sgh3   = sgp4ds.sgh3
    sgh4   = sgp4ds.sgh4
    sh2    = sgp4ds.sh2
    sh3    = sgp4ds.sh3
    zmol   = sgp4ds.zmol
    ee2    = sgp4ds.ee2
    e3     = sgp4ds.e3
    xi2    = sgp4ds.xi2
    xi3    = sgp4ds.xi3
    xl2    = sgp4ds.xl2
    xl3    = sgp4ds.xl3
    xl4    = sgp4ds.xl4
    xgh2   = sgp4ds.xgh2
    xgh3   = sgp4ds.xgh3
    xgh4   = sgp4ds.xgh4
    xh2    = sgp4ds.xh2
    xh3    = sgp4ds.xh3
    pe     = sgp4ds.pe
    pinc   = sgp4ds.pinc
    pgh    = sgp4ds.pgh
    ph     = sgp4ds.ph
    pl     = sgp4ds.pl
    pgh0   = sgp4ds.pgh0
    ph0    = sgp4ds.ph0
    pe0    = sgp4ds.pe0
    pinc0  = sgp4ds.pinc0
    pl0    = sgp4ds.pl0
    isynfl = sgp4ds.isynfl
    iresfl = sgp4ds.iresfl
    ilsz   = sgp4ds.ilsz

    # == Constants =========================================================================

    ZNS    = T(1.19459E-5)
    C1SS   = T(2.9864797e-6)
    ZES    = T(0.01675)
    ZNL    = T(1.5835218e-4)
    ZEL    = T(0.05490)
    C1L    = T(4.7968065e-7)
    ZSINIS = T(0.39785416)
    ZCOSIS = T(0.91744867)
    ZCOSGS = T(0.1945905)
    ZSINGS = T(-0.98088458)
    Q22    = T(1.7891679e-6)
    Q31    = T(2.1460748e-6)
    Q33    = T(2.2123015e-7)
    G22    = T(5.7686396)
    G32    = T(0.95240898)
    G44    = T(1.8014998)
    G52    = T(1.0508330)
    G54    = T(4.4108898)
    ROOT22 = T(1.7891679e-6)
    ROOT32 = T(3.7393792e-7)
    ROOT44 = T(7.3636953e-9)
    ROOT52 = T(1.1428639e-7)
    ROOT54 = T(2.1765803e-9)
    THDT   = T(4.37526908801129966e-3)

    # == Auxiliary Variables ===============================================================

    e₀²        = e₀ * e₀
    sqrt_1_e₀² = √(1 - e₀²)
    inv_all₀   = 1 / all₀
    inv_nll₀   = 1 / nll₀
    se         = T(0)
    si         = T(0)
    sl         = T(0)
    sgh        = T(0)
    shdq       = T(0)

    sin_i₀, cos_i₀ = sincos(i₀)
    sin_Ω₀, cos_Ω₀ = sincos(Ω₀)
    sin_ω₀, cos_ω₀ = sincos(ω₀)

    sin_i₀² = sin_i₀ * sin_i₀
    cos_i₀² = cos_i₀ * cos_i₀
    xpidot  = ∂ω + ∂Ω

    # == Initial Configuration =============================================================

    # Drop terms if inclination is smaller than 3 deg.
    ishq = (i₀ >= 3π / 180) ? true : false

    # Do not let `sin_i₀` be 0.
    abs(sin_i₀) < 1e-12 && (sin_i₀ = sign(sin_i₀) * T(1e-12))

    # Compute the Greenwich Mean Sidereal Time at epoch.
    gmst = T(jd_to_gmst(epoch))

    # == Initialize Lunar Solar Terms ======================================================

    # `day` is the number of days since Jan 0, 1900 at 12h.
    day = T(epoch - (_JD_1900 - 1))

    xnodce = mod(T(4.5236020) - T(9.2422029e-4) * day, T(2π))

    stem, ctem = sincos(xnodce)

    zcosil = T(0.91375164) - T(0.03568096) * ctem
    zsinil = sqrt(1 - zcosil^2)
    zsinhl = T(0.089683511) * stem / zsinil
    zcoshl = sqrt(1 - zsinhl^2)
    gam    = T(5.8351514) + T(0.0019443680) * day
    zx     = T(0.39785416) * stem / zsinil
    zy     = zcoshl * ctem + T(0.91744867) * zsinhl * stem
    zx     = atan(zx, zy)
    zx     = gam + zx - xnodce

    zsingl, zcosgl = sincos(zx)

    zmol = mod(T(4.7199672) + T(0.22997150)  * day - gam, T(2π))
    zmos = mod(T(6.2565837) + T(0.017201977) * day,       T(2π))

    # == Do Solar Terms ====================================================================

    zcosg = ZCOSGS
    zsing = ZSINGS
    zcosi = ZCOSIS
    zsini = ZSINIS
    zcosh = cos_Ω₀
    zsinh = sin_Ω₀
    cc    = C1SS
    zn    = ZNS
    ze    = ZES
    zmo   = zmos

    for ls = 0:1
        a1  = +zcosg * zcosh + zsing * zcosi * zsinh
        a3  = -zsing * zcosh + zcosg * zcosi * zsinh
        a7  = -zcosg * zsinh + zsing * zcosi * zcosh
        a8  = +zsing * zsini
        a9  = +zsing * zsinh + zcosg * zcosi * zcosh
        a10 = +zcosg * zsini
        a2  = +cos_i₀ * a7  + sin_i₀ * a8
        a4  = +cos_i₀ * a9  + sin_i₀ * a10
        a5  = -sin_i₀ * a7  + cos_i₀ * a8
        a6  = -sin_i₀ * a9  + cos_i₀ * a10

        x1 = +a1 * cos_ω₀ + a2 * sin_ω₀
        x2 = +a3 * cos_ω₀ + a4 * sin_ω₀
        x3 = -a1 * sin_ω₀ + a2 * cos_ω₀
        x4 = -a3 * sin_ω₀ + a4 * cos_ω₀
        x5 = +a5 * sin_ω₀
        x6 = +a6 * sin_ω₀
        x7 = +a5 * cos_ω₀
        x8 = +a6 * cos_ω₀

        z31 = 12x1^2    - 3x3^2
        z32 = 24x1 * x2 - 6x3 * x4
        z33 = 12x2^2    - 3x4^2
        z1  = 3(   a1^2 + a2^2   ) + z31 * e₀²
        z2  = 6(a1 * a3 + a2 * a4) + z32 * e₀²
        z3  = 3(   a3^2 + a4^2   ) + z33 * e₀²
        z11 = -6a1 * a5 + e₀² * (-24x1 * x7 - 6x3 * x5)
        z12 = -6(a1 * a6 + a3 * a5) + e₀² * (-24(x2 * x7 + x1 * x8) - 6(x3 * x6 + x4 * x5))
        z13 = -6a3 * a6 + e₀² * (-24x2 * x8 - 6x4 * x6)
        z21 = +6a2 * a5 + e₀² * (+24x1 * x5 - 6x3 * x7)
        z22 = +6(a4 * a5 + a2 * a6) + e₀² * (24(x2 * x5 + x1 * x6) - 6(x4 * x7 + x3 * x8) )
        z23 = +6a4 * a6 + e₀² * (24x2 * x6 - 6x4 * x8)
        z1  = +2z1 + (1 - e₀²) * z31
        z2  = +2z2 + (1 - e₀²) * z32
        z3  = +2z3 + (1 - e₀²) * z33
        s3  = +cc * inv_nll₀
        s2  = -T(0.5) * s3 / sqrt_1_e₀²
        s4  = +s3 * sqrt_1_e₀²
        s1  = -15e₀ * s4
        s5  = +x1 * x3 + x2 * x4
        s6  = +x2 * x3 + x1 * x4
        s7  = +x2 * x4 - x1 * x3
        se  = +s1 * zn * s5
        si  = +s2 * zn * (z11 + z13)
        sl  = -zn * s3 * (z1 + z3 - 14 - 6e₀²)
        sgh = +s4 * zn * (z31 + z33 - 6)

        shdq = zero(T)

        if ishq
            sh   = -zn * s2 * (z21 + z23);
            shdq = sh / sin_i₀;
        end

        ee2  =  +2s1 * s6
        e3   =  +2s1 * s7
        xi2  =  +2s2 * z12
        xi3  =  +2s2 * (z13 - z11)
        xl2  =  -2s3 * z2
        xl3  =  -2s3 * (z3 - z1)
        xl4  =  -2s3 * (-21 - 9e₀²) * ze
        xgh2 =  +2s4 * z32
        xgh3 =  +2s4 * (z33 - z31)
        xgh4 = -18s4 * ze
        xh2  =  -2s2 * z22
        xh3  =  -2s2 * (z23 - z21)

        ls == 1 && break

        # == Do Lunar Terms ================================================================

        sse   = se
        ssi   = si
        ssl   = sl
        ssh   = shdq
        ssg   = sgh - cos_i₀ * ssh
        se2   = ee2
        si2   = xi2
        sl2   = xl2
        sgh2  = xgh2
        sh2   = xh2
        se3   = e3
        si3   = xi3
        sl3   = xl3
        sgh3  = xgh3
        sh3   = xh3
        sl4   = xl4
        sgh4  = xgh4
        zcosg = zcosgl
        zsing = zsingl
        zcosi = zcosil
        zsini = zsinil
        zcosh = cos_Ω₀ * zcoshl + sin_Ω₀ * zsinhl
        zsinh = sin_Ω₀ * zcoshl - cos_Ω₀ * zsinhl
        zn    = ZNL
        cc    = C1L
        ze    = ZEL
        zmo   = zmol
    end

    sse += se
    ssi += si
    ssl += sl
    ssg += sgh - cos_i₀ * shdq
    ssh += shdq

    if (nll₀ < T(0.0052359877)) && (nll₀ > T(0.0034906585))
        # == 24h Synchronous Resonance Terms Initialization ================================

        iresfl = true;
        isynfl = true;

        g200    = e₀² * (T(0.8125) * e₀² - T(2.5)) + 1
        g310    = 2e₀² + 1
        g300    = e₀² * (T(6.60937) * e₀² - 6) + 1
        f220    = T(0.75) * (cos_i₀ + 1)^2
        f311    = T(0.9375) * (3cos_i₀ + 1) * sin_i₀^2 - T(0.75) * (cos_i₀ + 1)
        f330    = T(1.875) * (cos_i₀ + 1)^3
        del1    = 3(nll₀^2 * inv_all₀^2)
        del2    = 2del1 * f220 * g200 * Q22
        del3    = 3del1 * f330 * g300 * Q33 * inv_all₀
        del1    =  del1 * f311 * g310 * Q31 * inv_all₀
        fasx2   = T(0.13130908)
        fasx4   = T(2.8843198)
        fasx6   = T(0.37448087)
        xlamo   = mod(M₀ + Ω₀ + ω₀ - gmst, T(2π))
        bfact   = ∂M + xpidot - THDT + ssl + ssg + ssh

    elseif (nll₀ >= T(0.00826)) && (nll₀ <= T(0.00924)) && (e₀ >= T(0.5))
        # == Geopotential Resonance Initialization for 12 Hour Orbits ======================

        iresfl = true
        isynfl = false

        g201 = -T(0.306) - T(0.44) * (e₀ - T(0.64))

        if e₀ <= 0.65
            g211 = @evalpoly(e₀, +T( 3.6160), -T( 13.2470), +T( 16.29000))
            g310 = @evalpoly(e₀, -T(19.3020), +T(117.3900), -T(228.4190 ), +T( 156.5910))
            g322 = @evalpoly(e₀, -T(18.9068), +T(109.7927), -T(214.6334 ), +T( 146.5816))
            g410 = @evalpoly(e₀, -T(41.1220), +T(242.6940), -T(471.0940 ), +T( 313.9530))
            g422 = @evalpoly(e₀, -T(146.407), +T(841.8800), -T(1629.014 ), +T(1083.435 ))
            g520 = @evalpoly(e₀, -T(532.114), +T(3017.977), -T(5740.032 ), +T(3708.276 ))
        else
            g211 = @evalpoly(e₀, -  T(72.099), +T(  331.8190), -T( 508.7380), +T(  266.7240))
            g310 = @evalpoly(e₀, - T(346.844), +T( 1582.851 ), -T( 2415.925), +T( 1246.113 ))
            g322 = @evalpoly(e₀, - T(342.585), +T( 1554.908 ), -T( 2366.899), +T( 1215.972 ))
            g410 = @evalpoly(e₀, -T(1052.797), +T( 4758.686 ), -T( 7193.992), +T( 3651.957 ))
            g422 = @evalpoly(e₀, -T(3581.690), +T(16178.11  ), -T(24462.77 ), +T(12422.52  ))

            if e₀ <= T(0.715)
                g520 = @evalpoly(e₀, +T(1464.74), -T(4664.75), +T(3763.64))
            else
                g520 = @evalpoly(e₀, -T(5149.66), +T(29936.92), -T(54087.36), +T(31324.56))
            end
        end

        if e₀ < T(0.7)
            g533 = @evalpoly(e₀, -T(919.22770), +T(4988.6100), -T(9064.7700), +T(5542.210))
            g521 = @evalpoly(e₀, -T(822.71072), +T(4568.6173), -T(8491.4146), +T(5337.524))
            g532 = @evalpoly(e₀, -T(853.66600), +T(4690.2500), -T(8624.7700), +T(5341.400))
        else
            g533 = @evalpoly(e₀, -T(37995.780), +T(161616.52), -T(229838.20), +T(109377.94))
            g521 = @evalpoly(e₀, -T(51752.104), +T(218913.95), -T(309468.16), +T(146349.42))
            g532 = @evalpoly(e₀, -T(40023.880), +T(170470.89), -T(242699.48), +T(115605.82))
        end

        f220 = +T(0.75)  * (1 + 2cos_i₀ + cos_i₀²)
        f221 = +T(1.5)   * sin_i₀²
        f321 = +T(1.875) * sin_i₀ * (1 - 2cos_i₀ - 3cos_i₀²)
        f322 = -T(1.875) * sin_i₀ * (1 + 2cos_i₀ - 3cos_i₀²)
        f441 = +35sin_i₀² * f220
        f442 = +T(39.375) * sin_i₀²^2
        f522 = +T(9.84375) * sin_i₀ * (
            sin_i₀² * (+1 - 2cos_i₀ - 5cos_i₀²) +
            T(0.33333333) * (-2 + 4cos_i₀ + 6cos_i₀²)
        )
        f523 = sin_i₀ * (
            T(4.92187512) * sin_i₀² * (-2 - 4cos_i₀ + 10cos_i₀²) +
            T(6.56250012) * (+1 + 2cos_i₀ -  3cos_i₀²)
        )
        f542 = T(29.53125) * sin_i₀ * (
            +2 - 8cos_i₀ + cos_i₀² * (-12 + 8cos_i₀ + 10cos_i₀²)
        )
        f543 = T(29.53125) * sin_i₀ * (
            -2 - 8cos_i₀ + cos_i₀² * (+12 + 8cos_i₀ - 10cos_i₀²)
        )

        temp1   = 3 * (nll₀ * inv_all₀)^2
        temp0   = temp1 * ROOT22
        d2201   = temp0 * f220 * g201
        d2211   = temp0 * f221 * g211
        temp1  *= inv_all₀
        temp0   = temp1 * ROOT32
        d3210   = temp0 * f321 * g310
        d3222   = temp0 * f322 * g322
        temp1  *= inv_all₀
        temp0   = 2temp1 * ROOT44
        d4410   = temp0 * f441 * g410
        d4422   = temp0 * f442 * g422
        temp1  *= inv_all₀
        temp0   = temp1 * ROOT52
        d5220   = temp0 * f522 * g520
        d5232   = temp0 * f523 * g532
        temp0   = 2temp1 * ROOT54
        d5421   = temp0 * f542 * g521
        d5433   = temp0 * f543 * g533
        xlamo   = mod(M₀ + 2Ω₀ - 2gmst, T(2π))
        bfact   = ∂M + 2∂Ω - 2THDT + ssl + 2ssh
    else
        # == Non Resonant Orbits ===========================================================

        iresfl = false
        isynfl = false
    end

    if iresfl
        # == Initialize the Integrator =====================================================

        xfact = bfact - nll₀
        xli   = xlamo
        atime = T(0)

        # TODO: Check if this variable can be removed from Sgp4DeepSpace.
        xni   = nll₀

        # == Compute the "dot" Terms =======================================================

        if isynfl
            sin_1, cos_1 = sincos(  (xli - fasx2) )
            sin_2, cos_2 = sincos( 2(xli - fasx4) )
            sin_3, cos_3 = sincos( 3(xli - fasx6) )

            xndot = del1 * sin_1 +  del2 * sin_2 +  del3 * sin_3
            xnddt = del1 * cos_1 + 2del2 * cos_2 + 3del3 * cos_3
        else
            ω = ω₀ + ∂ω * atime

            sin_1,  cos_1  = sincos(2ω + xli  - G22)
            sin_2,  cos_2  = sincos(   + xli  - G22)
            sin_3,  cos_3  = sincos(+ω + xli  - G32)
            sin_4,  cos_4  = sincos(-ω + xli  - G32)
            sin_5,  cos_5  = sincos(+ω + xli  - G52)
            sin_6,  cos_6  = sincos(-ω + xli  - G52)
            sin_7,  cos_7  = sincos(2ω + 2xli - G44)
            sin_8,  cos_8  = sincos(     2xli - G44)
            sin_9,  cos_9  = sincos(+ω + 2xli - G54)
            sin_10, cos_10 = sincos(-ω + 2xli - G54)

            xndot = d2201 * sin_1 + d2211 * sin_2 + d3210 * sin_3 +
                    d3222 * sin_4 + d5220 * sin_5 + d5232 * sin_6 +
                    d4410 * sin_7 + d4422 * sin_8 + d5421 * sin_9 +
                    d5433 * sin_10

            xnddt =   d2201 * cos_1 + d2211 * cos_2 + d3210 * cos_3 +
                      d3222 * cos_4 + d5220 * cos_5 + d5232 * cos_6 +
                    2(d4410 * cos_7 + d4422 * cos_8 + d5421 * cos_9 +
                      d5433 * cos_10)
        end

        xldot  = xni + xfact
        xnddt *= xldot
    end

    # Set up for original mode (LS terms at epoch non-zero).
    pgh0 = ph0 = pe0 = pinc0 = pl0 = T(0)

    # Pack variables.
    sgp4ds.atime  = atime
    sgp4ds.xli    = xli
    sgp4ds.xni    = xni
    sgp4ds.xnq    = xnq
    sgp4ds.xfact  = xfact
    sgp4ds.ssl    = ssl
    sgp4ds.ssg    = ssg
    sgp4ds.ssh    = ssh
    sgp4ds.sse    = sse
    sgp4ds.ssi    = ssi
    sgp4ds.xlamo  = xlamo
    sgp4ds.omegaq = omegaq
    sgp4ds.omgdt  = omgdt
    sgp4ds.gmst   = gmst
    sgp4ds.del1   = del1
    sgp4ds.del2   = del2
    sgp4ds.del3   = del3
    sgp4ds.fasx2  = fasx2
    sgp4ds.fasx4  = fasx4
    sgp4ds.fasx6  = fasx6
    sgp4ds.d2201  = d2201
    sgp4ds.d2211  = d2211
    sgp4ds.d3210  = d3210
    sgp4ds.d3222  = d3222
    sgp4ds.d4410  = d4410
    sgp4ds.d4422  = d4422
    sgp4ds.d5220  = d5220
    sgp4ds.d5232  = d5232
    sgp4ds.d5421  = d5421
    sgp4ds.d5433  = d5433
    sgp4ds.xnddt  = xnddt
    sgp4ds.xndot  = xndot
    sgp4ds.xldot  = xldot
    sgp4ds.zmos   = zmos
    sgp4ds.se2    = se2
    sgp4ds.se3    = se3
    sgp4ds.si2    = si2
    sgp4ds.si3    = si3
    sgp4ds.sl2    = sl2
    sgp4ds.sl3    = sl3
    sgp4ds.sl4    = sl4
    sgp4ds.sgh2   = sgh2
    sgp4ds.sgh3   = sgh3
    sgp4ds.sgh4   = sgh4
    sgp4ds.sh2    = sh2
    sgp4ds.sh3    = sh3
    sgp4ds.zmol   = zmol
    sgp4ds.ee2    = ee2
    sgp4ds.e3     = e3
    sgp4ds.xi2    = xi2
    sgp4ds.xi3    = xi3
    sgp4ds.xl2    = xl2
    sgp4ds.xl3    = xl3
    sgp4ds.xl4    = xl4
    sgp4ds.xgh2   = xgh2
    sgp4ds.xgh3   = xgh3
    sgp4ds.xgh4   = xgh4
    sgp4ds.xh2    = xh2
    sgp4ds.xh3    = xh3
    sgp4ds.pe     = pe
    sgp4ds.pinc   = pinc
    sgp4ds.pgh    = pgh
    sgp4ds.ph     = ph
    sgp4ds.pl     = pl
    sgp4ds.pgh0   = pgh0
    sgp4ds.ph0    = ph0
    sgp4ds.pe0    = pe0
    sgp4ds.pinc0  = pinc0
    sgp4ds.pl0    = pl0
    sgp4ds.isynfl = isynfl
    sgp4ds.iresfl = iresfl
    sgp4ds.ilsz   = ilsz

    return nothing
end

"""
    _dssec!(sgp4ds::Sgp4DeepSpace{T}, nll₀::T, e₀::T, i₀::T, ω₀::T, Ω_k::T, ω_k::T, M_k::T, ∂ω::T, Δt::Number) where T<:Number

Compute the secular effects.

!!! note

    The internal values in `sgp4ds` will be modified.

# Arguments

- `sgp4ds::Sgp4DeepSpace`: Deep space structure (see [`Sgp4DeepSpace`](@ref)).
- `nll₀::Number`: Initial mean motion [rad/min].
- `e₀::Number`: Initial eccentricity.
- `i₀::Number`: Initial inclination [rad].
- `ω₀::Number`: Initial argument of perigee [rad].
- `Ω_k::Number`: Current right ascension of the ascending node [rad].
- `ω_k::Number`: Current argument of perigee [rad].
- `M_k::Number`: Current mean anomaly [rad].
- `∂ω::Number`: Time-derivative of the argument of perigee [rad/min].
- `Δt::Number`: Time interval since the epoch [min].

# Returns

The following elements perturbed by the secular effects:

- `T`: Mean motion [rad/min].
- `T`: Eccentricity.
- `T`: Inclination [rad].
- `T`: Right ascension of the ascending node [rad].
- `T`: Argument of perigee [rad].
- `T`: Mean anomaly [rad].
"""
function _dssec!(
    sgp4ds::Sgp4DeepSpace{T},
    nll₀::T,
    e₀::T,
    i₀::T,
    ω₀::T,
    Ω_k::T,
    ω_k::T,
    M_k::T,
    ∂ω::T,
    Δt::Number
) where T<:Number

    # Unpack variables.
    atime  = sgp4ds.atime
    xli    = sgp4ds.xli
    xni    = sgp4ds.xni
    xfact  = sgp4ds.xfact
    ssl    = sgp4ds.ssl
    ssg    = sgp4ds.ssg
    ssh    = sgp4ds.ssh
    sse    = sgp4ds.sse
    ssi    = sgp4ds.ssi
    xlamo  = sgp4ds.xlamo
    gmst   = sgp4ds.gmst
    del1   = sgp4ds.del1
    del2   = sgp4ds.del2
    del3   = sgp4ds.del3
    fasx2  = sgp4ds.fasx2
    fasx4  = sgp4ds.fasx4
    fasx6  = sgp4ds.fasx6
    d2201  = sgp4ds.d2201
    d2211  = sgp4ds.d2211
    d3210  = sgp4ds.d3210
    d3222  = sgp4ds.d3222
    d4410  = sgp4ds.d4410
    d4422  = sgp4ds.d4422
    d5220  = sgp4ds.d5220
    d5232  = sgp4ds.d5232
    d5421  = sgp4ds.d5421
    d5433  = sgp4ds.d5433
    xnddt  = sgp4ds.xnddt
    xndot  = sgp4ds.xndot
    xldot  = sgp4ds.xldot
    iresfl = sgp4ds.iresfl
    isynfl = sgp4ds.isynfl

    # == Constants =========================================================================

    STEP = T(720.0)
    G22  = T(5.7686396)
    G32  = T(0.95240898)
    G44  = T(1.8014998)
    G52  = T(1.0508330)
    G54  = T(4.4108898)
    THDT = T(4.37526908801129966e-3)

    # == Initialization ====================================================================

    M_sec = M_k + ssl * Δt
    e_sec = e₀  + sse * Δt
    i_sec = i₀  + ssi * Δt
    Ω_sec = Ω_k + ssh * Δt
    ω_sec = ω_k + ssg * Δt

    # TODO: Verify what this variable means. This is found in `dspace.m` of Vallado's
    # implementation [2].
    θ = mod(gmst + THDT * Δt, T(2π))

    # If the orbit is not resonant, then nothing more should be computed.
    !iresfl && return nll₀, e_sec, i_sec, Ω_sec, ω_sec, M_sec

    # == Update Resonances using Numerical (Euler-Maclaurin) Integration ===================

    # -- Epoch restart ---------------------------------------------------------------------

    # This verification is different between Vallado's [2] and [3]. We will use [2] since it
    # seems more recent.
    if  (atime == 0) || (Δt * atime <= 0) || (abs(Δt) < abs(atime))
        atime = T(0)
        xni   = nll₀
        xli   = xlamo
    end

    # -- Integration -----------------------------------------------------------------------

    ft = Δt - atime

    # In [3], the integration process is performed only if `ft` is larger than `STEP`.
    # However, Vallado's implementation [2] does not verify this and the integration is
    # performed every time. This behavior was chose because it seems that [3] is a more
    # recent version of the algorithm.

    # Check integration direction.
    delt = (Δt >= atime) ? STEP : -STEP

    # Perform the integration with step `delt` until the difference between the time `Δt`
    # and `atime` is less then `STEP`.
    while true
        # Compute the dot terms.
        if isynfl

            sin_1, cos_1 = sincos( (xli - fasx2))
            sin_2, cos_2 = sincos(2(xli - fasx4))
            sin_3, cos_3 = sincos(3(xli - fasx6))

            xndot = del1 * sin_1 +  del2 * sin_2 +  del3 * sin_3
            xnddt = del1 * cos_1 + 2del2 * cos_2 + 3del3 * cos_3
        else
            ω = ω₀ + ∂ω * atime

            sin_1,  cos_1  = sincos(2ω + xli  - G22)
            sin_2,  cos_2  = sincos(   + xli  - G22)
            sin_3,  cos_3  = sincos(+ω + xli  - G32)
            sin_4,  cos_4  = sincos(-ω + xli  - G32)
            sin_5,  cos_5  = sincos(+ω + xli  - G52)
            sin_6,  cos_6  = sincos(-ω + xli  - G52)
            sin_7,  cos_7  = sincos(2ω + 2xli - G44)
            sin_8,  cos_8  = sincos(     2xli - G44)
            sin_9,  cos_9  = sincos(+ω + 2xli - G54)
            sin_10, cos_10 = sincos(-ω + 2xli - G54)

            xndot = d2201 * sin_1 + d2211 * sin_2 + d3210 * sin_3 +
                    d3222 * sin_4 + d5220 * sin_5 + d5232 * sin_6 +
                    d4410 * sin_7 + d4422 * sin_8 + d5421 * sin_9 +
                    d5433 * sin_10

            xnddt =   d2201 * cos_1 + d2211 * cos_2 + d3210 * cos_3 +
                      d3222 * cos_4 + d5220 * cos_5 + d5232 * cos_6 +
                    2(d4410 * cos_7 + d4422 * cos_8 + d5421 * cos_9 +
                      d5433 * cos_10)
        end

        xldot  = xni + xfact
        xnddt *= xldot

        ft = Δt - atime
        (abs(ft) < STEP) && break

        # In Vallado's implementation [2], this is in the final of the loop instead of at
        # the beginning.
        xli   += delt * (xldot + delt * xndot / 2)
        xni   += delt * (xndot + delt * xnddt / 2)
        atime += delt
    end

    xl    = xli + ft * (xldot + ft * xndot / 2)
    n_sec = xni + ft * (xndot + ft * xnddt / 2)
    M_sec = !isynfl ? xl - 2Ω_sec + 2θ : xl - Ω_sec - ω_sec + θ

    # Pack variables.
    sgp4ds.atime = atime
    sgp4ds.xni   = xni
    sgp4ds.xli   = xli
    sgp4ds.xnddt = xnddt
    sgp4ds.xndot = xndot
    sgp4ds.xldot = xldot

    return n_sec, e_sec, i_sec, Ω_sec, ω_sec, M_sec
end

"""
    _dsper!(sgp4ds::Sgp4DeepSpace{T}, e_k::T, i_k::T, Ω_k::T, ω_k::T, M_k::T, Δt:Number) where T<:Number

Compute the effects caused by Lunar-Solar periodics.

!!! note

    The internal values in `sgp4ds` will be modified.

# Arguments

- `sgp4ds::Sgp4DeepSpace`: Deep space structure (see [`Sgp4DeepSpace`](@ref)).
- `e_k::Number`: Current eccentricity.
- `i_k::Number`: Current inclination [rad].
- `Ω_k::Number`: Current right ascension of the ascending node [rad].
- `ω_k::Number`: Current argument of perigee [rad].
- `M_k::Number`: Current mean anomaly [rad].
- `Δt::Number`: Time interval since the epoch [min].

# Returns

The following elements perturbed by lunar-solar periodics.

- `T`: Eccentricity.
- `T`: Inclination [rad].
- `T`: Right ascension of the ascending node [rad].
- `T`: Argument of perigee [rad].
- `T`: Mean anomaly [rad].
"""
function _dsper!(
    sgp4ds::Sgp4DeepSpace{T},
    e_k::T,
    i_k::T,
    Ω_k::T,
    ω_k::T,
    M_k::T,
    Δt::Number
) where T<:Number

    # Unpack variables.
    zmos = sgp4ds.zmos
    se2  = sgp4ds.se2
    se3  = sgp4ds.se3
    si2  = sgp4ds.si2
    si3  = sgp4ds.si3
    sl2  = sgp4ds.sl2
    sl3  = sgp4ds.sl3
    sl4  = sgp4ds.sl4
    sgh2 = sgp4ds.sgh2
    sgh3 = sgp4ds.sgh3
    sgh4 = sgp4ds.sgh4
    sh2  = sgp4ds.sh2
    sh3  = sgp4ds.sh3
    zmol = sgp4ds.zmol
    ee2  = sgp4ds.ee2
    e3   = sgp4ds.e3
    xi2  = sgp4ds.xi2
    xi3  = sgp4ds.xi3
    xl2  = sgp4ds.xl2
    xl3  = sgp4ds.xl3
    xl4  = sgp4ds.xl4
    xgh2 = sgp4ds.xgh2
    xgh3 = sgp4ds.xgh3
    xgh4 = sgp4ds.xgh4
    xh2  = sgp4ds.xh2
    xh3  = sgp4ds.xh3

    # == Constants =========================================================================

    ZNS  = T(1.19459E-5)
    ZES  = T(0.01675)
    ZNL  = T(1.5835218e-4)
    ZEL  = T(0.05490)

    # == Update Solar Terms ================================================================

    zm = zmos +  ZNS * Δt
    zf = zm   + 2ZES * sin(zm)

    sinzf, coszf = sincos(zf)

    f2   = +sinzf * sinzf / 2 - (1 // 4)
    f3   = -sinzf * coszf / 2
    ses  = se2 * f2 + se3 * f3
    sis  = si2 * f2 + si3 * f3
    sls  = sl2 * f2 + sl3 * f3 + sl4 * sinzf
    sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
    shs  = sh2  * f2 + sh3  * f3

    # == Update Lunar Terms ================================================================

    zm    = zmol +  ZNL * Δt
    zf    = zm   + 2ZEL * sin(zm)

    sinzf, coszf = sincos(zf)

    f2   = +sinzf * sinzf / 2 - (1 // 4)
    f3   = -sinzf * coszf / 2
    sel  = ee2 * f2 + e3 * f3
    sil  = xi2 * f2 + xi3 * f3
    sll  = xl2 * f2 + xl3 * f3 + xl4 * sinzf
    sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
    shl  = xh2  * f2 + xh3  * f3

    # == Save Computed Values ==============================================================

    pgh  = sghs + sghl
    ph   = shs  + shl
    pe   = ses  + sel
    pinc = sis  + sil
    pl   = sls  + sll

    # Update inclination and eccentricity.
    e_per = e_k + pe
    i_per = i_k + pinc

    sinis, cosis = sincos(i_per)

    # The original algorithm considered the original inclination to select the Lyddane
    # Lunar-Solar perturbations algorithm. However, Vallado's implementation [2] test the
    # perturbed inclination to select this. It is mentioned that this is the behavior
    # selected in GSFC source code.
    if i_per >= T(0.2)
        tmp_ph = ph / sinis;
        ω_per  = ω_k + pgh - cosis * tmp_ph;
        Ω_per  = Ω_k + tmp_ph;
        M_per  = M_k + pl;
    else
        sinok = sin(Ω_k)
        cosok = cos(Ω_k)

        #                     |----------    dalf     ----------|
        alfdp = sinis * sinok + ph * cosok + pinc * cosis * sinok
        #                     |----------    dbet     ----------|
        betdp = sinis * cosok - ph * sinok + pinc * cosis * cosok

        # For the following computation, in which `Ω_per` is used without a trigonometric
        # function, it is advisable to make sure that it stays in the interval [0, 2π].
        Ω_per = mod(Ω_k, T(2π))

        #                                 |----------    dls    ----------|
        xls   = M_k + ω_k + cosis * Ω_per + pl + pgh - pinc * Ω_per * sinis
        Ω_aux = Ω_per
        Ω_per = mod(atan(alfdp, betdp), T(2π))

        if abs(Ω_aux - Ω_per) > π
            Ω_per = (Ω_per < Ω_aux) ? Ω_per + T(2π) : Ω_per - T(2π)
        end

        M_per = M_k + pl;
        ω_per = xls - M_per - cosis * Ω_per
    end

    # Pack variables.
    sgp4ds.pgh  = pgh
    sgp4ds.ph   = ph
    sgp4ds.pe   = pe
    sgp4ds.pinc = pinc
    sgp4ds.pl   = pl

    return e_per, i_per, Ω_per, ω_per, M_per
end
