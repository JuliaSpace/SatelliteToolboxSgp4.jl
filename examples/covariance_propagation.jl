# ==========================================================================================
#              Covariance Propagation via SGP4 + ForwardDiff Jacobian
# ==========================================================================================
#
# Demonstrates linear covariance propagation through the SGP4 orbit propagator using
# ForwardDiff to compute the state transition Jacobian.
#
# Given an initial covariance P₀ in mean-element space, the propagated covariance at time
# Δt is:
#
#   P(Δt) = J(Δt) P₀ J(Δt)ᵀ
#
# where J = ∂[r; v]/∂x₀ is the 6×N Jacobian of the SGP4 mapping from mean elements to
# the Cartesian TEME state.
#
# ==========================================================================================

using SatelliteToolboxSgp4
using DiffResults
using ForwardDiff
using LinearAlgebra
using StaticArrays
using Printf

# ------------------------------------------------------------------------------------------
# TLE
# ------------------------------------------------------------------------------------------

tle_amazonia = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
    """

# ------------------------------------------------------------------------------------------
# Extract mean elements from the TLE into a state vector x₀
# ------------------------------------------------------------------------------------------

d2r   = π / 180
epoch = Float64(tle_epoch(tle_amazonia))
n₀    = Float64(tle_amazonia.mean_motion) * 2π / 1440   # [rad/min]
e₀    = Float64(tle_amazonia.eccentricity)
i₀    = Float64(tle_amazonia.inclination) * d2r          # [rad]
Ω₀    = Float64(tle_amazonia.raan) * d2r                 # [rad]
ω₀    = Float64(tle_amazonia.argument_of_perigee) * d2r  # [rad]
M₀    = Float64(tle_amazonia.mean_anomaly) * d2r         # [rad]
bstar = Float64(tle_amazonia.bstar)

# State vector: [n₀, e₀, i₀, Ω₀, ω₀, M₀, B*]
x₀ = [n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar]

# ------------------------------------------------------------------------------------------
# Define the SGP4 map: x₀ → [r; v] at time Δt
# ------------------------------------------------------------------------------------------

function sgp4_map(x, Δt, epoch_jd)
    r, v, _ = sgp4(Δt, epoch_jd, x[1], x[2], x[3], x[4], x[5], x[6], x[7])
    return vcat(SVector{3}(r), SVector{3}(v))
end

"""
    sgp4_value_and_jacobian(x₀, Δt, epoch_jd) -> (y, J)

Compute both the TEME state vector `y = [r; v]` and the Jacobian `J = ∂y/∂x₀` in a single
ForwardDiff pass using `DiffResults`.
"""
function sgp4_value_and_jacobian(x₀, Δt, epoch_jd)
    f = x -> sgp4_map(x, Δt, epoch_jd)
    result = DiffResults.JacobianResult(zeros(6), x₀)
    ForwardDiff.jacobian!(result, f, x₀)
    return DiffResults.value(result), DiffResults.jacobian(result)
end

# ------------------------------------------------------------------------------------------
# Fabricate an initial covariance P₀ in mean-element space
# ------------------------------------------------------------------------------------------
#
# In practice, P₀ would come from an orbit determination solution (e.g. a least-squares
# fit or Kalman filter). Here we use representative diagonal uncertainties for
# demonstration.
#
#   σ_n₀     = 1e-8  rad/min    (mean motion)
#   σ_e₀     = 1e-6             (eccentricity)
#   σ_i₀     = 1e-5  rad        (inclination)
#   σ_Ω₀     = 1e-4  rad        (RAAN)
#   σ_ω₀     = 1e-4  rad        (argument of perigee)
#   σ_M₀     = 1e-4  rad        (mean anomaly)
#   σ_bstar  = 1e-5             (B* drag term)

σ = [1e-8, 1e-6, 1e-5, 1e-4, 1e-4, 1e-4, 1e-5]
P₀ = Diagonal(σ .^ 2)

# ------------------------------------------------------------------------------------------
# Propagate covariance at several times
# ------------------------------------------------------------------------------------------

println("=" ^ 80)
println("Covariance Propagation via SGP4 + ForwardDiff Jacobian")
println("=" ^ 80)
println()
@printf("  Satellite: %s\n", tle_amazonia.name)
@printf("  Epoch:     JD %.8f\n", epoch)
@printf("  State:     x₀ = [n₀, e₀, i₀, Ω₀, ω₀, M₀, B*]\n")
println()

println("── Propagated 1σ uncertainties in TEME ────────────────────────────────")
println()
@printf(
    "  %8s  %12s  %12s  %12s  %12s  %12s  %12s\n",
    "Δt [min]",
    "σ_x [km]",
    "σ_y [km]",
    "σ_z [km]",
    "σ_vx [km/s]",
    "σ_vy [km/s]",
    "σ_vz [km/s]"
)
@printf(
    "  %8s  %12s  %12s  %12s  %12s  %12s  %12s\n",
    "--------",
    "------------",
    "------------",
    "------------",
    "------------",
    "------------",
    "------------"
)

for Δt in [0.0, 10.0, 30.0, 60.0, 120.0, 360.0, 720.0, 1440.0]
    _, J = sgp4_value_and_jacobian(x₀, Δt, epoch)

    P_t = J * P₀ * J'

    σ_r = sqrt.(diag(P_t)[1:3])
    σ_v = sqrt.(diag(P_t)[4:6])

    @printf(
        "  %8.1f  %12.6f  %12.6f  %12.6f  %12.6e  %12.6e  %12.6e\n",
        Δt,
        σ_r[1],
        σ_r[2],
        σ_r[3],
        σ_v[1],
        σ_v[2],
        σ_v[3]
    )
end

# ------------------------------------------------------------------------------------------
# Detailed view at Δt = 60 min
# ------------------------------------------------------------------------------------------

Δt_detail = 60.0

println()
println("── Detailed view at Δt = $Δt_detail min ──────────────────────────────")
println()

y, J = sgp4_value_and_jacobian(x₀, Δt_detail, epoch)

@printf("  Position TEME:  [%+14.6f, %+14.6f, %+14.6f] km\n", y[1], y[2], y[3])
@printf("  Velocity TEME:  [%+14.8f, %+14.8f, %+14.8f] km/s\n", y[4], y[5], y[6])
println()

println("  Jacobian J = ∂[r;v]/∂x₀  (6×7):")
for i in 1:6
    @printf("    [")
    for j in 1:7
        @printf(" %+12.4e", J[i, j])
        j < 7 && @printf(",")
    end
    @printf(" ]\n")
end

P_t = J * P₀ * J'

println()
println("  Propagated covariance P(Δt) (6×6):")
for i in 1:6
    @printf("    [")
    for j in 1:6
        @printf(" %+12.4e", P_t[i, j])
        j < 6 && @printf(",")
    end
    @printf(" ]\n")
end

σ_pos = sqrt(tr(P_t[1:3, 1:3]))
σ_vel = sqrt(tr(P_t[4:6, 4:6]))
println()
@printf("  RSS position uncertainty: %.6f km\n", σ_pos)
@printf("  RSS velocity uncertainty: %.6e km/s\n", σ_vel)

println()
println("=" ^ 80)
