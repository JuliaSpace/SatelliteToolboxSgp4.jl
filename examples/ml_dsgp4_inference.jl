# ==========================================================================================
#                    ML-∂SGP4 Inference Example
# ==========================================================================================
#
# Shows how to load a pre-trained ML-∂SGP4 model and run inference. Mirrors the standard
# SGP4 workflow:
#
#   SGP4:       sgp4d   = sgp4_init(tle; sgp4c = sgp4c_wgs84)
#               r, v    = sgp4!(sgp4d, Δt)
#
#   ML-∂SGP4:   mlsgp4d = ml_dsgp4_init(tle; model = model)
#               r, v    = ml_dsgp4!(mlsgp4d, Δt)
#
# Requires a saved model file (produced by ml_dsgp4_training.jl or your own training run).
#
# ==========================================================================================

using SatelliteToolboxSgp4
using Lux, Zygote, Optimisers
using StaticArrays
using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------------
# TLE to propagate
# ------------------------------------------------------------------------------------------

tle_amazonia = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
    """

# ------------------------------------------------------------------------------------------
# Load model (≈ Sgp4Constants), then bind to TLE (≈ sgp4_init)
# ------------------------------------------------------------------------------------------

model_path = joinpath(@__DIR__, "ml_dsgp4_model.bin")

if !isfile(model_path)
    error("No saved model found at $(model_path). Run ml_dsgp4_training.jl first.")
end

model   = ml_dsgp4_load(model_path)
mlsgp4d = ml_dsgp4_init(tle_amazonia; model = model)

println("=" ^ 72)
println("ML-∂SGP4 Inference Example")
println("=" ^ 72)
println()
@printf("  Loaded model from: %s\n", model_path)
@printf("  Hidden layers:     %s\n", string(model.config.hidden_layers))
println()

# ------------------------------------------------------------------------------------------
# In-place propagation:  ml_dsgp4!(mlsgp4d, Δt)  — mirrors sgp4!(sgp4d, Δt)
# ------------------------------------------------------------------------------------------

println("── In-place propagation (ml_dsgp4!) ──────────────────────────────────")
println()

@printf("  %8s  %14s  %14s  %14s  %14s  %14s  %14s\n",
    "Δt [min]", "x [km]", "y [km]", "z [km]", "vx [km/s]", "vy [km/s]", "vz [km/s]")
@printf("  %8s  %14s  %14s  %14s  %14s  %14s  %14s\n",
    "--------", "--------------", "--------------", "--------------",
    "--------------", "--------------", "--------------")

for Δt in 0.0:30.0:360.0
    r, v = ml_dsgp4!(mlsgp4d, Δt)
    @printf("  %8.1f  %14.4f  %14.4f  %14.4f  %14.6f  %14.6f  %14.6f\n",
        Δt, r[1], r[2], r[3], v[1], v[2], v[3])
end

# ------------------------------------------------------------------------------------------
# Convenience API:  ml_dsgp4(Δt, tle; model=...)  — mirrors sgp4(Δt, tle; sgp4c=...)
# ------------------------------------------------------------------------------------------

println()
println("── Convenience propagation (ml_dsgp4) ────────────────────────────────")
println()

r, v, _ = ml_dsgp4(120.0, tle_amazonia; model = model)
@printf("  ml_dsgp4(120.0, tle; model = model)\n")
@printf("    r = [%+.6f, %+.6f, %+.6f] km\n", r[1], r[2], r[3])
@printf("    v = [%+.8f, %+.8f, %+.8f] km/s\n", v[1], v[2], v[3])

# ------------------------------------------------------------------------------------------
# Side-by-side comparison with plain SGP4
# ------------------------------------------------------------------------------------------

println()
println("── ML-∂SGP4 vs plain SGP4 ───────────────────────────────────────────")
println()

sgp4d = sgp4_init(tle_amazonia)

@printf("  %8s  %16s  %16s\n", "Δt [min]", "|Δr| [km]", "|Δv| [km/s]")
@printf("  %8s  %16s  %16s\n", "--------", "----------------", "----------------")

for Δt in 0.0:60.0:720.0
    r_ml, v_ml = ml_dsgp4!(mlsgp4d, Δt)
    r_sg, v_sg = sgp4!(sgp4d, Δt)
    Δr = norm(collect(r_ml) .- collect(r_sg))
    Δv = norm(collect(v_ml) .- collect(v_sg))
    @printf("  %8.1f  %16.6f  %16.6e\n", Δt, Δr, Δv)
end

println()
println("=" ^ 72)
