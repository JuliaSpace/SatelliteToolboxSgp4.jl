# ==========================================================================================
#                    ML-∂SGP4 Training Example
# ==========================================================================================
#
# Demonstrates the SatelliteToolboxSgp4 Lux extension for neural-network-corrected SGP4,
# based on the ML-dSGP4 paradigm from:
#
#   Acciarini, G., Baydin, A. G., & Izzo, D. (2025). "Closing the gap between SGP4 and
#   high-precision propagation via differentiable programming." Acta Astronautica, 226,
#   694-701.  https://doi.org/10.1016/j.actaastro.2024.10.063
#
# ==========================================================================================

using SatelliteToolboxSgp4
using Lux, Zygote, Optimisers
using StaticArrays
using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------------
# Data Generation
# ------------------------------------------------------------------------------------------

"""
    generate_training_data(tle::TLE, time_range_min) -> (vjd, vr_teme, vv_teme)

Generate reference ephemeris data by propagating a TLE through SGP4.

In a real application, replace this with numerically integrated ephemerides (e.g. from a
high-fidelity propagator) or observed tracking data. Here we use SGP4 itself as "truth"
for demonstration purposes.
"""
function generate_training_data(tle::TLE, time_range_min)
    sgp4d = sgp4_init(tle)
    epoch = sgp4d.epoch

    tsinces = collect(Float64, time_range_min)
    vjd     = epoch .+ tsinces ./ 1440.0

    vr_teme = Vector{SVector{3, Float64}}(undef, length(tsinces))
    vv_teme = Vector{SVector{3, Float64}}(undef, length(tsinces))

    for (k, t) in enumerate(tsinces)
        r, v = sgp4!(sgp4d, t)
        vr_teme[k] = r
        vv_teme[k] = v
    end

    return vjd, vr_teme, vv_teme
end

# ------------------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------------------

tle_amazonia = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
    """

println("=" ^ 80)
println("ML-∂SGP4 Training Example (Lux backend)")
println("=" ^ 80)
println()
println("NOTE: This example uses SGP4 output as 'truth' for demonstration.")
println("In practice, replace generate_training_data with high-precision")
println("numerically-integrated ephemerides or observed tracking data.")
println()

# -- Step 1: Generate training data --------------------------------------------------------

println("Generating training data...")
vjd, vr_teme, vv_teme = generate_training_data(tle_amazonia, 0.0:2.0:1440.0)
println("  $(length(vjd)) observations over 1 day\n")

# -- Step 2: Train (returns MLdSGP4Model, analogous to Sgp4Constants) ----------------------

config = MLdSGP4Config(hidden_layers = [64, 64])
println("Config:  $(config)\n")

println("Training...")
model = ml_dsgp4_train(
    tle_amazonia, vjd, vr_teme, vv_teme;
    config        = config,
    epochs        = 50,
    batch_size    = 32,
    learning_rate = 1e-3,
    verbose       = true,
)

# -- Step 3: Create propagator and evaluate (mirrors sgp4_init / sgp4!) --------------------

mlsgp4d = ml_dsgp4_init(tle_amazonia; model = model)

println("\nEvaluation at t = 0, 60, 120 min:")
sgp4d = sgp4_init(tle_amazonia)
for t in [0.0, 60.0, 120.0]
    r_ml, v_ml = ml_dsgp4!(mlsgp4d, t)
    r_sg, v_sg = sgp4!(sgp4d, t)
    Δr = norm(collect(r_ml) .- collect(r_sg))
    Δv = norm(collect(v_ml) .- collect(v_sg))
    @printf("  t = %6.1f min  |Δr| = %10.4f km  |Δv| = %12.6e km/s\n", t, Δr, Δv)
end

# -- Step 4: Save and reload ---------------------------------------------------------------

model_path = joinpath(@__DIR__, "ml_dsgp4_model.bin")
ml_dsgp4_save(model_path, model)
println("\nModel saved to: $(model_path)")

loaded  = ml_dsgp4_load(model_path)
mlsgp4d_loaded = ml_dsgp4_init(tle_amazonia; model = loaded)
r1, v1  = ml_dsgp4!(mlsgp4d, 60.0)
r2, v2  = ml_dsgp4!(mlsgp4d_loaded, 60.0)
@printf("Round-trip check (t=60 min): |Δr| = %.2e km\n", norm(collect(r1) .- collect(r2)))

# -- Step 5: Convenience API (mirrors sgp4(Δt, tle; sgp4c=...)) ----------------------------

r, v, mlsgp4d_new = ml_dsgp4(60.0, tle_amazonia; model = model)
@printf("\nml_dsgp4(60.0, tle; model=...) -> r = [%.3f, %.3f, %.3f] km\n", r[1], r[2], r[3])

println("\n" * "=" ^ 80)
