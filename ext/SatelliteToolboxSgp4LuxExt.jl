module SatelliteToolboxSgp4LuxExt

using SatelliteToolboxSgp4
using Lux
using Optimisers
using Zygote

import Printf: @printf, @sprintf
import Random
import Random: seed!, randperm
import Serialization: serialize, deserialize

using StaticArrays

import SatelliteToolboxSgp4: ml_dsgp4_init, ml_dsgp4, ml_dsgp4!, ml_dsgp4_train,
                              ml_dsgp4_save, ml_dsgp4_load

# ==========================================================================================
#                                         Types
# ==========================================================================================
#
# The design mirrors the SGP4 type hierarchy exactly:
#
#   Sgp4Constants   ←→  MLdSGP4Model       (the "model": immutable config / weights)
#   Sgp4Propagator  ←→  MLdSGP4            (the "propagator": model + TLE + workspace)
#
# Neural-network-corrected SGP4 following the ML-∂SGP4 paradigm:
#
#   Acciarini, G., Baydin, A. G., & Izzo, D. (2025). "Closing the gap between SGP4 and
#   high-precision propagation via differentiable programming." Acta Astronautica, 226,
#   694-701.
#
#   Input correction:   x_corrected = x ⊙ (1 + α ⊙ tanh(NN_in(x)))
#   Output correction:  y_corrected = y ⊙ (1 + β ⊙ tanh(NN_out(y)))
#
# ==========================================================================================

"""
    MLdSGP4Model{C<:MLdSGP4Config, M, Ps, St, Tα<:Number, Tβ<:Number}

Trained correction weights for ML-∂SGP4.  Analogous to `Sgp4Constants`: this holds the
learned parameters that define **how** corrections are applied, independent of any specific
satellite.

Created by [`ml_dsgp4_train`](@ref), saved/loaded by [`ml_dsgp4_save`](@ref) /
[`ml_dsgp4_load`](@ref), and passed as a keyword to [`ml_dsgp4_init`](@ref) and
[`ml_dsgp4`](@ref).
"""
struct MLdSGP4Model{C<:MLdSGP4Config, M, Ps, St, Tα<:Number, Tβ<:Number}
    config::C
    input_net::M
    output_net::M
    input_ps::Ps
    output_ps::Ps
    input_st::St
    output_st::St
    α::Vector{Tα}
    β::Vector{Tβ}
end

"""
    MLdSGP4{C<:MLdSGP4Config, M, Ps, St, SPs, Tα<:Number, Tβ<:Number}

ML-corrected SGP4 propagator.  Analogous to `Sgp4Propagator`: this bundles a trained
[`MLdSGP4Model`](@ref) together with TLE-specific orbital data and an embedded
`Sgp4Propagator` workspace.

At construction time the Lux parameters are converted to `SMatrix`/`SVector` so that the
forward pass through the Lux `Chain` is entirely stack-allocated and allocation-free.

Created by [`ml_dsgp4_init`](@ref), propagated by [`ml_dsgp4!`](@ref).
"""
struct MLdSGP4{C<:MLdSGP4Config, TT<:Number, BT<:Number, ET<:Number, SgT<:Number, SgT2<:Number, M, Ps, St, SPs, Tα<:Number, Tβ<:Number}
    model::MLdSGP4Model{C, M, Ps, St, Tα, Tβ}
    tle_elements::SVector{6, TT}
    bstar::BT
    epoch::ET
    sgp4d::Sgp4Propagator{SgT, SgT2}
    input_ps::SPs
    output_ps::SPs
    α::SVector{6, Tα}
    β::SVector{6, Tβ}
end

# ==========================================================================================
#                                    Internal helpers
# ==========================================================================================

function _build_network(config::MLdSGP4Config)
    hl = config.hidden_layers
    layers = Any[Dense(6 => hl[1], leakyrelu)]
    for i in 2:length(hl)
        push!(layers, Dense(hl[i-1] => hl[i], leakyrelu))
    end
    push!(layers, Dense(hl[end] => 6))
    return Chain(layers...)
end

function _extract_tle_elements(tle::TLE)
    d2r = π / 180
    n₀ = Float64(tle.mean_motion) * (2π / 1440)
    e₀ = Float64(tle.eccentricity)
    i₀ = Float64(tle.inclination) * d2r
    Ω₀ = Float64(tle.raan) * d2r
    ω₀ = Float64(tle.argument_of_perigee) * d2r
    M₀ = Float64(tle.mean_anomaly) * d2r
    bstar = Float64(tle.bstar)
    epoch = Float64(tle_epoch(tle))
    return [e₀, ω₀, i₀, M₀, n₀, Ω₀], bstar, epoch
end

function _new_sgp4_workspace()
    sgp4d = Sgp4Propagator{Float64, Float64}()
    sgp4d.sgp4c  = sgp4c_wgs84
    sgp4d.sgp4ds = SatelliteToolboxSgp4.Sgp4DeepSpace{Float64}()
    return sgp4d
end

function _staticify_dense_params(lp)
    w = lp.weight
    b = lp.bias
    M, N = size(w)
    return (weight = SMatrix{M, N}(w), bias = SVector{M}(b))
end

_staticify_params(ps) = map(_staticify_dense_params, ps)

function _setup_model(config::MLdSGP4Config)
    rng = Random.default_rng()

    input_net  = _build_network(config)
    output_net = _build_network(config)

    input_ps, input_st   = Lux.setup(rng, input_net)
    output_ps, output_st = Lux.setup(rng, output_net)

    input_ps  = Lux.f64(input_ps)
    output_ps = Lux.f64(output_ps)
    input_st  = Lux.f64(input_st)
    output_st = Lux.f64(output_st)

    return input_net, output_net, input_ps, output_ps, input_st, output_st
end

# Random init with nonzero α/β — suitable starting point for training.
function _create_model(config::MLdSGP4Config)
    input_net, output_net, input_ps, output_ps, input_st, output_st =
        _setup_model(config)

    α = fill(Float64(config.input_correction), 6)
    β = fill(Float64(config.output_correction), 6)

    return MLdSGP4Model(
        config, input_net, output_net,
        input_ps, output_ps, input_st, output_st,
        α, β,
    )
end

# Zero α/β — corrections vanish so ML-∂SGP4 behaves identically to plain SGP4.
function _zero_model(config::MLdSGP4Config = MLdSGP4Config())
    input_net, output_net, input_ps, output_ps, input_st, output_st =
        _setup_model(config)

    return MLdSGP4Model(
        config, input_net, output_net,
        input_ps, output_ps, input_st, output_st,
        zeros(Float64, 6), zeros(Float64, 6),
    )
end

# ==========================================================================================
#                                      ml_dsgp4_init
# ==========================================================================================

"""
    ml_dsgp4_init(tle::TLE; model::MLdSGP4Model = _zero_model()) -> MLdSGP4

Create and initialize an ML-corrected SGP4 propagator from a TLE.

Mirrors `sgp4_init(tle; sgp4c=...)`: the TLE data is baked into the returned struct so that
subsequent calls to `ml_dsgp4!` only require the elapsed time.  The `model` keyword is the
ML equivalent of `sgp4c` — it provides the learned correction parameters.

When `model` is omitted a zero-correction model is used, so the propagator behaves
identically to plain SGP4.
"""
function ml_dsgp4_init(tle::TLE; model::MLdSGP4Model = _zero_model())
    tle_vec, bstar, epoch = _extract_tle_elements(tle)
    return MLdSGP4(
        model,
        SVector{6}(tle_vec),
        bstar,
        epoch,
        _new_sgp4_workspace(),
        _staticify_params(model.input_ps),
        _staticify_params(model.output_ps),
        SVector{6}(model.α),
        SVector{6}(model.β),
    )
end

# ==========================================================================================
#                                  Internal forward pass
# ==========================================================================================

function _ml_dsgp4_forward(mlsgp4d::MLdSGP4, Δt::Number)
    m   = mlsgp4d.model
    x   = mlsgp4d.tle_elements
    α   = mlsgp4d.α
    β   = mlsgp4d.β
    nR  = m.config.normalization_R
    nV  = m.config.normalization_V

    nn_in, _ = m.input_net(x, mlsgp4d.input_ps, m.input_st)

    x_corr = SVector{6}(
        x[1] * (1 + α[1] * tanh(nn_in[1])),
        x[2] * (1 + α[2] * tanh(nn_in[2])),
        x[3] * (1 + α[3] * tanh(nn_in[3])),
        x[4] * (1 + α[4] * tanh(nn_in[4])),
        x[5] * (1 + α[5] * tanh(nn_in[5])),
        x[6] * (1 + α[6] * tanh(nn_in[6])),
    )

    sgp4_init!(mlsgp4d.sgp4d, mlsgp4d.epoch, x_corr[5], x_corr[1], x_corr[3],
               x_corr[6], x_corr[2], x_corr[4], mlsgp4d.bstar)
    r_teme, v_teme = sgp4!(mlsgp4d.sgp4d, Δt)

    y_norm = SVector{6}(
        r_teme[1] / nR, r_teme[2] / nR, r_teme[3] / nR,
        v_teme[1] / nV, v_teme[2] / nV, v_teme[3] / nV,
    )

    nn_out, _ = m.output_net(y_norm, mlsgp4d.output_ps, m.output_st)

    r = SVector{3}(
        y_norm[1] * (1 + β[1] * tanh(nn_out[1])) * nR,
        y_norm[2] * (1 + β[2] * tanh(nn_out[2])) * nR,
        y_norm[3] * (1 + β[3] * tanh(nn_out[3])) * nR,
    )
    v = SVector{3}(
        y_norm[4] * (1 + β[4] * tanh(nn_out[4])) * nV,
        y_norm[5] * (1 + β[5] * tanh(nn_out[5])) * nV,
        y_norm[6] * (1 + β[6] * tanh(nn_out[6])) * nV,
    )

    return r, v
end

# ==========================================================================================
#                                     ml_dsgp4! / ml_dsgp4
# ==========================================================================================

"""
    ml_dsgp4!(mlsgp4d::MLdSGP4, Δt::Number) -> (r_teme, v_teme)

Propagate the ML-corrected SGP4 propagator to time `Δt` [min] from the TLE epoch.

Mirrors `sgp4!(sgp4d, Δt)`.

# Returns

- `SVector{3, Float64}`: Position vector [km] in the TEME frame.
- `SVector{3, Float64}`: Velocity vector [km/s] in the TEME frame.
"""
function ml_dsgp4!(mlsgp4d::MLdSGP4, Δt::Number)
    return _ml_dsgp4_forward(mlsgp4d, Δt)
end

"""
    ml_dsgp4(Δt::Number, tle::TLE; model::MLdSGP4Model = _zero_model()) -> (r_teme, v_teme, mlsgp4d)

Initialize an ML-corrected SGP4 propagator and propagate to time `Δt` [min].

Mirrors `sgp4(Δt, tle; sgp4c=...)`.

# Returns

- `SVector{3, Float64}`: Position vector [km] in the TEME frame.
- `SVector{3, Float64}`: Velocity vector [km/s] in the TEME frame.
- `MLdSGP4`: The initialized propagator (can be reused with `ml_dsgp4!`).
"""
function ml_dsgp4(Δt::Number, tle::TLE; model::MLdSGP4Model = _zero_model())
    mlsgp4d = ml_dsgp4_init(tle; model = model)
    r_teme, v_teme = ml_dsgp4!(mlsgp4d, Δt)
    return r_teme, v_teme, mlsgp4d
end

# ==========================================================================================
#                                      ml_dsgp4_train
# ==========================================================================================

"""
    ml_dsgp4_train(tle, vjd, vr_teme, vv_teme; kwargs...) -> MLdSGP4Model

Train ML correction weights against reference ephemeris data for the given TLE.

Returns a new [`MLdSGP4Model`](@ref) with updated parameters.

# Arguments

- `tle::TLE`: The satellite TLE to correct.
- `vjd`: Julian Day timestamps of the reference observations.
- `vr_teme`: Position vectors [km] in TEME at each timestamp.
- `vv_teme`: Velocity vectors [km/s] in TEME at each timestamp.

# Keywords

- `model::Union{MLdSGP4Model, Nothing} = nothing`: Starting model weights.  If `nothing`,
  a fresh model is created from `config`.
- `config::MLdSGP4Config = MLdSGP4Config()`: Architecture config (used only when
  `model === nothing`).
- `epochs::Int = 50`: Number of training epochs.
- `batch_size::Int = 32`: Mini-batch size.
- `learning_rate = 1e-3`: Adam learning rate.
- `verbose::Bool = true`: Print training progress.
- `seed::Int = 42`: RNG seed for reproducibility.
"""
function ml_dsgp4_train(
    tle::TLE,
    vjd::AbstractVector,
    vr_teme::AbstractVector,
    vv_teme::AbstractVector;
    model::Union{MLdSGP4Model, Nothing} = nothing,
    config::MLdSGP4Config               = MLdSGP4Config(),
    epochs::Int                         = 50,
    batch_size::Int                     = 32,
    learning_rate                       = 1e-3,
    verbose::Bool                       = true,
    seed::Int                           = 42,
)
    seed!(seed)

    mdl = isnothing(model) ? _create_model(config) : model
    cfg = mdl.config

    tle_elements, bstar_val, epoch_jd = _extract_tle_elements(tle)
    norm_R = cfg.normalization_R
    norm_V = cfg.normalization_V

    N = length(vjd)
    tsinces = [(vjd[k] - epoch_jd) * 1440.0 for k in 1:N]
    targets = [vcat(
        collect(vr_teme[k]) ./ norm_R,
        collect(vv_teme[k]) ./ norm_V,
    ) for k in 1:N]

    sgp4d = _new_sgp4_workspace()

    ps = (
        input  = mdl.input_ps,
        output = mdl.output_ps,
        α      = copy(mdl.α),
        β      = copy(mdl.β),
    )
    st = (
        input  = mdl.input_st,
        output = mdl.output_st,
    )

    input_net  = mdl.input_net
    output_net = mdl.output_net

    opt_state = Optimisers.setup(Optimisers.Adam(learning_rate), ps)

    if verbose
        @printf("  %5s  %14s  %14s  %14s\n",
            "Epoch", "Train Loss", "Pos RMSE [km]", "Vel RMSE [km/s]")
        @printf("  %5s  %14s  %14s  %14s\n",
            "-----", "--------------", "--------------", "--------------")
    end

    best_loss = Inf

    for ep in 1:epochs
        indices    = randperm(N)
        epoch_loss = 0.0
        n_batches  = 0

        for batch_start in 1:batch_size:N
            batch_end = min(batch_start + batch_size - 1, N)
            batch_idx = indices[batch_start:batch_end]

            loss, grads = Zygote.withgradient(ps) do p
                bl = 0.0
                for k in batch_idx
                    nn_in, _ = input_net(tle_elements, p.input, st.input)
                    δ_in  = p.α .* tanh.(nn_in)
                    x_corr = tle_elements .* (1 .+ δ_in)

                    sgp4_init!(sgp4d, epoch_jd, x_corr[5], x_corr[1], x_corr[3],
                               x_corr[6], x_corr[2], x_corr[4], bstar_val)
                    r_teme, v_teme = sgp4!(sgp4d, tsinces[k])

                    y_norm = vcat(
                        collect(r_teme) ./ norm_R,
                        collect(v_teme) ./ norm_V,
                    )

                    nn_out, _ = output_net(y_norm, p.output, st.output)
                    δ_out  = p.β .* tanh.(nn_out)
                    y_corr = y_norm .* (1 .+ δ_out)

                    bl += sum((y_corr .- targets[k]).^2)
                end
                bl / length(batch_idx)
            end

            opt_state, ps = Optimisers.update(opt_state, ps, grads[1])
            epoch_loss += loss
            n_batches  += 1
        end

        avg_loss = epoch_loss / n_batches
        avg_loss < best_loss && (best_loss = avg_loss)

        if verbose
            eval_mdl = MLdSGP4Model(
                cfg, input_net, output_net,
                ps.input, ps.output, st.input, st.output,
                ps.α, ps.β,
            )
            eval_prop = MLdSGP4(
                eval_mdl,
                SVector{6}(tle_elements),
                bstar_val,
                epoch_jd,
                sgp4d,
                _staticify_params(ps.input),
                _staticify_params(ps.output),
                SVector{6}(ps.α),
                SVector{6}(ps.β),
            )
            pos_sse = 0.0
            vel_sse = 0.0
            for k in 1:N
                r, v = _ml_dsgp4_forward(eval_prop, tsinces[k])
                dr = collect(r) .- collect(vr_teme[k])
                dv = collect(v) .- collect(vv_teme[k])
                pos_sse += sum(dr .^ 2)
                vel_sse += sum(dv .^ 2)
            end
            @printf("  %5d  %14.6e  %14.6e  %14.6e\n",
                ep, avg_loss, sqrt(pos_sse / N), sqrt(vel_sse / N))
        end
    end

    verbose && @printf("\n  Training complete. Best loss: %.6e\n", best_loss)

    return MLdSGP4Model(
        cfg, input_net, output_net,
        ps.input, ps.output, st.input, st.output,
        ps.α, ps.β,
    )
end

# ==========================================================================================
#                                   ml_dsgp4_save / ml_dsgp4_load
# ==========================================================================================

"""
    ml_dsgp4_save(path::AbstractString, model::MLdSGP4Model)

Serialize an [`MLdSGP4Model`](@ref) to `path`.
"""
function ml_dsgp4_save(path::AbstractString, model::MLdSGP4Model)
    state = (;
        config    = model.config,
        input_ps  = model.input_ps,
        output_ps = model.output_ps,
        input_st  = model.input_st,
        output_st = model.output_st,
        α         = model.α,
        β         = model.β,
    )
    open(path, "w") do io
        serialize(io, state)
    end
    return nothing
end

"""
    ml_dsgp4_load(path::AbstractString) -> MLdSGP4Model

Deserialize an [`MLdSGP4Model`](@ref) from `path`.

Bind the returned model to a TLE with `ml_dsgp4_init(tle; model=...)` before propagating.
"""
function ml_dsgp4_load(path::AbstractString)
    state = open(deserialize, path, "r")

    cfg = state.config
    input_net  = _build_network(cfg)
    output_net = _build_network(cfg)

    return MLdSGP4Model(
        cfg, input_net, output_net,
        state.input_ps, state.output_ps,
        state.input_st, state.output_st,
        copy(state.α), copy(state.β),
    )
end

end # module
