## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

using ForwardDiff
using StaticArrays

@testset "Aqua.jl" begin
    Aqua.test_all(
        SatelliteToolboxSgp4;
        ambiguities  = (recursive = false),
        deps_compat  = (check_extras = false),
    )
end

if VERSION >= v"1.12"
    @warn "JET.jl test skipped on Julia 1.12+ due to MethodTableView incompatibility"
else
    @testset "JET Testing" begin
        rep = JET.test_package(SatelliteToolboxSgp4; toplevel_logger=nothing, target_modules=(@__MODULE__,))
    end
end


if VERSION >= v"1.12"
    @warn "Allocation Check skipped on Julia 1.12+ as it is falsely flagging rem2pi internals"
else
    @testset "Allocation Check" begin
        # -- sgp4_init! (from orbital elements) ------------------------------------------------

        @test length(
            check_allocs(
                (sgp4d, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar) -> begin
                    sgp4_init!(sgp4d, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)
                end,
                (Sgp4Propagator{Float64, Float64}, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64)
            )
        ) == 0

        # -- sgp4! (propagation step) ----------------------------------------------------------

        @test length(
            check_allocs(
                (sgp4d, t) -> begin
                    sgp4!(sgp4d, t)
                end,
                (Sgp4Propagator{Float64, Float64}, Float64)
            )
        ) == 0

        # -- TLE Fitting: _init_sgp4_with_state_vector! ----------------------------------------

        @test length(
            check_allocs(
                (sgp4d, sv, epoch) -> begin
                    SatelliteToolboxSgp4._init_sgp4_with_state_vector!(sgp4d, sv, epoch)
                end,
                (Sgp4Propagator{Float64, Float64}, SVector{7, Float64}, Float64)
            )
        ) == 0

        # -- TLE Fitting: _sgp4_jacobian (FiniteDiffJacobian) ----------------------------------

        @test length(
            check_allocs(
                (sgp4d, Δt, x₁, y₁) -> begin
                    SatelliteToolboxSgp4._sgp4_jacobian(
                        FiniteDiffJacobian(), sgp4d, Δt, x₁, y₁
                    )
                end,
                (Sgp4Propagator{Float64, Float64}, Float64, SVector{7, Float64}, SVector{6, Float64})
            )
        ) == 0

        # -- TLE Fitting: _sgp4_fwd_jacobian_eval ----------------------------------------------

        _D = ForwardDiff.Dual{ForwardDiff.Tag{Nothing, Float64}, Float64, 7}

        @test length(
            check_allocs(
                (sgp4d_ad, epoch, Δt, x₁) -> begin
                    SatelliteToolboxSgp4._sgp4_fwd_jacobian_eval(sgp4d_ad, epoch, Δt, x₁)
                end,
                (Sgp4Propagator{Float64, _D}, Float64, Float64, SVector{7, Float64})
            )
        ) == 0

        # -- TLE Fitting: _sgp4_jacobian (ForwardDiffJacobian) with pre-allocated propagator -----

        @test length(
            check_allocs(
                (sgp4d, sgp4d_ad, Δt, x₁, y₁) -> begin
                    SatelliteToolboxSgp4._sgp4_jacobian(
                        ForwardDiffJacobian(), sgp4d, Δt, x₁, y₁;
                        sgp4d_ad = sgp4d_ad
                    )
                end,
                (Sgp4Propagator{Float64, Float64}, Sgp4Propagator{Float64, _D}, Float64, SVector{7, Float64}, SVector{6, Float64})
            )
        ) == 0

        # -- TLE Fitting: fit_sgp4_tle! (FiniteDiffJacobian) -----------------------------------
        # fit_sgp4_tle! inherently allocates. We use a regression bound here.

        @test length(
            check_allocs(
                (sgp4d, vjd, vr_teme, vv_teme) -> begin
                    fit_sgp4_tle!(
                        sgp4d, vjd, vr_teme, vv_teme;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    Sgp4Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                )
            )
        ) <= 45

        # -- TLE Fitting: fit_sgp4_tle! (ForwardDiffJacobian) ----------------------------------

        @test length(
            check_allocs(
                (sgp4d, vjd, vr_teme, vv_teme) -> begin
                    fit_sgp4_tle!(
                        sgp4d, vjd, vr_teme, vv_teme;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    Sgp4Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                )
            )
        ) <= 47

        # -- ML-dSGP4: ml_dsgp4! (zero model) -------------------------------------------------

        _tle_perf = tle"""
            AMAZONIA 1
            1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
            2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
            """

        _mlsgp4d_zero = ml_dsgp4_init(_tle_perf)
        _T_ml = typeof(_mlsgp4d_zero)

        @test length(
            check_allocs(
                (prop, t) -> begin
                    ml_dsgp4!(prop, t)
                end,
                (_T_ml, Float64)
            )
        ) == 0

        # -- ML-dSGP4: ml_dsgp4! (trained model) -----------------------------------------------

        _sgp4d_ref = sgp4_init(_tle_perf)
        _epoch_ref = _sgp4d_ref.epoch
        _times_ref = [0.0, 60.0, 120.0]
        _vjd_ref   = _epoch_ref .+ _times_ref ./ 1440.0
        _vr_ref    = Vector{SVector{3, Float64}}(undef, length(_times_ref))
        _vv_ref    = Vector{SVector{3, Float64}}(undef, length(_times_ref))
        for (k, t) in enumerate(_times_ref)
            r, v = sgp4!(_sgp4d_ref, t)
            _vr_ref[k] = r
            _vv_ref[k] = v
        end

        _trained_model = ml_dsgp4_train(
            _tle_perf, _vjd_ref, _vr_ref, _vv_ref;
            epochs  = 2,
            verbose = false,
        )
        _mlsgp4d_trained = ml_dsgp4_init(_tle_perf; model = _trained_model)
        _T_ml_trained = typeof(_mlsgp4d_trained)

        @test length(
            check_allocs(
                (prop, t) -> begin
                    ml_dsgp4!(prop, t)
                end,
                (_T_ml_trained, Float64)
            )
        ) == 0

    end
end
