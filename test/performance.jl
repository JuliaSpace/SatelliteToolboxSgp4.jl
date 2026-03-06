## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

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

    end
end
