## Description #############################################################################
#
# Tests related to memory allocations.
#
############################################################################################

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
                (
                    Sgp4Propagator{Float64, Float64},
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                ),
            ),
        ) == 0

        # -- sgp4! (propagation step) ----------------------------------------------------------

        @test length(
            check_allocs(
                (sgp4d, t) -> begin
                    sgp4!(sgp4d, t)
                end,
                (Sgp4Propagator{Float64, Float64}, Float64),
            ),
        ) == 0

        # -- TLE Fitting: _init_sgp4_with_state_vector! ----------------------------------------

        @test length(
            check_allocs(
                (sgp4d, sv, epoch) -> begin
                    SatelliteToolboxSgp4._init_sgp4_with_state_vector!(sgp4d, sv, epoch)
                end,
                (Sgp4Propagator{Float64, Float64}, SVector{7, Float64}, Float64),
            ),
        ) == 0

        # -- Fitting: _sgp4_jacobian! (FiniteDiffJacobian) -------------------------------------

        @test length(
            check_allocs(
                (vJ, sgp4d, vjd, epoch, x₁, vŷ) -> begin
                    SatelliteToolboxSgp4._sgp4_jacobian!(
                        FiniteDiffJacobian(), vJ, sgp4d, nothing, vjd, epoch, x₁, vŷ
                    )
                end,
                (
                    Array{Float64, 3},
                    Sgp4Propagator{Float64, Float64},
                    Vector{Float64},
                    Float64,
                    SVector{7, Float64},
                    Vector{SVector{6, Float64}},
                ),
            ),
        ) == 0

        # -- Fitting: _sgp4_jacobian! (ForwardDiffJacobian) ------------------------------------

        _D = ForwardDiff.Dual{ForwardDiff.Tag{Nothing, Float64}, Float64, 7}

        @test length(
            check_allocs(
                (vJ, sgp4d, sgp4d_ad, vjd, epoch, x₁, vŷ) -> begin
                    SatelliteToolboxSgp4._sgp4_jacobian!(
                        ForwardDiffJacobian(), vJ, sgp4d, sgp4d_ad, vjd, epoch, x₁, vŷ
                    )
                end,
                (
                    Array{Float64, 3},
                    Sgp4Propagator{Float64, Float64},
                    Sgp4Propagator{Float64, _D},
                    Vector{Float64},
                    Float64,
                    SVector{7, Float64},
                    Vector{SVector{6, Float64}},
                ),
            ),
        ) == 0

        # -- Fitting: fit_sgp4_mean_elements! (FiniteDiffJacobian) -----------------------------
        # fit_sgp4_mean_elements! inherently allocates. We use a regression bound here.

        @test length(
            check_allocs(
                (sgp4d, vjd, vr_teme, vv_teme) -> begin
                    fit_sgp4_mean_elements!(
                        sgp4d,
                        TLE,
                        vjd,
                        vr_teme,
                        vv_teme;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    Sgp4Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 47

        # -- Fitting: fit_sgp4_mean_elements! (ForwardDiffJacobian) ----------------------------

        @test length(
            check_allocs(
                (sgp4d, vjd, vr_teme, vv_teme) -> begin
                    fit_sgp4_mean_elements!(
                        sgp4d,
                        TLE,
                        vjd,
                        vr_teme,
                        vv_teme;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    Sgp4Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 49
    end
end
