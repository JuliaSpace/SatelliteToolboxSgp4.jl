## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

using ForwardDiff
using StaticArrays

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxSgp4; ambiguities=(recursive = false), deps_compat=(check_extras = false))
end

@testset "JET Testing" begin
    rep = JET.test_package(SatelliteToolboxSgp4; toplevel_logger=nothing, target_modules=(@__MODULE__,))
end

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
    # Called 7× per Jacobian evaluation in the inner loop of fit_sgp4_tle!.

    @test length(
        check_allocs(
            (sgp4d, sv, epoch) -> begin
                SatelliteToolboxSgp4._init_sgp4_with_state_vector!(sgp4d, sv, epoch)
            end,
            (Sgp4Propagator{Float64, Float64}, SVector{7, Float64}, Float64)
        )
    ) == 0

    # -- TLE Fitting: _sgp4_jacobian (FiniteDiffJacobian) ----------------------------------
    # The MMatrix{6, 7} scratch buffer is the sole expected allocation source.

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
    # With a pre-allocated Dual-typed propagator, the ForwardDiff Jacobian evaluation is
    # allocation-free — matching the FiniteDiff path (minus the MMatrix scratch buffer).

    _D = ForwardDiff.Dual{ForwardDiff.Tag{Nothing, Float64}, Float64, 7}

    @test length(
        check_allocs(
            (sgp4d_ad, epoch, Δt, x₁) -> begin
                SatelliteToolboxSgp4._sgp4_fwd_jacobian_eval(sgp4d_ad, epoch, Δt, x₁)
            end,
            (Sgp4Propagator{Float64, _D}, Float64, Float64, SVector{7, Float64})
        )
    ) == 0
end
