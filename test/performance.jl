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
    end
end

@testset "Fitting Allocations" begin
    # The fitting allocates only a fixed set of workspaces. Hence, the number of allocations
    # must not depend on the number of measurements or iterations. We check this at runtime
    # because the static analysis of AllocCheck also follows the code paths used only when
    # `verbose = true`, whose allocations inside Base I/O internals vary with the platform,
    # the loaded packages, and the AllocCheck version.

    tle_input = tle"""
        AMAZONIA 1
        1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
        2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
        """

    function fit_allocations(
        ::Type{S}, jacobian_method, num_measurements, num_iterations
    ) where {S}
        sgp4d   = sgp4_init(tle_input)
        vt      = range(0, 200; length = num_measurements)
        ret     = map(t -> sgp4!(sgp4d, t), vt)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = collect(sgp4d.epoch .+ vt ./ 1440)

        # Setting `atol` and `rtol` to 0 disables the convergence check, so the algorithm
        # performs exactly `num_iterations` iterations.
        fit!() = fit_sgp4_mean_elements!(
            sgp4d,
            S,
            vjd,
            vr_teme,
            vv_teme;
            atol            = 0,
            rtol            = 0,
            jacobian_method = jacobian_method,
            max_iterations  = num_iterations,
            verbose         = false,
        )

        # Compile the function before measuring the allocations.
        fit!()

        return @allocations fit!()
    end

    for S in (TLE, OrbitMeanElementsMessage),
        jacobian_method in (FiniteDiffJacobian(), ForwardDiffJacobian())

        num_allocations = fit_allocations(S, jacobian_method, 50, 2)

        @test fit_allocations(S, jacobian_method, 200, 2) == num_allocations
        @test fit_allocations(S, jacobian_method, 50, 6) == num_allocations
    end
end
