## Description #############################################################################
#
# Tests related to the SGP4 initialization using Orbit Mean-Elements Messages (OMMs).
#
############################################################################################

# The fixture is a public message distributed by Space-Track for the satellite AMAZONIA 1.
const _OMM_FIXTURE = "./omm_tests/amazonia_1.xml"

@testset "Initialization Using OMM" verbose = true begin
    omm = read_omm(_OMM_FIXTURE)
    tle = convert(TLE, omm)

    @testset "Mean Elements" begin
        sgp4d_omm = sgp4_init(omm)
        sgp4d_tle = sgp4_init(tle)

        @test sgp4d_omm isa Sgp4Propagator{Float64, Float64}

        # The epoch in the OMM has microsecond precision, whereas the TLE epoch is limited
        # by the floating-point representation of the day fraction.
        expected_epoch =
            datetime2julian(DateTime(2025, 12, 30, 18, 12, 4, 533)) + 984e-6 / 86400

        @test sgp4d_omm.epoch ≈ expected_epoch atol = 1e-9
        @test sgp4d_omm.epoch ≈ sgp4d_tle.epoch atol = 1e-8

        # The remaining elements are converted with the same expressions.
        @test sgp4d_omm.n₀ == sgp4d_tle.n₀
        @test sgp4d_omm.e₀ == sgp4d_tle.e₀
        @test sgp4d_omm.i₀ == sgp4d_tle.i₀
        @test sgp4d_omm.Ω₀ == sgp4d_tle.Ω₀
        @test sgp4d_omm.ω₀ == sgp4d_tle.ω₀
        @test sgp4d_omm.M₀ == sgp4d_tle.M₀
        @test sgp4d_omm.bstar == sgp4d_tle.bstar
        @test sgp4d_omm.bstar == 0.0001533
    end

    @testset "Propagation" begin
        sgp4d_omm = sgp4_init(omm)
        sgp4d_tle = sgp4_init(tle)

        # Both propagators must provide the same state. The epochs can differ by one unit
        # in the last place of the Julian Day (about 40 μs) due to the different rounding
        # paths, which is covered by the tolerances.
        for t in (0.0, 100.0, 1440.0)
            r_omm, v_omm = sgp4!(sgp4d_omm, t)
            r_tle, v_tle = sgp4!(sgp4d_tle, t)

            @test r_omm ≈ r_tle atol = 1e-3
            @test v_omm ≈ v_tle atol = 1e-6

            # The simultaneous initialization and propagation must match.
            r, v, sgp4d = sgp4(t, omm)

            @test r == r_omm
            @test v == v_omm
            @test sgp4d.epoch == sgp4d_omm.epoch
        end
    end

    @testset "In-place Initialization" begin
        sgp4d = sgp4_init(tle)
        sgp4d_omm = sgp4_init(omm)

        # Modify the propagator and initialize it again using the OMM.
        sgp4!(sgp4d, 100.0)
        sgp4_init!(sgp4d, omm)

        @test sgp4d.epoch == sgp4d_omm.epoch
        @test sgp4d.n₀ == sgp4d_omm.n₀
        @test sgp4d.e₀ == sgp4d_omm.e₀
        @test sgp4d.i₀ == sgp4d_omm.i₀
        @test sgp4d.Ω₀ == sgp4d_omm.Ω₀
        @test sgp4d.ω₀ == sgp4d_omm.ω₀
        @test sgp4d.M₀ == sgp4d_omm.M₀
        @test sgp4d.bstar == sgp4d_omm.bstar
        @test sgp4!(sgp4d, 10.0) == sgp4!(sgp4d_omm, 10.0)
    end

    @testset "Mean Motion From Semi-Major Axis" begin
        a  = 7134.084
        GM = 398600.4418

        omm_a = OrbitMeanElementsMessage(
            omm; mean_motion = nothing, semi_major_axis = a, GM = GM
        )

        sgp4d = sgp4_init(omm_a)

        @test sgp4d.n₀ ≈ √(GM / a^3) * 60 atol = 1e-12
    end

    @testset "Missing Drag Term" begin
        omm_no_bstar = OrbitMeanElementsMessage(
            omm;
            ephemeris_type      = nothing,
            classification_type = nothing,
            norad_cat_id        = nothing,
            element_set_number  = nothing,
            rev_at_epoch        = nothing,
            bstar               = nothing,
            mean_motion_dot     = nothing,
            mean_motion_ddot    = nothing,
        )

        sgp4d = sgp4_init(omm_no_bstar)

        @test sgp4d.bstar == 0
    end

    @testset "Other Constants" begin
        sgp4d = sgp4_init(omm; sgp4c = Sgp4Constants{Float32}(SGP4C_WGS84))

        @test sgp4d isa Sgp4Propagator{Float64, Float32}
        @test sgp4d.epoch isa Float64

        r, v, sgp4d = sgp4(10.0, omm; sgp4c = SGP4C_WGS72)

        @test sgp4d.sgp4c === SGP4C_WGS72
        @test eltype(r) === Float64
    end

    @testset "Errors" begin
        # == Mean Element Theory Different From SGP4 =======================================

        omm_theory = OrbitMeanElementsMessage(omm; mean_element_theory = "SPECIAL")

        @test_throws ArgumentError sgp4_init(omm_theory)
        @test_throws ArgumentError sgp4(10.0, omm_theory)

        # == Semi-Major Axis Without the Gravitational Coefficient =========================

        omm_no_gm = OrbitMeanElementsMessage(
            omm; mean_motion = nothing, semi_major_axis = 7134.084, GM = nothing
        )

        @test_throws ArgumentError sgp4_init(omm_no_gm)
    end
end
