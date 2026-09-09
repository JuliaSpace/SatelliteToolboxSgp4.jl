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

@testset "Fitting Using OMM" verbose = true begin
    omm_input = read_omm(_OMM_FIXTURE)
    tle_input = convert(TLE, omm_input)

    # Generate the osculating elements (TEME).
    sgp4d   = sgp4_init(omm_input)
    ret     = map(t -> sgp4!(sgp4d, t), 0:0.2:200)
    vr_teme = first.(ret)
    vv_teme = last.(ret)
    vjd     = sgp4d.epoch .+ (0:0.2:200) ./ 1440

    @testset "Without Template" begin
        omm, P, stats = fit_sgp4_mean_elements(
            OrbitMeanElementsMessage,
            vjd,
            vr_teme,
            vv_teme;
            atol                = 1e-10,
            rtol                = 1e-10,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 1000,
            verbose             = false,
        )

        @test omm isa OrbitMeanElementsMessage
        @test P isa SMatrix{7, 7, Float64}

        # The fit must converge within the allowed iterations.
        @test stats.converged
        @test 1 < stats.iterations < 1000
        @test stats.position_rmse isa Float64
        @test stats.velocity_rmse isa Float64
        @test stats.total_rmse isa Float64
        @test 0 < stats.position_rmse < 1e-3
        @test 0 < stats.velocity_rmse < 1e-6
        @test stats.total_rmse ≈ √(stats.position_rmse^2 + stats.velocity_rmse^2)

        # The default metadata mirrors the default TLE fields.
        @test omm.originator == "SatelliteToolboxSgp4.jl"
        @test omm.object_name == "UNDEFINED"
        @test omm.center_name == "EARTH"
        @test omm.ref_frame == "TEME"
        @test omm.time_system == "UTC"
        @test omm.mean_element_theory == "SGP4"
        @test omm.classification_type == 'U'
        @test omm.norad_cat_id == 9999
        @test omm.mean_motion_dot == 0
        @test omm.mean_motion_ddot == 0

        # The fitted elements must match the input message.
        @test omm.bstar ≈ omm_input.bstar atol = 1e-6
        @test omm.eccentricity ≈ omm_input.eccentricity atol = 1e-7
        @test omm.inclination ≈ omm_input.inclination atol = 1e-4
        @test omm.raan ≈ omm_input.raan atol = 1e-4
        @test omm.arg_of_pericenter ≈ omm_input.arg_of_pericenter atol = 1e-4
        @test omm.mean_anomaly ≈ omm_input.mean_anomaly atol = 1e-4
        @test omm.mean_motion ≈ omm_input.mean_motion atol = 1e-7

        # The epoch is stored with millisecond precision.
        @test abs(Dates.value(DateTime(omm.epoch) - DateTime(omm_input.epoch))) <= 1

        # The covariance section must contain the position and velocity block of `P`.
        cov = omm.covariance_matrix

        @test !isnothing(cov)
        @test cov.cov_ref_frame == "TEME"
        @test cov.cx_x == P[1, 1]
        @test cov.cy_x == P[2, 1]
        @test cov.cz_z == P[3, 3]
        @test cov.cx_dot_x == P[4, 1]
        @test cov.cy_dot_y_dot == P[5, 5]
        @test cov.cz_dot_z_dot == P[6, 6]

        # The message must be usable to initialize the propagator.
        sgp4d_fit = sgp4_init(omm)

        @test sgp4d_fit.bstar ≈ sgp4d.bstar atol = 1e-6
    end

    @testset "With Template and Without Covariance" begin
        omm, ~ = fit_sgp4_mean_elements(
            OrbitMeanElementsMessage,
            vjd,
            vr_teme,
            vv_teme;
            atol                = 1e-10,
            rtol                = 1e-10,
            include_covariance  = false,
            initial_guess       = omm_input,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 10,
            template            = omm_input,
            verbose             = false,
        )

        # The metadata must be copied from the template.
        @test omm.originator == omm_input.originator
        @test omm.object_name == omm_input.object_name
        @test omm.object_id == omm_input.object_id
        @test omm.norad_cat_id == omm_input.norad_cat_id
        @test omm.element_set_number == omm_input.element_set_number
        @test omm.rev_at_epoch == omm_input.rev_at_epoch
        @test omm.mean_motion_dot == 0

        @test isnothing(omm.covariance_matrix)

        @test omm.eccentricity ≈ omm_input.eccentricity atol = 1e-7
        @test omm.mean_motion ≈ omm_input.mean_motion atol = 1e-7
        @test omm.mean_anomaly ≈ omm_input.mean_anomaly atol = 1e-4

        # Both representations must provide the same mean elements.
        tle, ~ = fit_sgp4_mean_elements(
            TLE,
            vjd,
            vr_teme,
            vv_teme;
            atol                = 1e-10,
            rtol                = 1e-10,
            initial_guess       = omm_input,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 10,
            template            = tle_input,
            verbose             = false,
        )

        @test omm.mean_motion ≈ tle.mean_motion atol = 1e-10
        @test omm.eccentricity ≈ tle.eccentricity atol = 1e-10
        @test omm.inclination ≈ tle.inclination atol = 1e-10
        @test omm.raan ≈ tle.raan atol = 1e-10
        @test omm.arg_of_pericenter ≈ tle.argument_of_perigee atol = 1e-10
        @test omm.mean_anomaly ≈ tle.mean_anomaly atol = 1e-10
        @test omm.bstar ≈ tle.bstar atol = 1e-12
    end

    @testset "Default Sink" begin
        kwargs = (;
            atol                = 1e-10,
            rtol                = 1e-10,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 1000,
            verbose             = false,
        )

        # Without a sink type, the mean elements must be returned as an OMM.
        omm_ref, P_ref = fit_sgp4_mean_elements(
            OrbitMeanElementsMessage, vjd, vr_teme, vv_teme; kwargs...
        )
        omm, P = fit_sgp4_mean_elements(vjd, vr_teme, vv_teme; kwargs...)

        @test omm isa OrbitMeanElementsMessage
        @test P == P_ref
        @test omm.mean_motion == omm_ref.mean_motion
        @test omm.mean_anomaly == omm_ref.mean_anomaly
        @test omm.covariance_matrix == omm_ref.covariance_matrix

        sgp4d  = Sgp4Propagator{Float64}(SGP4C_WGS84)
        omm, P = fit_sgp4_mean_elements!(sgp4d, vjd, vr_teme, vv_teme; kwargs...)

        @test omm isa OrbitMeanElementsMessage
        @test P == P_ref
        @test sgp4d.epoch ≈ vjd[begin] atol = 1e-9
    end

    @testset "Statistics" begin
        # A single iteration cannot converge, since the convergence is checked against the
        # residue of the previous iteration.
        omm, P, stats = fit_sgp4_mean_elements(
            vjd,
            vr_teme,
            vv_teme;
            mean_elements_epoch = vjd[begin],
            max_iterations      = 1,
            verbose             = false,
        )

        @test !stats.converged
        @test stats.iterations == 1
        @test stats.total_rmse > 0

        # Stopping by the iteration limit must also be reported when it is reached later.
        omm, P, stats = fit_sgp4_mean_elements(
            vjd,
            vr_teme,
            vv_teme;
            atol                = 0,
            rtol                = 0,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 3,
            verbose             = false,
        )

        @test !stats.converged
        @test stats.iterations == 3

        # The number type of the statistics must follow the constants. The drag term is not
        # estimated and the Jacobian is obtained by automatic differentiation, since the
        # finite differences are not accurate enough in `Float32` to keep the normal
        # equations well conditioned.
        omm, P, stats = fit_sgp4_mean_elements(
            vjd,
            vr_teme,
            vv_teme;
            sgp4c               = Sgp4Constants{Float32}(SGP4C_WGS84),
            estimate_bstar      = false,
            jacobian_method     = ForwardDiffJacobian(),
            mean_elements_epoch = vjd[begin],
            max_iterations      = 2,
            verbose             = false,
        )

        @test P isa SMatrix{7, 7, Float32}
        @test stats.total_rmse isa Float32
    end

    @testset "Constants" begin
        kwargs = (;
            atol                = 1e-10,
            rtol                = 1e-10,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 1000,
            verbose             = false,
        )

        # The non-allocating version must use the constants of the propagator.
        sgp4d = Sgp4Propagator(SGP4C_WGS72)
        omm_ref, P_ref, stats_ref = fit_sgp4_mean_elements!(
            sgp4d, vjd, vr_teme, vv_teme; kwargs...
        )

        @test sgp4d.sgp4c === SGP4C_WGS72

        # The allocating version must accept the constants through the keyword.
        omm, P, stats = fit_sgp4_mean_elements(
            vjd, vr_teme, vv_teme; sgp4c = SGP4C_WGS72, kwargs...
        )

        @test omm.mean_motion == omm_ref.mean_motion
        @test omm.mean_anomaly == omm_ref.mean_anomaly
        @test P == P_ref
        @test stats == stats_ref

        # Different constants must lead to different mean elements.
        omm_wgs84, ~ = fit_sgp4_mean_elements(vjd, vr_teme, vv_teme; kwargs...)

        @test omm_wgs84.mean_motion != omm.mean_motion
    end

    @testset "DateTime Epoch" begin
        kwargs = (; atol = 1e-10, rtol = 1e-10, max_iterations = 1000, verbose = false)

        epoch_jd = vjd[begin]
        epoch_dt = julian2datetime(epoch_jd)

        omm_ref, ~ = fit_sgp4_mean_elements(
            vjd, vr_teme, vv_teme; mean_elements_epoch = epoch_jd, kwargs...
        )

        omm, ~ = fit_sgp4_mean_elements(
            vjd, vr_teme, vv_teme; mean_elements_epoch = epoch_dt, kwargs...
        )

        # The `DateTime` conversion truncates the epoch to milliseconds.
        @test DateTime(omm.epoch) == DateTime(omm_ref.epoch)
        @test omm.mean_motion ≈ omm_ref.mean_motion atol = 1e-10
        @test omm.mean_anomaly ≈ omm_ref.mean_anomaly atol = 1e-6
    end

    @testset "NamedTuple Template" begin
        kwargs = (;
            atol                = 1e-10,
            rtol                = 1e-10,
            mean_elements_epoch = vjd[begin],
            max_iterations      = 1000,
            verbose             = false,
        )

        omm, ~ = fit_sgp4_mean_elements(
            vjd,
            vr_teme,
            vv_teme;
            template = (;
                object_name     = "AMAZONIA 1",
                object_id       = "2021-015A",
                norad_cat_id    = 47699,
                originator      = "INPE",
                header_comments = ["Fitted from osculating data."],
            ),
            kwargs...,
        )

        # The provided fields must override the defaults, whereas the others are kept.
        @test omm.object_name == "AMAZONIA 1"
        @test omm.object_id == "2021-015A"
        @test omm.norad_cat_id == 47699
        @test omm.originator == "INPE"
        @test omm.header.comments == ["Fitted from osculating data."]
        @test omm.center_name == "EARTH"
        @test omm.mean_element_theory == "SGP4"
        @test omm.mean_motion ≈ omm_input.mean_motion atol = 1e-7

        # The same applies to the TLE.
        tle, ~ = fit_sgp4_mean_elements(
            TLE,
            vjd,
            vr_teme,
            vv_teme;
            template = (; name = "AMAZONIA 1", satellite_number = 47699),
            kwargs...,
        )

        @test tle.name == "AMAZONIA 1"
        @test tle.satellite_number == 47699
        @test tle.classification == 'U'
        @test tle.international_designator == "999999"
        @test tle.mean_motion ≈ omm_input.mean_motion atol = 1e-7

        # Fields set by the fit cannot be overridden.
        @test_throws ArgumentError fit_sgp4_mean_elements(
            vjd, vr_teme, vv_teme; template = (; eccentricity = 0.1), kwargs...
        )
        @test_throws ArgumentError fit_sgp4_mean_elements(
            vjd, vr_teme, vv_teme; template = (; mean_element_theory = "DSST"), kwargs...
        )
        @test_throws ArgumentError fit_sgp4_mean_elements(
            TLE, vjd, vr_teme, vv_teme; template = (; bstar = 0.1), kwargs...
        )
    end

    @testset "Epoch Update" begin
        new_epoch = DateTime(omm_input.epoch) + Day(1)

        omm = update_sgp4_mean_elements_epoch(omm_input, new_epoch; verbose = false)
        tle = update_sgp4_mean_elements_epoch(tle_input, new_epoch; verbose = false)

        @test omm.object_name == omm_input.object_name
        @test omm.norad_cat_id == omm_input.norad_cat_id
        @test DateTime(omm.epoch) == new_epoch
        @test isnothing(omm.covariance_matrix)

        # The update must not depend on the representation. The angles differ slightly
        # because the TLE epoch is rounded to the resolution of its day fraction, whereas
        # the OMM epoch keeps the microseconds.
        @test omm.mean_motion ≈ tle.mean_motion atol = 1e-10
        @test omm.eccentricity ≈ tle.eccentricity atol = 1e-10
        @test omm.inclination ≈ tle.inclination atol = 1e-10
        @test omm.raan ≈ tle.raan atol = 1e-5
        @test omm.arg_of_pericenter ≈ tle.argument_of_perigee atol = 1e-5
        @test omm.mean_anomaly ≈ tle.mean_anomaly atol = 1e-5
        @test omm.bstar == omm_input.bstar

        # The in-place version must initialize the propagator with the updated message.
        sgp4d = Sgp4Propagator{Float64}(SGP4C_WGS84)
        omm!  = update_sgp4_mean_elements_epoch!(sgp4d, omm_input, new_epoch; verbose = false)

        @test sgp4d.epoch ≈ datetime2julian(new_epoch) atol = 1e-9
        @test omm!.mean_motion == omm.mean_motion

        # The Jacobian can also be computed by automatic differentiation.
        omm_ad = update_sgp4_mean_elements_epoch(
            omm_input, new_epoch; jacobian_method = ForwardDiffJacobian(), verbose = false
        )

        @test omm_ad.mean_motion ≈ omm.mean_motion atol = 1e-10
        @test omm_ad.mean_anomaly ≈ omm.mean_anomaly atol = 1e-6

        # The finite-difference settings must be forwarded.
        omm_fd = update_sgp4_mean_elements_epoch(
            omm_input,
            new_epoch;
            jacobian_perturbation     = 1e-4,
            jacobian_perturbation_tol = 1e-8,
            verbose                   = false,
        )

        @test omm_fd.mean_motion ≈ omm.mean_motion atol = 1e-10
        @test omm_fd.mean_anomaly ≈ omm.mean_anomaly atol = 1e-6

        # The constants can be selected in the allocating version.
        omm_wgs72 = update_sgp4_mean_elements_epoch(
            omm_input, new_epoch; sgp4c = SGP4C_WGS72, verbose = false
        )

        sgp4d = Sgp4Propagator(SGP4C_WGS72)
        omm_wgs72! = update_sgp4_mean_elements_epoch!(
            sgp4d, omm_input, new_epoch; verbose = false
        )

        @test omm_wgs72.mean_motion == omm_wgs72!.mean_motion
        @test omm_wgs72.mean_motion != omm.mean_motion
    end
end
