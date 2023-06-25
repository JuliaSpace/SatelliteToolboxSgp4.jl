# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Test related to the SGP4 TLE.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@testset "Fitting SGP4 TLEs" verbose = true begin
    # We just need to run the SGP4, obtain the osculating elements, convert them to TLE, and
    # compare with the original TLE.

    @testset "Without Initial Guess" begin
        # Scenario 1: Low Earth Orbit
        # ==================================================================================

        tle_input = tle"""
            AMAZONIA 1
            1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
            2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
            """

        # Generate the osculating elements (TEME).
        sgp4d   = sgp4_init(tle_input)
        ret     = map(t -> sgp4!(sgp4d, t), 0:0.2:200)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:0.2:200) ./ 1440

        # Obtain the mean elements.
        tle, ~ = fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            atol                     = 1e-10,
            rtol                     = 1e-10,
            element_set_number       = 999,
            international_designator = "21015A",
            mean_elements_epoch      = vjd[begin],
            name                     = "AMAZONIA 1",
            revolution_number        = 3043,
            satellite_number         = 47699,
            max_iterations           = 1000,
            verbose                  = false,
        )

        # Compare.
        @test tle.classification           == tle_input.classification
        @test tle.element_set_number       == tle_input.element_set_number
        @test tle.epoch_year               == tle_input.epoch_year
        @test tle.international_designator == tle_input.international_designator
        @test tle.name                     == tle_input.name
        @test tle.revolution_number        == tle_input.revolution_number
        @test tle.satellite_number         == tle_input.satellite_number

        @test tle.bstar               ≈  tle_input.bstar               atol = 1e-6
        @test tle.eccentricity        ≈  tle_input.eccentricity        atol = 1e-7
        @test tle.epoch_day           ≈  tle_input.epoch_day           atol = 1e-8
        @test tle.inclination         ≈  tle_input.inclination         atol = 1e-4
        @test tle.mean_anomaly        ≈  tle_input.mean_anomaly        atol = 1e-4
        @test tle.mean_motion         ≈  tle_input.mean_motion         atol = 1e-7
        @test tle.raan                ≈  tle_input.raan                atol = 1e-4
        @test tle.argument_of_perigee ≈  tle_input.argument_of_perigee atol = 1e-4

        # Scenario 2: Geostationary Orbit
        # ==================================================================================

        tle_input = tle"""
            MOLNIYA 1-83
            1 21897U 92011A   06176.02341244 -.00001273  00000-0 -13525-3 0  3044
            2 21897  62.1749 198.0096 7421690 253.0462  20.1561  2.01269994104880
            """

        # Generate the osculating elements (TEME).
        sgp4d   = sgp4_init(tle_input )
        ret     = map(t -> sgp4!(sgp4d, t), 0:10:2880)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:10:2880) ./ 1440

        # Obtain the mean elements.
        tle, ~ = fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            atol                     = 5e-5,
            rtol                     = 5e-5,
            element_set_number       = 304,
            international_designator = "92011A",
            mean_elements_epoch      = vjd[begin],
            name                     = "MOLNIYA 1-83",
            revolution_number        = 10488,
            satellite_number         = 21897,
            max_iterations           = 7000,
            verbose                  = false,
        )

        # Compare.
        @test tle.classification           == tle_input.classification
        @test tle.element_set_number       == tle_input.element_set_number
        @test tle.epoch_year               == tle_input.epoch_year
        @test tle.international_designator == tle_input.international_designator
        @test tle.name                     == tle_input.name
        @test tle.revolution_number        == tle_input.revolution_number
        @test tle.satellite_number         == tle_input.satellite_number

        @test tle.bstar               ≈  tle_input.bstar               atol = 1e-6
        @test tle.eccentricity        ≈  tle_input.eccentricity        atol = 1e-7
        @test tle.epoch_day           ≈  tle_input.epoch_day           atol = 1e-8
        @test tle.inclination         ≈  tle_input.inclination         atol = 1e-4
        @test tle.mean_anomaly        ≈  tle_input.mean_anomaly        atol = 1e-4
        @test tle.mean_motion         ≈  tle_input.mean_motion         atol = 1e-7
        @test tle.raan                ≈  tle_input.raan                atol = 1e-4
        @test tle.argument_of_perigee ≈  tle_input.argument_of_perigee atol = 1e-4

        # Scenario 3: TLE epoch outside the interval
        # ==================================================================================

        tle_input = tle"""
            AMAZONIA 1
            1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
            2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
            """

        # Generate the osculating elements (TEME).
        sgp4d   = sgp4_init(tle_input)
        ret     = map(t -> sgp4!(sgp4d, t), 0:0.2:200)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:0.2:200) ./ 1440

        # Obtain the mean elements.
        tle, ~ = fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            atol                     = 1e-10,
            rtol                     = 1e-10,
            element_set_number       = 999,
            international_designator = "21015A",
            mean_elements_epoch      = vjd[begin] - 1,
            name                     = "AMAZONIA 1",
            revolution_number        = 3043,
            satellite_number         = 47699,
            max_iterations           = 1000,
            verbose                  = false,
        )

        # Compare.
        @test tle.classification           == tle_input.classification
        @test tle.element_set_number       == tle_input.element_set_number
        @test tle.epoch_year               == tle_input.epoch_year
        @test tle.international_designator == tle_input.international_designator
        @test tle.name                     == tle_input.name
        @test tle.revolution_number        == tle_input.revolution_number
        @test tle.satellite_number         == tle_input.satellite_number

        @test tle.bstar        ≈  tle_input.bstar               atol = 1e-6
        @test tle.eccentricity ≈  tle_input.eccentricity        atol = 1e-6
        @test tle.epoch_day    ≈  tle_input.epoch_day - 1       atol = 1e-8
        @test tle.inclination  ≈  tle_input.inclination         atol = 1e-4
        @test tle.raan         ≈  tle_input.raan - 0.9856002605 atol = 7e-3
    end

    @testset "TLE as Initial Guess" verbose = true begin
        # Scenario 1: Low Earth Orbit
        # ==================================================================================

        tle_input = tle"""
            AMAZONIA 1
            1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
            2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
            """

        # Generate the osculating elements (TEME).
        sgp4d   = sgp4_init(tle_input)
        ret     = map(t -> sgp4!(sgp4d, t), 0:0.2:200)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:0.2:200) ./ 1440

        # Obtain the mean elements.
        tle, ~ = fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            element_set_number       = 999,
            international_designator = "21015A",
            initial_guess            = tle_input,
            mean_elements_epoch      = vjd[begin],
            name                     = "AMAZONIA 1",
            revolution_number        = 3043,
            satellite_number         = 47699,
            max_iterations           = 10,
            verbose                  = false,
        )

        # Compare.
        @test tle.classification           == tle_input.classification
        @test tle.element_set_number       == tle_input.element_set_number
        @test tle.epoch_day                ≈  tle_input.epoch_day atol = 1e-8
        @test tle.epoch_year               == tle_input.epoch_year
        @test tle.international_designator == tle_input.international_designator
        @test tle.name                     == tle_input.name
        @test tle.revolution_number        == tle_input.revolution_number
        @test tle.satellite_number         == tle_input.satellite_number

        @test tle.bstar               ≈  tle_input.bstar               atol = 1e-6
        @test tle.eccentricity        ≈  tle_input.eccentricity        atol = 1e-7
        @test tle.epoch_day           ≈  tle_input.epoch_day           atol = 1e-8
        @test tle.inclination         ≈  tle_input.inclination         atol = 1e-4
        @test tle.mean_anomaly        ≈  tle_input.mean_anomaly        atol = 1e-4
        @test tle.mean_motion         ≈  tle_input.mean_motion         atol = 1e-7
        @test tle.raan                ≈  tle_input.raan                atol = 1e-4
        @test tle.argument_of_perigee ≈  tle_input.argument_of_perigee atol = 1e-4

        # Scenario 2: Geostationary Orbit
        # ==================================================================================

        tle_input = tle"""
            MOLNIYA 1-83
            1 21897U 92011A   06176.02341244 -.00001273  00000-0 -13525-3 0  3044
            2 21897  62.1749 198.0096 7421690 253.0462  20.1561  2.01269994104880
            """

        # Generate the osculating elements (TEME).
        sgp4d   = sgp4_init(tle_input )
        ret     = map(t -> sgp4!(sgp4d, t), 0:10:2880)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:10:2880) ./ 1440

        # Obtain the mean elements.
        tle, ~ = fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            element_set_number       = 304,
            international_designator = "92011A",
            initial_guess            = tle_input,
            mean_elements_epoch      = vjd[begin],
            name                     = "MOLNIYA 1-83",
            revolution_number        = 10488,
            satellite_number         = 21897,
            max_iterations           = 10,
            verbose                  = false,
        )

        # Compare.
        @test tle.classification           == tle_input.classification
        @test tle.element_set_number       == tle_input.element_set_number
        @test tle.epoch_day                ≈  tle_input.epoch_day atol = 1e-8
        @test tle.epoch_year               == tle_input.epoch_year
        @test tle.international_designator == tle_input.international_designator
        @test tle.name                     == tle_input.name
        @test tle.revolution_number        == tle_input.revolution_number
        @test tle.satellite_number         == tle_input.satellite_number

        @test tle.bstar               ≈  tle_input.bstar               atol = 1e-6
        @test tle.eccentricity        ≈  tle_input.eccentricity        atol = 1e-7
        @test tle.epoch_day           ≈  tle_input.epoch_day           atol = 1e-8
        @test tle.inclination         ≈  tle_input.inclination         atol = 1e-4
        @test tle.mean_anomaly        ≈  tle_input.mean_anomaly        atol = 1e-4
        @test tle.mean_motion         ≈  tle_input.mean_motion         atol = 1e-7
        @test tle.raan                ≈  tle_input.raan                atol = 1e-4
        @test tle.argument_of_perigee ≈  tle_input.argument_of_perigee atol = 1e-4

        # Scenario 3: TLE epoch outside the interval
        # ==================================================================================

        tle_input = tle"""
            AMAZONIA 1
            1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
            2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
            """

        # Generate the osculating elements (TEME).
        sgp4d   = sgp4_init(tle_input)
        ret     = map(t -> sgp4!(sgp4d, t), 0:0.2:200)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:0.2:200) ./ 1440

        # Obtain the mean elements.
        tle, ~ = fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            element_set_number       = 999,
            initial_guess            = tle_input,
            international_designator = "21015A",
            mean_elements_epoch      = vjd[begin] - 1,
            name                     = "AMAZONIA 1",
            revolution_number        = 3043,
            satellite_number         = 47699,
            max_iterations           = 1000,
            verbose                  = false,
        )

        # Compare.
        @test tle.classification           == tle_input.classification
        @test tle.element_set_number       == tle_input.element_set_number
        @test tle.epoch_year               == tle_input.epoch_year
        @test tle.international_designator == tle_input.international_designator
        @test tle.name                     == tle_input.name
        @test tle.revolution_number        == tle_input.revolution_number
        @test tle.satellite_number         == tle_input.satellite_number

        @test tle.bstar               ≈  tle_input.bstar               atol = 1e-6
        @test tle.eccentricity        ≈  tle_input.eccentricity        atol = 1e-6
        @test tle.epoch_day           ≈  tle_input.epoch_day - 1       atol = 1e-8
        @test tle.inclination         ≈  tle_input.inclination         atol = 1e-4
        @test tle.raan                ≈  tle_input.raan - 0.9856002605 atol = 7e-3
    end

    @testset "Errors" begin
        tle_input = tle"""
            AMAZONIA 1
            1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
            2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
            """

        sgp4d   = sgp4_init(tle_input)
        ret     = map(t -> sgp4!(sgp4d, t), 0:0.2:200)
        vr_teme = first.(ret)
        vv_teme = last.(ret)
        vjd     = sgp4d.epoch .+ (0:0.2:200) ./ 1440

        # Wrong dimensions in the input vectors
        # ==================================================================================

        @test_throws ArgumentError fit_sgp4_tle(vjd[1:end-1], vr_teme, vv_teme)
        @test_throws ArgumentError fit_sgp4_tle(vjd, vr_teme[1:end-1], vv_teme)
        @test_throws ArgumentError fit_sgp4_tle(vjd, vr_teme, vv_teme[1:end-1])

        # Wrong dimensions in the weight vector
        # ==================================================================================

        @test_throws ArgumentError fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            weight_vector = [1, 2, 3, 4, 5]
        )

        # Wrong dimensions in the initial guess vector
        # ==================================================================================

        @test_throws ArgumentError fit_sgp4_tle(
            vjd,
            vr_teme,
            vv_teme;
            initial_guess = [1, 2, 3, 4, 5, 6]
        )
    end
end

@testset "Updating SGP4 TLE epochs" verbose = true begin
    tle_input = tle"""
        AMAZONIA 1
        1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
        2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
        """

    tle = update_sgp4_tle_epoch(
        tle_input,
        DateTime("2021-09-28T11:40:12.955");
        verbose = false
    )

    # Compare.
    @test tle.classification           == tle_input.classification
    @test tle.element_set_number       == tle_input.element_set_number
    @test tle.epoch_year               == tle_input.epoch_year
    @test tle.international_designator == tle_input.international_designator
    @test tle.name                     == tle_input.name
    @test tle.revolution_number        == tle_input.revolution_number
    @test tle.satellite_number         == tle_input.satellite_number

    @test tle.bstar        ≈  tle_input.bstar               atol = 1e-6
    @test tle.eccentricity ≈  tle_input.eccentricity        atol = 1e-6
    @test tle.epoch_day    ≈  tle_input.epoch_day + 1       atol = 1e-8
    @test tle.inclination  ≈  tle_input.inclination         atol = 1e-4
    @test tle.raan         ≈  tle_input.raan + 0.9856002605 atol = 7e-3
end
