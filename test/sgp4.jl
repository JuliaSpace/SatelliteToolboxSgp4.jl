## Description #############################################################################
#
# Test SGP4 algorithm. All tests are based on [1].
#
## References ##############################################################################
#
# [1] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting Spacetrack
#     Report #3: Rev1. AIAA.
#
############################################################################################

@testset "Constructors" begin
    sgp4ds = SatelliteToolboxSgp4.Sgp4DeepSpace{Float64}(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        false,
        false,
    )

    sgp4c = Sgp4Propagator{Float64, Float64}(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        :sdp4,
        sgp4c_wgs84,
        sgp4ds,
    )

    # Some random tests.
    @test sgp4c.epoch == 0.0
    @test sgp4c.a_k == 0.0
    @test sgp4c.∂M == 0.0
    @test sgp4c.algorithm == :sdp4
    @test sgp4c.sgp4ds == sgp4ds
end

@testset "Tests from the Paper AIAA 2006-6753" verbose = true begin
    # Read all TLEs that will be used to test.
    tles = read_tles_from_file("./sgp4_tests/sgp4_tests.tle")

    @testset "Default" begin
        for tle in tles
            filename = @sprintf(
                "./sgp4_tests/aiaa-2006-6753/sgp4_tle_%d_result.txt", tle.satellite_number
            )
            SGP4_results = readdlm(filename; comments = true)

            # Initialize the orbit propagator.
            sgp4d = sgp4_init(tle; sgp4c = sgp4c_wgs72)

            t = SGP4_results[:, 1]

            @inbounds for k in 1:length(t)

                # Propagate the orbit.
                r_teme, v_teme = sgp4!(sgp4d, t[k])

                # Assemble the result vector.
                st_sgp4_result = vcat(t[k], r_teme, v_teme)

                # Compare the values.
                @test t[k] == SGP4_results[k, 1]
                @test r_teme[1] ≈ SGP4_results[k, 2] atol=1e-8
                @test r_teme[2] ≈ SGP4_results[k, 3] atol=1e-8
                @test r_teme[3] ≈ SGP4_results[k, 4] atol=1e-8
                @test v_teme[1] ≈ SGP4_results[k, 5] atol=1e-9
                @test v_teme[2] ≈ SGP4_results[k, 6] atol=1e-9
                @test v_teme[3] ≈ SGP4_results[k, 7] atol=1e-9
            end
        end
    end

    @testset "In-place Initialization" begin
        # First, we create a dummy SGP4 structure but with the correct constants and epoch
        # type.
        sgp4d = sgp4_init(0.0, 0, 0, 0, 0, 0, 0, 0; sgp4c = sgp4c_wgs72)

        for tle in tles
            filename = @sprintf(
                "./sgp4_tests/aiaa-2006-6753/sgp4_tle_%d_result.txt", tle.satellite_number
            )
            SGP4_results = readdlm(filename; comments = true)

            # Initialize the orbit propagator.
            sgp4_init!(sgp4d, tle)

            t = SGP4_results[:, 1]

            @inbounds for k in 1:length(t)

                # Propagate the orbit.
                r_teme, v_teme = sgp4!(sgp4d, t[k])

                # Assemble the result vector.
                st_sgp4_result = vcat(t[k], r_teme, v_teme)

                # Compare the values.
                @test t[k] == SGP4_results[k, 1]
                @test r_teme[1] ≈ SGP4_results[k, 2] atol=1e-8
                @test r_teme[2] ≈ SGP4_results[k, 3] atol=1e-8
                @test r_teme[3] ≈ SGP4_results[k, 4] atol=1e-8
                @test v_teme[1] ≈ SGP4_results[k, 5] atol=1e-9
                @test v_teme[2] ≈ SGP4_results[k, 6] atol=1e-9
                @test v_teme[3] ≈ SGP4_results[k, 7] atol=1e-9
            end
        end
    end

    @testset "Simultaneous Creation and Propagation" begin
        # Read all TLEs that will be used to test.
        tles = read_tles_from_file("./sgp4_tests/sgp4_tests.tle")

        for tle in tles
            filename = @sprintf(
                "./sgp4_tests/aiaa-2006-6753/sgp4_tle_%d_result.txt", tle.satellite_number
            )
            SGP4_results = readdlm(filename; comments = true)
            t = SGP4_results[:, 1]

            # Initialize the orbit propagator.
            r_teme, v_teme, sgp4d = sgp4(t[end], tle; sgp4c = sgp4c_wgs72)

            # We test just the final instant to save computational burden.
            @test t[end] == SGP4_results[end, 1]
            @test r_teme[1] ≈ SGP4_results[end, 2] atol=1e-8
            @test r_teme[2] ≈ SGP4_results[end, 3] atol=1e-8
            @test r_teme[3] ≈ SGP4_results[end, 4] atol=1e-8
            @test v_teme[1] ≈ SGP4_results[end, 5] atol=1e-9
            @test v_teme[2] ≈ SGP4_results[end, 6] atol=1e-9
            @test v_teme[3] ≈ SGP4_results[end, 7] atol=1e-9
        end
    end
end

@testset "Errors" begin
    tle = tle"""
       AMAZONIA 1
       1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
       2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652
       """

    sgp4d = sgp4_init(tle)
    sgp4d.algorithm = :any
    @test_throws ErrorException sgp4!(sgp4d, 10)
end

@testset "Retrograde Equatorial Orbit" begin
    # An inclination of 180° leads to a division by zero in the long-period periodic term
    # unless the denominator is clamped as in Vallado's implementation.
    n₀ = 15 * 2π / 1440
    sgp4d_pro = sgp4_init(2.46e6, n₀, 0.001, 0.0, 0.0, 0.0, 0.0, 1e-4)
    sgp4d_ret = sgp4_init(2.46e6, n₀, 0.001, Float64(π), 0.0, 0.0, 0.0, 1e-4)

    r_pro, v_pro = sgp4!(sgp4d_pro, 100.0)
    r_ret, v_ret = sgp4!(sgp4d_ret, 100.0)

    @test !any(isnan, r_ret)
    @test !any(isnan, v_ret)

    # The retrograde orbit must mirror the prograde one about the x-axis.
    @test r_ret[1] ≈ +r_pro[1] atol = 1e-6
    @test r_ret[2] ≈ -r_pro[2] atol = 1e-6
    @test r_ret[3] ≈ 0 atol = 1e-6
    @test v_ret[1] ≈ +v_pro[1] atol = 1e-9
    @test v_ret[2] ≈ -v_pro[2] atol = 1e-9
    @test v_ret[3] ≈ 0 atol = 1e-9
end
