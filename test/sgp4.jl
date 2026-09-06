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
        0,
        0,
        0,
        :sdp4,
        SGP4C_WGS84,
        sgp4ds,
    )

    # Some random tests.
    @test sgp4c.epoch == 0.0
    @test sgp4c.a_k == 0.0
    @test sgp4c.∂M == 0.0
    @test sgp4c.algorithm == :sdp4
    @test sgp4c.sgp4ds == sgp4ds

    # The zero constructor of the deep space structure must match the positional one.
    @test SatelliteToolboxSgp4.Sgp4DeepSpace{Float64}() == sgp4ds

    # Constructor that only sets the gravitational constants.
    sgp4d = Sgp4Propagator{Float64}(SGP4C_WGS72)

    @test sgp4d isa Sgp4Propagator{Float64, Float64}
    @test sgp4d.sgp4c === SGP4C_WGS72
    @test sgp4d.sgp4ds == SatelliteToolboxSgp4.Sgp4DeepSpace{Float64}()

    sgp4d = Sgp4Propagator{Float32}(Sgp4Constants{Float32}(SGP4C_WGS72))

    @test sgp4d isa Sgp4Propagator{Float32, Float32}
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
            sgp4d = sgp4_init(tle; sgp4c = SGP4C_WGS72)

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
        sgp4d = sgp4_init(0.0, 0, 0, 0, 0, 0, 0, 0; sgp4c = SGP4C_WGS72)

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
            r_teme, v_teme, sgp4d = sgp4(t[end], tle; sgp4c = SGP4C_WGS72)

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

@testset "Number Type" begin
    jd = 2.46e6

    # The number type of the propagator must always be the number type of the constants,
    # regardless of the input types.
    sgp4d = sgp4_init(jd, 0.06, 0.001, 1.7, 2.8, 2.4, 3.9, 4e-5; sgp4c = SGP4C_WGS72)
    @test sgp4d isa Sgp4Propagator{Float64, Float64}

    sgp4c_f32 = Sgp4Constants{Float32}(SGP4C_WGS72)

    sgp4d = sgp4_init(jd, 0.06, 0.001, 1.7, 2.8, 2.4, 3.9, 4e-5; sgp4c = sgp4c_f32)
    @test sgp4d isa Sgp4Propagator{Float64, Float32}

    sgp4d = sgp4_init(jd, 0.06, 0.001, 1.7, 2.8, 2.4, 3.9, 0; sgp4c = sgp4c_f32)
    @test sgp4d isa Sgp4Propagator{Float64, Float32}

    sgp4d = sgp4_init(jd, 0.06f0, 0.001f0, 1.7f0, 2.8f0, 2.4f0, 3.9f0, 4.0f-5)
    @test sgp4d isa Sgp4Propagator{Float64, Float64}

    r, v, sgp4d = sgp4(10, jd, 0.06, 0.001, 1.7, 2.8, 2.4, 3.9, 4e-5; sgp4c = sgp4c_f32)
    @test sgp4d isa Sgp4Propagator{Float64, Float32}
    @test eltype(r) === Float32
    @test eltype(v) === Float32

    # The epoch type is independent from the propagation type.
    sgp4d = sgp4_init(2.46f6, 0.06, 0.001, 1.7, 2.8, 2.4, 3.9, 4e-5)
    @test sgp4d isa Sgp4Propagator{Float32, Float64}
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

@testset "Show" begin
    tle = tle"""
        AMAZONIA 1
        1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
        2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""

    sgp4d = sgp4_init(tle)
    sgp4!(sgp4d, 10.0)

    # == Compact ===========================================================================

    @test repr(sgp4d) ==
        "Sgp4Propagator{Float64, Float64} (SGP4): Epoch = 2.46003e6 (2023-03-24T16:28:40.388)"

    # == Multi-line ========================================================================

    expected = join(
        (
            "Sgp4Propagator{Float64, Float64} (SGP4):",
            "  Epoch            : 2.46003e6 (2023-03-24T16:28:40.388)",
            "  Last Propagation : 10.0 min",
            "  ├─ Mean Elements",
            "  │    Semi-Major Axis    : 7133.946314 km",
            "  │    Mean Motion        : 14.40814394 rev/day",
            "  │    Eccentricity       : 0.0001247",
            "  │    Inclination        : 98.4304°",
            "  │    RA of Asc. Node    : 162.1097°",
            "  │    Arg. of Pericenter : 136.2017°",
            "  │    Mean Anomaly       : 223.9283°",
            "  │    B*                 : 4.3e-5 1/er",
            "  └─ Constants",
            "       R₀  : 6378.137 km",
            "       XKE : 0.07436685317 er^(3/2)/min",
            "       J₂  : 0.001082629989",
            "       J₃  : -2.53215306e-6",
            "       J₄  : -1.61098761e-6",
        ),
        '\n',
    )

    str = sprint(show, MIME("text/plain"), sgp4d)
    @test str == expected

    # The body can be printed under another header.
    str = sprint(SatelliteToolboxBase.print_tree_body, sgp4d)
    @test str == expected[(length("Sgp4Propagator{Float64, Float64} (SGP4):\n") + 1):end]

    # The decorations must not change the text when the output supports colors.
    str_color = sprint(show, MIME("text/plain"), sgp4d; context = :color => true)

    @test occursin("\e[1m", str_color)
    @test occursin("\e[90mkm\e[39m", str_color)
    @test replace(str_color, r"\e\[[0-9;]*m" => "") == expected

    # == Other Algorithms ==================================================================

    tle_ds = tle"""
        1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
        2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656"""

    @test occursin("(SDP4): ", repr(sgp4_init(tle_ds)))

    tle_lp = tle"""
        1 28872U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1534
        2 28872  96.4736 157.9986 0303955 244.0492 110.6523 16.46015938 10708"""

    @test occursin("(SGP4 (low perigee)): ", repr(sgp4_init(tle_lp)))

    # == Uninitialized Propagator ==========================================================

    sgp4d = Sgp4Propagator{Float64}(SGP4C_WGS84)

    @test repr(sgp4d) == "Sgp4Propagator{Float64, Float64} (not initialized)"
    expected = "Sgp4Propagator{Float64, Float64}:\n  Status : not initialized"
    @test sprint(show, MIME("text/plain"), sgp4d) == expected
end
