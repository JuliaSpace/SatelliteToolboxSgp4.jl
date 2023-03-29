# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Test SGP4 algorithm. All tests are based on [1].
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Test

using DelimitedFiles
using Printf
using SatelliteToolboxTle
using SatelliteToolboxSgp4

@testset "Tests from the paper AIAA 2006-6753" begin
    # Read all TLEs that will be used to test.
    tles = read_tles_from_file("./sgp4_tests/sgp4_tests.tle")

    for tle in tles
        filename = @sprintf(
            "./sgp4_tests/aiaa-2006-6753/sgp4_tle_%d_result.txt",
            tle.satellite_number
        )
        SGP4_results = readdlm(filename; comments = true)

        # Initialize the orbit propagator.
        sgp4d = sgp4_init(tle; sgp4c = sgp4c_wgs72)

        t = SGP4_results[:,1]
        @inbounds for k = 1:length(t)

            # Propagate the orbit.
            r_TEME, v_TEME = sgp4!(sgp4d, t[k])

            # Assemble the result vector.
            st_sgp4_result = vcat(t[k], r_TEME, v_TEME)

            # Compare the values.
            @test t[k]      == SGP4_results[k,1]
            @test r_TEME[1] ≈  SGP4_results[k,2] atol=1e-8
            @test r_TEME[2] ≈  SGP4_results[k,3] atol=1e-8
            @test r_TEME[3] ≈  SGP4_results[k,4] atol=1e-8
            @test v_TEME[1] ≈  SGP4_results[k,5] atol=1e-9
            @test v_TEME[2] ≈  SGP4_results[k,6] atol=1e-9
            @test v_TEME[3] ≈  SGP4_results[k,7] atol=1e-9
        end
    end
end
