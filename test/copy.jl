## Description #############################################################################
#
# Test structure copying.
#
############################################################################################

@testset "Sgp4Propagator and Spg4DeepSpace" begin
    for sgp4c in (sgp4c_wgs84, sgp4c_wgs84_f32)

        sgp4d = sgp4_init(
            tle"""
            1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600
            2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119
            """;
            sgp4c = sgp4c
        )

        # We need to propagate the orbit to initialize all the internal terms related to the
        # deep space structure.
        sgp4!(sgp4d, 0)

        new_sgp4d = copy(sgp4d)

        @test typeof(new_sgp4d) == typeof(sgp4d)

        for f in fieldnames(typeof(sgp4d))
            if f != :sgp4ds
                @test getfield(new_sgp4d, f) == getfield(sgp4d, f)
            end
        end

        for f in fieldnames(typeof(sgp4d.sgp4ds))
            # Some fields in the deep space structure can be `NaN` because they was not
            # initialized for this orbit. We should skip them.
            if !isnan(getfield(sgp4d.sgp4ds, f))
                @test getfield(new_sgp4d.sgp4ds, f) == getfield(sgp4d.sgp4ds, f)
            end
        end

        # Test if both structures does not share the same memory region.
        new_sgp4d.epoch = 100
        @test new_sgp4d.epoch != sgp4d.epoch

        new_sgp4d.sgp4ds.atime = 123
        @test new_sgp4d.sgp4ds.atime != sgp4d.sgp4ds.atime
    end
end
