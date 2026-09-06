## Description #############################################################################
#
# Code quality tests using Aqua.jl and JET.jl.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxSgp4; ambiguities = (recursive = false))
end

if isempty(VERSION.prerelease)
    @testset "JET.jl" verbose = true begin
        # == Runtime Errors ================================================================

        @testset "Package Analysis" begin
            JET.test_package(
                SatelliteToolboxSgp4;
                toplevel_logger = nothing,
                target_modules = (SatelliteToolboxSgp4,),
            )
        end

        # == Optimization Analysis of the Propagation Kernel ===============================

        @testset "Propagation Kernel" begin
            tle = tle"""
                AMAZONIA 1
                1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
                2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""

            # Deep space orbit to exercise the SDP4 branches.
            tle_ds = tle"""
                MOLNIYA 2-14
                1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
                2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656"""

            for (sgp4c, T) in ((sgp4c_wgs84, Float64), (sgp4c_wgs84_f32, Float32))
                sgp4d = sgp4_init(tle; sgp4c = sgp4c)

                @test_opt target_modules = (SatelliteToolboxSgp4,) sgp4_init!(
                    sgp4d,
                    tle_epoch(tle),
                    T(0.06),
                    T(0.001),
                    T(1.7),
                    T(2.8),
                    T(2.4),
                    T(3.9),
                    T(4.3e-5),
                )
                @test_opt target_modules = (SatelliteToolboxSgp4,) sgp4_init!(sgp4d, tle)
                @test_opt target_modules = (SatelliteToolboxSgp4,) sgp4!(sgp4d, T(10))

                sgp4_init!(sgp4d, tle_ds)
                @test_opt target_modules = (SatelliteToolboxSgp4,) sgp4!(sgp4d, T(10))
            end
        end
    end
else
    @warn "JET.jl tests are skipped on Julia pre-release versions."
end
