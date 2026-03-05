using Test

using Dates
using DelimitedFiles
using Printf
using SatelliteToolboxTle
using SatelliteToolboxSgp4

@testset "SGP4 Propagator" verbose = true begin
    include("./sgp4.jl")
end

@testset "SGP4 TLEs" verbose = true begin
    include("./tle.jl")
end

@testset "Copy Structures" verbose = true begin
    include("./copy.jl")
end

if isempty(VERSION.prerelease)
    using Pkg

    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")

    using JET
    using AllocCheck
    using Aqua

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end
else
    @warn "Performance checks not guaranteed to work on julia-nightly, skipping tests"
end
