using Test

using Aqua
using Dates
using DelimitedFiles
using JET
using Printf
using SatelliteToolboxOrbitDataMessages
using SatelliteToolboxTle
using StaticArrays
using SatelliteToolboxSgp4

@testset "SGP4 Propagator" verbose = true begin
    include("./sgp4.jl")
end

@testset "SGP4 TLEs" verbose = true begin
    include("./tle.jl")
end

@testset "SGP4 OMMs" verbose = true begin
    include("./omm.jl")
end

@testset "Copy Structures" verbose = true begin
    include("./copy.jl")
end

@testset "Quality Tests" verbose = true begin
    include("./quality.jl")
end

if isempty(VERSION.prerelease)
    using Pkg

    Pkg.add("AllocCheck")

    using AllocCheck

    using ForwardDiff

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end
else
    @warn "Performance checks not guaranteed to work on julia-nightly, skipping tests"
end
