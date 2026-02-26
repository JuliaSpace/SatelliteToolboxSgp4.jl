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
    Pkg.add("Random")
    Pkg.add("Serialization")
    Pkg.add("Lux")
    Pkg.add("Optimisers")
    Pkg.add("Zygote")

    using JET
    using AllocCheck
    using Aqua
    using Random
    using Serialization
    using Lux
    using Optimisers
    using Zygote

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end

    @testset "Lux Extension (ML-∂SGP4)" verbose = true begin
        include("./extension.jl")
    end
else
    @warn "Performance checks and extensions not guaranteed to work on julia-nightly, skipping tests"
end
