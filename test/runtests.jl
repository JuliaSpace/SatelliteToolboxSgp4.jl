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
