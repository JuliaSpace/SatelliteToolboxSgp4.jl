using Lux, Optimisers, Zygote
using StaticArrays
using LinearAlgebra

const _tle_amazonia = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
    """

# ------------------------------------------------------------------------------------------
# Generate a small reference dataset from plain SGP4.
# ------------------------------------------------------------------------------------------

function _generate_ref_data(tle, times)
    sgp4d = sgp4_init(tle)
    epoch = sgp4d.epoch
    vjd   = epoch .+ times ./ 1440.0

    vr = Vector{SVector{3, Float64}}(undef, length(times))
    vv = Vector{SVector{3, Float64}}(undef, length(times))
    for (k, t) in enumerate(times)
        r, v = sgp4!(sgp4d, t)
        vr[k] = r
        vv[k] = v
    end
    return vjd, vr, vv
end

times_train = collect(0.0:10.0:1440.0)
vjd, vr, vv = _generate_ref_data(_tle_amazonia, times_train)

# ==========================================================================================

@testset "Zero model matches plain SGP4" begin
    mlsgp4d = ml_dsgp4_init(_tle_amazonia)
    sgp4d   = sgp4_init(_tle_amazonia)

    for Δt in [0.0, 60.0, 360.0, 1440.0]
        r_ml, v_ml = ml_dsgp4!(mlsgp4d, Δt)
        r_sg, v_sg = sgp4!(sgp4d, Δt)

        @test r_ml isa SVector{3, Float64}
        @test v_ml isa SVector{3, Float64}
        @test r_ml ≈ r_sg rtol = 1e-14
        @test v_ml ≈ v_sg rtol = 1e-14
    end
end

@testset "ml_dsgp4 convenience API" begin
    r, v, mlsgp4d = ml_dsgp4(120.0, _tle_amazonia)

    @test r isa SVector{3, Float64}
    @test v isa SVector{3, Float64}
    @test norm(r) > 6000.0

    r2, v2 = ml_dsgp4!(mlsgp4d, 120.0)
    @test r2 ≈ r rtol = 1e-14
    @test v2 ≈ v rtol = 1e-14
end


@testset "Save / load round trip" begin
    config = MLdSGP4Config(hidden_layers = [16, 16])

    model = ml_dsgp4_train(
        _tle_amazonia, vjd, vr, vv;
        config  = config,
        epochs  = 3,
        verbose = false,
    )

    path = tempname() * ".bin"
    try
        ml_dsgp4_save(path, model)
        loaded = ml_dsgp4_load(path)

        @test loaded.config.hidden_layers == model.config.hidden_layers
        @test loaded.α ≈ model.α
        @test loaded.β ≈ model.β

        mlsgp4d_orig   = ml_dsgp4_init(_tle_amazonia; model = model)
        mlsgp4d_loaded = ml_dsgp4_init(_tle_amazonia; model = loaded)

        for Δt in [0.0, 60.0, 720.0]
            r1, v1 = ml_dsgp4!(mlsgp4d_orig, Δt)
            r2, v2 = ml_dsgp4!(mlsgp4d_loaded, Δt)
            @test r1 ≈ r2 atol = 1e-12
            @test v1 ≈ v2 atol = 1e-14
        end
    finally
        rm(path; force = true)
    end
end

@testset "ml_dsgp4 with trained model" begin
    config = MLdSGP4Config(hidden_layers = [16, 16])

    model = ml_dsgp4_train(
        _tle_amazonia, vjd, vr, vv;
        config  = config,
        epochs  = 3,
        verbose = false,
    )

    r, v, mlsgp4d = ml_dsgp4(120.0, _tle_amazonia; model = model)

    r2, v2 = ml_dsgp4!(mlsgp4d, 120.0)
    @test r ≈ r2 atol = 1e-12
    @test v ≈ v2 atol = 1e-14
end