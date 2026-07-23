using SatelliteToolboxSgp4
using BenchmarkTools
using Printf
using Dates

# ==========================================================================================
#                    Finite-Difference vs ForwardDiff Jacobian Comparison
# ==========================================================================================

const REPORT_PATH = joinpath(@__DIR__, "finite_diff_vs_autodiff_report.md")

const ACCURACY_FIELDS = [
    :bstar,
    :eccentricity,
    :inclination,
    :raan,
    :argument_of_perigee,
    :mean_anomaly,
    :mean_motion,
]

function generate_osc_data(tle_input::TLE, time_range)
    sgp4d   = sgp4_init(tle_input)
    ret     = map(t -> sgp4!(sgp4d, t), time_range)
    vr_teme = first.(ret)
    vv_teme = last.(ret)
    vjd     = sgp4d.epoch .+ collect(time_range) ./ 1440
    return vjd, vr_teme, vv_teme
end

function run_scenario(sc)
    vjd, vr_teme, vv_teme = generate_osc_data(sc.tle_input, sc.time_range)
    kw = merge(sc.fit_kwargs, (; mean_elements_epoch = vjd[begin]))

    println("\n  Benchmarking FiniteDiffJacobian...")
    b_fd = @benchmark fit_sgp4_tle(
        $vjd, $vr_teme, $vv_teme; jacobian_method = FiniteDiffJacobian(), $kw...
    )
    display(b_fd)

    println("\n  Benchmarking ForwardDiffJacobian...")
    b_ad = @benchmark fit_sgp4_tle(
        $vjd, $vr_teme, $vv_teme; jacobian_method = ForwardDiffJacobian(), $kw...
    )
    display(b_ad)

    t_fd     = median(b_fd).time / 1e6
    t_ad     = median(b_ad).time / 1e6
    alloc_fd = median(b_fd).allocs
    alloc_ad = median(b_ad).allocs
    mem_fd   = median(b_fd).memory / 1024
    mem_ad   = median(b_ad).memory / 1024

    @printf("\n  Median: FD = %.1f ms, AD = %.1f ms\n\n", t_fd, t_ad)

    tle_fd, _ = fit_sgp4_tle(
        vjd, vr_teme, vv_teme; jacobian_method = FiniteDiffJacobian(), kw...
    )
    tle_ad, _ = fit_sgp4_tle(
        vjd, vr_teme, vv_teme; jacobian_method = ForwardDiffJacobian(), kw...
    )

    errors = Dict{Symbol, NamedTuple{(:ref, :fd, :ad), Tuple{Float64, Float64, Float64}}}()
    for f in ACCURACY_FIELDS
        r = Float64(getfield(sc.tle_input, f))
        e_fd = abs(Float64(getfield(tle_fd, f)) - r)
        e_ad = abs(Float64(getfield(tle_ad, f)) - r)
        errors[f] = (; ref = r, fd = e_fd, ad = e_ad)
    end

    return (; name = sc.name, t_fd, t_ad, alloc_fd, alloc_ad, mem_fd, mem_ad, errors)
end

# ------------------------------------------------------------------------------------------
# TLE definitions
# ------------------------------------------------------------------------------------------

tle_amazonia = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436
    """

tle_molniya = tle"""
    MOLNIYA 1-83
    1 21897U 92011A   06176.02341244 -.00001273  00000-0 -13525-3 0  3044
    2 21897  62.1749 198.0096 7421690 253.0462  20.1561  2.01269994104880
    """

# ------------------------------------------------------------------------------------------
# Tolerance tiers
# ------------------------------------------------------------------------------------------

tolerance_tiers = [
    (
        label = "Default Tolerances",
        scenarios = [
            (
                name       = "LEO (AMAZONIA 1) — no initial guess",
                tle_input  = tle_amazonia,
                time_range = 0:0.2:200,
                fit_kwargs = (atol = 1e-10, rtol = 1e-10, max_iterations = 1000, verbose = false, mean_elements_epoch = nothing),
            ),
            (
                name       = "HEO (MOLNIYA 1-83) — no initial guess",
                tle_input  = tle_molniya,
                time_range = 0:10:2880,
                fit_kwargs = (atol = 5e-5, rtol = 5e-5, max_iterations = 7000, verbose = false, mean_elements_epoch = nothing),
            ),
            (
                name       = "LEO (AMAZONIA 1) — TLE initial guess",
                tle_input  = tle_amazonia,
                time_range = 0:0.2:200,
                fit_kwargs = (max_iterations = 10, verbose = false, initial_guess = tle_amazonia, mean_elements_epoch = nothing),
            ),
            (
                name       = "HEO (MOLNIYA 1-83) — TLE initial guess",
                tle_input  = tle_molniya,
                time_range = 0:10:2880,
                fit_kwargs = (max_iterations = 10, verbose = false, initial_guess = tle_molniya, mean_elements_epoch = nothing),
            ),
        ],
    ),
    (
        label = "Tight Tolerances (atol = rtol = 1e-14)",
        scenarios = [
            (
                name       = "LEO (AMAZONIA 1) — no initial guess",
                tle_input  = tle_amazonia,
                time_range = 0:0.2:200,
                fit_kwargs = (atol = 1e-14, rtol = 1e-14, max_iterations = 5000, verbose = false, mean_elements_epoch = nothing),
            ),
            (
                name       = "HEO (MOLNIYA 1-83) — no initial guess",
                tle_input  = tle_molniya,
                time_range = 0:10:2880,
                fit_kwargs = (atol = 1e-14, rtol = 1e-14, max_iterations = 10000, verbose = false, mean_elements_epoch = nothing),
            ),
            (
                name       = "LEO (AMAZONIA 1) — TLE initial guess",
                tle_input  = tle_amazonia,
                time_range = 0:0.2:200,
                fit_kwargs = (atol = 1e-14, rtol = 1e-14, max_iterations = 5000, verbose = false, initial_guess = tle_amazonia, mean_elements_epoch = nothing),
            ),
            (
                name       = "HEO (MOLNIYA 1-83) — TLE initial guess",
                tle_input  = tle_molniya,
                time_range = 0:10:2880,
                fit_kwargs = (atol = 1e-14, rtol = 1e-14, max_iterations = 10000, verbose = false, initial_guess = tle_molniya, mean_elements_epoch = nothing),
            ),
        ],
    ),
]

# ------------------------------------------------------------------------------------------
# Warmup pass — compile all code paths across every scenario
# ------------------------------------------------------------------------------------------

println("Warmup pass: compiling all code paths...")

for tier in tolerance_tiers
    for sc in tier.scenarios
        vjd, vr_teme, vv_teme = generate_osc_data(sc.tle_input, sc.time_range)
        kw = merge(sc.fit_kwargs, (; mean_elements_epoch = vjd[begin]))
        fit_sgp4_tle(vjd, vr_teme, vv_teme; jacobian_method = FiniteDiffJacobian(), kw...)
        fit_sgp4_tle(vjd, vr_teme, vv_teme; jacobian_method = ForwardDiffJacobian(), kw...)
    end
end

println("Warmup complete.\n")

# ------------------------------------------------------------------------------------------
# Run all tiers
# ------------------------------------------------------------------------------------------

all_tier_results = []

for tier in tolerance_tiers
    println("#" ^ 96)
    println("# $(tier.label)")
    println("#" ^ 96)

    tier_results = []
    for (i, sc) in enumerate(tier.scenarios)
        println("\n" * "=" ^ 96)
        println("Scenario $i: $(sc.name)")
        println("=" ^ 96)
        push!(tier_results, run_scenario(sc))
    end
    push!(all_tier_results, (; label = tier.label, results = tier_results))
end

# ------------------------------------------------------------------------------------------
# Generate Markdown report
# ------------------------------------------------------------------------------------------

function write_tier(io, label, results)
    println(io, "## $label")
    println(io)

    println(io, "### Performance")
    println(io)
    println(
        io,
        "| Scenario | FD Median (ms) | AD Median (ms) | Speedup | FD Allocs | AD Allocs | FD Mem (KiB) | AD Mem (KiB) |",
    )
    println(
        io,
        "|:---------|---------------:|---------------:|--------:|----------:|----------:|-------------:|-------------:|",
    )

    for r in results
        ratio = r.t_fd / r.t_ad
        arrow = ratio > 1 ? "AD" : "FD"
        @printf(
            io,
            "| %s | %.1f | %.1f | %.2fx %s | %d | %d | %.0f | %.0f |\n",
            r.name,
            r.t_fd,
            r.t_ad,
            max(ratio, 1/ratio),
            arrow,
            r.alloc_fd,
            r.alloc_ad,
            r.mem_fd,
            r.mem_ad
        )
    end
    println(io)

    println(io, "### Accuracy (Absolute Error vs Reference TLE)")
    println(io)

    for r in results
        println(io, "#### $(r.name)")
        println(io)
        println(io, "| Field | Reference | FD Error | AD Error | Winner |")
        println(io, "|:------|----------:|---------:|---------:|:------:|")

        for f in ACCURACY_FIELDS
            e = r.errors[f]
            winner = e.fd < e.ad ? "FD" : (e.ad < e.fd ? "AD" : "Tie")
            @printf(
                io, "| `%s` | %.12e | %.4e | %.4e | %s |\n", f, e.ref, e.fd, e.ad, winner
            )
        end
        println(io)
    end
end

open(REPORT_PATH, "w") do io
    println(io, "# Finite-Difference vs ForwardDiff Jacobian — Benchmark Report")
    println(io)
    println(io, "> Auto-generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "> Julia $(VERSION) — $(Sys.MACHINE)")
    println(io)

    for tier in all_tier_results
        write_tier(io, tier.label, tier.results)
        println(io, "---")
        println(io)
    end
end

println("=" ^ 96)
println("Report written to: $REPORT_PATH")
println("=" ^ 96)
