# Repository Guide

## Package Structure

- Single Julia package (`SatelliteToolboxSgp4`): the SGP4/SDP4 orbit propagator. Requires Julia 1.10 or newer (`[compat] julia = "1.10, 1.11, 1.12"`).
- `src/SatelliteToolboxSgp4.jl` is the module entrypoint and controls source include order: `types.jl` → `copy.jl` → `sgp4_model.jl` → `tle.jl` → `precompile.jl`. New code must respect this order when referencing symbols across files.
- `@reexport using SatelliteToolboxBase` and `@reexport using SatelliteToolboxTle` re-export those packages' public API, so their symbols are part of this package's public surface.
- `src/precompile.jl` holds a `PrecompileTools.@compile_workload` covering orbit propagation, TLE fitting, TLE epoch update, and structure copying (Float64 and Float32); update it when the precompiled public API changes.
- Tests are wired from `test/runtests.jl`, which `include`s `sgp4.jl`, `tle.jl`, and `copy.jl` unconditionally, then `performance.jl` only on non-prerelease Julia. `src/types.jl` and `src/precompile.jl` have no dedicated test file; `test/performance.jl` has no `src/` counterpart.
- Test-only dependencies are declared via `[extras]` + `[targets]` in `Project.toml` (Test, DelimitedFiles, Pkg, Printf — all stdlibs). There is no `test/Project.toml`.
- No `Manifest.toml` is committed; the dependency tree is not pinned.
- No package extensions (`[weakdeps]`/`[extensions]` absent; no `ext/`).

## Commands

- Instantiate: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`. A first run precompiles dependencies and can take minutes while printing little; use generous timeouts and do not assume a hang.
- Focused test file (core tests): `julia --project=. -e 'using SatelliteToolboxSgp4, SatelliteToolboxTle, Test, Dates, DelimitedFiles, Printf; include("test/sgp4.jl")'` — swap the path for `test/tle.jl` or `test/copy.jl`. The `using` list mirrors `test/runtests.jl`; the included test files rely on those names being in scope.
- `test/performance.jl` is not runnable via a focused include: `runtests.jl` `Pkg.add`s JET, AllocCheck, and Aqua at runtime and gates the whole block on non-prerelease Julia. Run it through the full `Pkg.test()`.
- There is no test-name selector; `runtests.jl` does not read `ARGS` and does not use TestItemRunner/ReTestItems.
- CI (`.github/workflows/ci.yml`) tests Julia 1.10 (oldest supported) and latest stable 1.x on Ubuntu (x64), macOS (arm64), and Windows (x64); it builds via `julia-actions/julia-buildpkg` before `julia-actions/julia-runtest`, then uploads coverage to Codecov. `ci-nightly.yml` repeats the matrix on Julia nightly. `Pkg.build()` is a no-op (no `deps/build.jl`). CompatHelper and TagBot workflows also run.

## Code Style

- No formatter is configured (no `.JuliaFormatter.toml`); CI does not run a format check. Do not invent a formatting step.

## Behavioral Constraints

- SGP4/SDP4 propagation correctness is validated against Vallado et al., "Revisiting Spacetrack Report #3" (see the `test/sgp4.jl` header). Preserve numerical fidelity when touching the algorithm in `src/sgp4_model.jl`.
- Hot paths are allocation-checked in `test/performance.jl`: `sgp4!`, `sgp4_init!`, and the internal TLE-fitting Jacobian functions (`_init_sgp4_with_state_vector!`, `_sgp4_jacobian`, `_sgp4_fwd_jacobian_eval`) must allocate zero. `fit_sgp4_tle!` has explicit regression bounds (≤45 allocs with `FiniteDiffJacobian`, ≤47 with `ForwardDiffJacobian`). Do not introduce allocations in those paths.
- JET and AllocCheck tests are skipped on Julia 1.12+ (see `test/performance.jl`); Aqua always runs on non-prerelease Julia.
- New tests follow the `@testset "Name" verbose = true begin ... end` pattern used in `test/runtests.jl`; match it when adding coverage.

## Not Configured

- No formatter, linter, pre-commit hooks, or generated docs are configured; do not invent them from README badges or Julia convention.
