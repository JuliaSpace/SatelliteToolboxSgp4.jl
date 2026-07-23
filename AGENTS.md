# AGENTS.md

Instructions for coding agents working in this repository. `CLAUDE.md` imports this file via `@AGENTS.md`.

## Package Structure

`SatelliteToolboxSgp4` is a pure-Julia implementation of the SGP4/SDP4 orbit propagator. It is part of the SatelliteToolbox ecosystem and re-exports `SatelliteToolboxBase` and `SatelliteToolboxTle`.

- **Julia compat:** `[compat] julia = "1.10, 1.11, 1.12"` — supported on Julia 1.10, 1.11, and 1.12 only (not an open `^1.10` range; do not assume 1.13+ works). Nightly is exercised in CI but is not a supported target.
- **Module entrypoint:** `src/SatelliteToolboxSgp4.jl`. Include order is fixed and load-order-sensitive: `types.jl` → `copy.jl` → `sgp4_model.jl` → `tle.jl` → `precompile.jl`. New code must respect this order; symbols defined in an earlier file are visible to later ones, not vice versa.
- **Re-exports:** `@reexport using SatelliteToolboxBase` and `@reexport using SatelliteToolboxTle` — public API of those packages is part of this package's public surface.
- **Runtime deps (`[deps]`):** Crayons, Dates, ForwardDiff, LinearAlgebra, PrecompileTools, Printf, Reexport, SatelliteToolboxBase, SatelliteToolboxTle, StaticArrays. `PrecompileTools` is used for precompile workloads in `src/precompile.jl`.
- **No package extensions:** no `[weakdeps]`, `[extensions]`, or `ext/` directory.
- **No build script:** no `deps/build.jl`; `Pkg.build()` is a no-op.
- **Test deps:** declared via `[extras]` + `[targets] test = ["Test", "DelimitedFiles", "Pkg", "Printf"]` (no `test/Project.toml`). The performance test additionally `Pkg.add`s `JET`, `AllocCheck`, and `Aqua` at runtime inside `test/runtests.jl`, so a first test run on a non-prerelease Julia will hit the network.
- **Test wiring (`test/runtests.jl`):** three unconditional `@testset`s include `sgp4.jl`, `tle.jl`, `copy.jl`. A fourth "Performance Tests" `@testset` includes `performance.jl` only when `isempty(VERSION.prerelease)` (i.e. skipped on nightly); it runs Aqua, JET (skipped on Julia 1.12+), and AllocCheck (skipped on Julia 1.12+).
- **Test ↔ src mapping:** `test/copy.jl` ↔ `src/copy.jl`; `test/sgp4.jl` ↔ `src/sgp4_model.jl`; `test/tle.jl` ↔ `src/tle.jl`. `src/types.jl` and `src/precompile.jl` have no dedicated test file (covered indirectly). Test fixtures live in `test/sgp4_tests/` (AIAA 2006-6753 reference TLEs and expected results).
- **Examples:** `examples/` is a separate environment with its own `Project.toml` (`[sources] SatelliteToolboxSgp4 = {path = ".."}`); run scripts there with `julia --project=examples`.
- **Manifest:** `Manifest.toml` is gitignored and not committed.

## Commands

- **Instantiate:** `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- **Full test suite (CI-equivalent):** `julia --project=. -e 'using Pkg; Pkg.test()'` — this is the canonical command; CI runs `julia-buildpkg` then `julia-runtest`, which is equivalent. Use generous timeouts: first run triggers precompilation and the performance test `Pkg.add`s three packages, so it can take several minutes with little output. Slow startup is precompilation, not a hang.
- **Focused test file:** `julia --project=. -e 'using SatelliteToolboxSgp4, Test, Dates, DelimitedFiles, Printf, SatelliteToolboxTle; include("test/<file>.jl")'` — extend the `using` list with any extra deps the target file needs (e.g. `ForwardDiff, StaticArrays` for `performance.jl`). There is no `@testitem`/test-name selector; `Pkg.test(test_args=...)` is not read by `runtests.jl`.
- **Examples:** `julia --project=examples examples/<script>.jl` (instantiate first: `julia --project=examples -e 'using Pkg; Pkg.instantiate()'`).
- **Format check (no `--project=.`):** `julia -e 'using JuliaFormatter; format(".")'` then `git diff --exit-code` — `format(".")` returns `true` when nothing changed. JuliaFormatter must be installed in the default env or a shared env like `@format`.

## Code Style

- **Formatter:** JuliaFormatter with `style = "blue"` plus alignment options in `.JuliaFormatter.toml`. This is the source of truth for formatting; run `format(".")` and match its output.
- **CI does not enforce formatting:** no format-check job exists in `.github/workflows/`. Formatting is a convention only; still apply it before committing.
- **No linter is configured** beyond what `performance.jl` runs (Aqua/JET/AllocCheck) inside the test suite.

## Behavioral Constraints

- Preserve the `src/` include order in `src/SatelliteToolboxSgp4.jl`; do not reorder includes or move symbols between files without checking load dependencies.
- `SatelliteToolboxBase` and `SatelliteToolboxTle` are re-exported — changes to their public types (e.g. `Orbit`, TLE structures) propagate into this package's public API. Treat them as part of the surface, not internal.
- The SGP4/SDP4 numerics in `src/sgp4_model.jl` follow Vallado et al., *Revisiting Spacetrack Report #3* (AIAA 2006-6753). Test fixtures in `test/sgp4_tests/aiaa-2006-6753/` are the regression baseline; do not modify expected-result files. If numerics change, expect those fixtures to flag it.
- Match the existing `@testset "..." verbose = true begin ... end` convention (see `test/runtests.jl`, `test/tle.jl`) when adding tests.
- The performance `@testset` is skipped on prerelease builds (`isempty(VERSION.prerelease)` guard) and JET/AllocCheck are additionally skipped on Julia 1.12+. Keep those guards when editing `test/performance.jl`.
- `Manifest.toml` is gitignored — never commit it.

## CI

`.github/workflows/ci.yml` — matrix: Julia `1.10` and `1` (latest stable 1.x), on `ubuntu-latest` / `macos-latest` / `windows-latest`, arch `x64` / `arm64`, with exclusions (no `ubuntu-latest/arm64`, no `macos-latest/x64`, no `windows-latest/arm64`). Steps: `julia-buildpkg` → `julia-runtest` → `julia-processcoverage` → `codecov`. No format-check or docs job.

`.github/workflows/ci-nightly.yml` — Julia `nightly`, same OS/arch matrix, `julia-buildpkg` → `julia-runtest` only (no coverage).

`CompatHelper.yml` and `TagBot.yml` run on schedule for compat bumps and release tagging.

## Not Configured

State these explicitly so they are not invented:
- **No Documenter docs build:** `docs/` exists but contains only `src/assets/logo.png`; there is no `docs/make.jl` or `docs/Project.toml`. Do not add a docs build pipeline without explicit request.
- **No format-check in CI.**
- **No `deps/build.jl`** — `Pkg.build()` is a no-op; `Pkg.test()` alone reproduces CI.
- **No package extensions / weakdeps.**
- **No `test/Project.toml`** — test deps come from `[extras]`/`[targets]`.
- **No `.pre-commit-config.yaml`.**
- **No committed `Manifest.toml`.**
