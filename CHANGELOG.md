SatelliteToolboxSgp4.jl Changelog
=================================

Version 3.0.0
-------------

- ![BREAKING][badge-breaking] Rename the constants `sgp4c_wgs84` and `sgp4c_wgs72` to
  `SGP4C_WGS84` and `SGP4C_WGS72`, following the naming pattern of the SatelliteToolbox
  ecosystem. The `Float32` variants `sgp4c_wgs84_f32` and `sgp4c_wgs72_f32` were removed
  since they can be obtained with the new converting constructor
  `Sgp4Constants{Float32}(SGP4C_WGS84)`.
- ![BREAKING][badge-breaking] The number type of the propagator created by `sgp4_init` and
  `sgp4` is now always the number type of the constants `sgp4c`, as already happened in
  `sgp4_init!`. Previously, non-floating-point inputs (e.g. integers or dual numbers)
  promoted the propagator type, so the same call could return a `Float32` or a `Float64`
  propagator depending on whether an integer literal was used. Other number types must now
  be selected by converting the constants with `Sgp4Constants{T}(sgp4c)`.
- ![BREAKING][badge-breaking] The internal structure `Sgp4DeepSpace` is now immutable and
  stored inline in `Sgp4Propagator`, which gained the resonance integrator state fields
  `atime`, `xli`, and `xni`. The new constructor `Sgp4Propagator{Tepoch}(sgp4c)` returns a
  propagator ready to be initialized by `sgp4_init!`, so the field `sgp4ds` no longer needs
  to be set manually. Measured with BenchmarkTools.jl on an Apple M-series CPU (Julia
  1.12.6), `sgp4!` became 3% to 12% faster for deep space orbits (e.g. 407 ns to 358 ns for
  the 12 h resonant case from AIAA 2006-6753), `copy` became 28% faster (27 ns to 19 ns),
  and `sgp4_init` performs one allocation instead of two, while the initialization time
  and the near-Earth propagation are unchanged.
- ![BREAKING][badge-breaking] Replace `fit_sgp4_tle`, `fit_sgp4_tle!`,
  `update_sgp4_tle_epoch`, and `update_sgp4_tle_epoch!` by `fit_sgp4_mean_elements`,
  `fit_sgp4_mean_elements!`, `update_sgp4_mean_elements_epoch`, and
  `update_sgp4_mean_elements_epoch!`, which select the representation of the mean elements
  through a sink type that can be `TLE` or `OrbitMeanElementsMessage`, defaulting to the
  latter when the sink type is omitted. The metadata of the output is now obtained from
  the new keyword `template`, which can be an object of the output type or a `NamedTuple`
  with keywords of its constructor, instead of the six TLE-specific keywords. When the
  output is an OMM, the position and velocity block of the fit covariance is stored in the
  covariance matrix section of the message, unless the new keyword `include_covariance` is
  `false`. The initial guess can also be an OMM. The fitting now throws the new exception
  `Sgp4FitDivergenceError` instead of an `ErrorException` when the iterations diverge.
- ![BREAKING][badge-breaking] Remove the fields `AE`, `θ²`, and `k₄` from
  `Sgp4Propagator` and the fields `xnddt`, `xndot`, `xldot`, `pe`, `pinc`, `pgh`, `ph`,
  and `pl` from the internal structure `Sgp4DeepSpace`, since they were never read after
  being stored. The fields `sin_M₀`, `cos_M₀`, and `cos_ω₀` were added to
  `Sgp4Propagator` to avoid evaluating trigonometric functions of the initial elements in
  every propagation. Hence, the positional constructor of `Sgp4Propagator` changed.
- ![Feature][badge-feature] The SGP4 propagator can now be initialized using an Orbit
  Mean-Elements Message (OMM) from **SatelliteToolboxOrbitDataMessages.jl** through the
  new methods `sgp4_init(omm)`, `sgp4_init!(sgp4d, omm)`, and `sgp4(Δt, omm)`. The
  package **SatelliteToolboxOrbitDataMessages.jl** is now re-exported.
- ![Feature][badge-feature] Add `show` methods to `Sgp4Propagator`, which print the
  algorithm, the epoch, and the last propagation instant followed by the sections with the
  mean elements, including the semi-major axis recovered from the mean motion and B*, and
  the gravitational constants, instead of every internal field. They use the tree layout of
  **SatelliteToolboxBase.jl** v2.1, which follows the orbit data messages, and overload
  `SatelliteToolboxBase.print_tree_body` so that the wrappers of the propagator can print
  the same body under their own header.
- ![Enhancement][badge-enhancement] The mean elements fitting initializes the propagator
  once per mean state vector and propagates it to all measurements, instead of
  initializing it for every measurement, since the initialization is more expensive than
  the propagation. Measured with BenchmarkTools.jl on an Apple M-series CPU (Julia
  1.12.6), fitting 1001 osculating state vectors of a LEO satellite took 147 ms instead of
  340 ms with `FiniteDiffJacobian()` and 80 ms instead of 159 ms with
  `ForwardDiffJacobian()`. The Jacobians of all measurements are now stored in a buffer
  allocated once per fit (336 KiB for 1001 measurements).
- ![Enhancement][badge-enhancement] Replace **Crayons.jl** with **StyledStrings.jl** to
  decorate the output of the TLE fitting algorithm.
- ![Enhancement][badge-enhancement] The progress line of the TLE fitting algorithm is
  updated in place only when `stdout` supports colors. Otherwise, each iteration prints a
  new line, keeping the output readable when it is redirected to a file.
- ![Enhancement][badge-enhancement] Improve the performance of `sgp4_init!` and `sgp4!` by
  using `cbrt` to compute the semi-major axis and by precomputing the trigonometric
  functions of the initial elements.
- ![Enhancement][badge-enhancement] Simplify the propagation code by removing dead code and
  by adding a converting constructor to `Sgp4Constants`.
- ![Enhancement][badge-enhancement] Allow **SatelliteToolboxBase.jl** v2.
- ![Enhancement][badge-enhancement] **Aqua.jl** and **JET.jl** are now declared as test
  dependencies and run in the new test file `test/quality.jl`, which also verifies with
  JET that the propagation kernel is free of dynamic dispatch.
- ![Bugfix][badge-bugfix] Fix the selection between the SGP4 and SDP4 algorithms, which
  compared the 225 min period threshold against the Kozai mean motion instead of the mean
  motion without the Kozai correction, as in Vallado's implementation.
- ![Bugfix][badge-bugfix] Fix a division by zero in the long-period periodic term for
  orbits with 180° inclination, which returned NaN position and velocity vectors.
- ![Bugfix][badge-bugfix] Fix the deep space initialization, which did not drop the
  lunar-solar node term for retrograde orbits with inclination above 177°, as in Vallado's
  implementation.
- ![Bugfix][badge-bugfix] Fix the clamp that keeps `sin(i₀)` away from zero in the deep
  space initialization, which had no effect when the inclination was exactly 0 or π.
- ![Bugfix][badge-bugfix] Fix `fit_sgp4_tle!`, which threw `UndefVarError` when
  `max_iterations` was lower than 1. An `ArgumentError` is now thrown.
- ![Bugfix][badge-bugfix] Fix wrong and stale comments in the SGP4 model.

Version 2.5.0
-------------

- ![Enhancement][badge-enhancement] Remove the unused fields `xnq`, `omegaq`, `omgdt`,
  `ilsz`, `pgh0`, `ph0`, `pe0`, `pinc0`, and `pl0` from the internal structure
  `Sgp4DeepSpace`.
- ![Enhancement][badge-enhancement] Reduce allocations in `fit_sgp4_tle!` by converting the
  measurements to static vectors before the least-square iterations.
- ![Bugfix][badge-bugfix] Fix the TLE epoch year computed in `fit_sgp4_tle!`, which was
  negative for epochs between 1980 and 1999.
- ![Bugfix][badge-bugfix] Fix the correction limiting algorithm in `fit_sgp4_tle!`, which
  froze any state component that was exactly zero.
- ![Bugfix][badge-bugfix] Fix the Kepler solver exit tolerance in `sgp4!`, which was
  unreachable for types with lower precision than `Float64`, causing the solver to always
  run all iterations when using `Float32`.
- ![Bugfix][badge-bugfix] Fix an access to undefined fields in the deep space
  initialization, which threw `UndefRefError` when using non-isbits number types.

Version 2.4.2
-------------

- ![Bugfix][badge-bugfix] Fix type-stability regressions in `sgp4_init!`, `sgp4!`, and
  `_dsper!` so the SGP4/SDP4 propagation hot path stays type-stable.
- ![Enhancement][badge-enhancement] Reduce allocations in `fit_sgp4_tle!` by replacing the
  diagonal weight matrix with a weight vector and updating the per-iteration accumulators.

Version 2.4.1
-------------

- ![Enhancement][badge-enhancement] Use the types related to the Jacobian definition from
  **SatelliteToolboxBase.jl**.

Version 2.4.0
-------------

- ![Feature][badge-feature] The SGP4 propagator now supports differentiability. Hence, the
  user can now also uses **ForwardDiff.jl** to compute the Jacobian in `fit_sgp4_tle`.

Version 2.3.0
-------------

- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).

Version 2.2.0
-------------

- ![Feature][badge-feature] We implemented `copy` for all structures.

Version 2.1.4
-------------

- ![Enhancement][badge-enhancement] We reduced the allocations in function `fit_sgp4_tle!`.

Version 2.1.3
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.

Version 2.1.2
-------------

- ![Enhancement][badge-enhancement] We updated the dependency compatibility bounds.

Version 2.1.1
-------------

- ![Enhancement][badge-enhancement] **SnoopPrecompile.jl** was replaced by
  **PrecompileTools.jl**.

Version 2.1.0
-------------

- ![Feature][badge-feature] We added the functions `update_sgp4_tle_epoch` and
  `update_sgp4_tle_epoch!` to update the epoch of a SGP4 TLE.
- ![Feature][badge-feature] We added the functions `fit_sgp4_tle` and `fit_sgp4_tle!` to fit
  a SGP4 TLE using a set of osculating state vectors represented in the TEME reference
  frame.

Version 2.0.0
-------------

- ![BREAKING][badge-breaking] We removed the field `Ω1` from the structure `Sgp4Propagator`
  because it was not being used in the propagation.
- ![BREAKING][badge-breaking] The structure `Sgp4Propagator` is not a `Base.@kwdef` anymore.
  We also added custom constructors to help initialize an instance with uninitiated fields.
  Hence, if one creates the SGP4 structure directly, i.e., without using `sgp4_init`, this
  version is breaking.
- ![Feature][badge-feature] We added the function `sgp4_init!` to initialize a SGP4
  propagator in-place, avoiding unnecessary allocations.
- ![Enhancement][badge-enhancement] The code was slightly improved, leading to a 5% speed
  gain in initialization and 4% speed gain in propagation.

Version 1.0.1
-------------

- ![Enhancement][badge-enhancement] Bump the version of **SatelliteToolboxBase.jl**.

Version 1.0.0
-------------

- ![Enhancement][badge-enhancement] After all testing and source-code cleaning, we can mark
  this package as stable, reaching v1.
- ![Enhancement][badge-enhancement] We use
  [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl) as
  dependencies to provide some functionalities for the algorithm.

Version 0.1.1
-------------

- ![Enhancement][badge-enhancement] We added precompilation statements to improve
  performance.
- ![Enhancement][badge-enhancement] The code was refactored to follow BlueStyle, and
  line-width was increase to 92, leading to a better source-code organization.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the submodule in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square
