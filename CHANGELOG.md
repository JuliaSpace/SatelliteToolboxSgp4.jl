SatelliteToolboxSgp4.jl Changelog
=================================

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

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg
