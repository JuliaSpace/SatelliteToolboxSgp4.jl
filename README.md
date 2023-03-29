SatelliteToolboxSgp4.jl
=======================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxSgp4.jl/branch/main/graph/badge.svg?token=480UYDX6H5)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxSgp4.jl)

This package contains the implementation of the
[SGP4/SDP4](https://en.wikipedia.org/wiki/Simplified_perturbations_models) orbit
propagator for the Julia language.

## Installation

``` julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxSgp4")
```

## Usage

First, we need to initialize the structure that contains the information to
propagate the orbit using the function `sgp4_init`. Usually, we pass a
[TLE](https://github.com/JuliaSpace/SatelliteToolboxTle.jl) to initialize the
SGP4 algorithm:

```julia
julia> using SatelliteToolboxTle

julia> tle = tle"""
       AMAZONIA 1
       1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
       2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652
       """

julia> sgp4d = sgp4_init(tle)
```

`sgp4_init` supports the keyword `sgp4c` to select the constants used to
propagate the orbit. It must be an object of type `Sgp4Constants`. The following
constants are already defined in this package:

- `sgp4c_wgs84`: (**DEFAULT**) Constants based on WGS84 using `Float64`.
- `sgp4c_wgs84_f32`: Constants based on WGS84 using `Float32`.
- `sgp4c_wgs72`: Constants based on WGS72 using `Float64`.
- `sgp4c_wgs72_f32`: Constants based on WGS72 using `Float32`.

> **Note**
> The propagator will use the same type of object `sgp4c` to propagate the
> orbit. Hence, if one selects `sgp4c_wgs84_f32`, the SGP4 will compute
> everything considering `Float32` numbers.

The SGP4 can also be initialized by passing the mean elements directly. For more
information, see the documentation of the function `sgp4_init`.

Afterward, we can propagate the orbit using the function `sgp4!(sgp4d, t)` that
propagates the mean elements defined in `sgp4d` by `t` minutes. This function
returns the position [km] and velocity [km/s] vectors represented in the
True Equator, Mean Equinox (TEME) reference frame.

```julia
# Propagate the orbit for 10 minutes.
julia> r_teme, v_teme = sgp4!(sgp4d, 10)
([-5300.1473032595195, 2356.780136349037, 4149.0611906521035], [4.464838382952148, -0.5106103512199875, 5.9760603775620815])
```

> **Warning**
> We do not use SI units here to keep consistency with the original SGP4/SDP4
> algorithms.

The function `sgp4(t, args...)` creates the propagator and propagates the orbit
defined in `args...` by `t` minutes. It returns the same information as the
function `sgp4!` and the initialized propagator structure. `args...` must be the
same arguments supported by `sgp4!`.

``` julia
julia> r_teme, v_teme, sgp4d = sgp4(10, tle)

julia> r_teme
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -5300.1473032595195
  2356.780136349037
  4149.0611906521035

julia> v_teme
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  4.464838382952148
 -0.5106103512199875
  5.9760603775620815
```

`sgp4` also supports the same keywords arguments as `sgp4!`.

## References

The code in this package was built using the following references:

- **[1]** Hoots, F. R., Roehrich, R. L (1980). *Models for Propagation of NORAD
  Elements Set*. **Spacetrack Report No. 3**.
- **[2]** Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006).
  *Revisiting Spacetrack Report #3: Rev1*. **AIAA**.
- **[3]** SGP4 Source code of [STRF](https://github.com/cbassa/strf), which the
  C code was converted by Paul. S. Crawford and Andrew R. Brooks.
