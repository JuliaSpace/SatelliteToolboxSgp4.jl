SatelliteToolboxSgp4.jl
=======================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxSgp4.jl/branch/main/graph/badge.svg?token=480UYDX6H5)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxSgp4.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package contains the implementation of the
[SGP4/SDP4](https://en.wikipedia.org/wiki/Simplified_perturbations_models) orbit propagator
for the Julia language.

## Installation

``` julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxSgp4")
```

## Usage

### Orbit Propagation

First, we need to initialize the structure that contains the information to propagate the
orbit using the function `sgp4_init`. Usually, we pass a
[TLE](https://github.com/JuliaSpace/SatelliteToolboxTle.jl) to initialize the SGP4
algorithm:

```julia
julia> using SatelliteToolboxTle

julia> tle = tle"""
       AMAZONIA 1
       1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
       2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652
       """

julia> sgp4d = sgp4_init(tle)
```

`sgp4_init` supports the keyword `sgp4c` to select the constants used to propagate the
orbit. It must be an object of type `Sgp4Constants`. The following constants are already
defined in this package:

- `sgp4c_wgs84`: (**DEFAULT**) Constants based on WGS84 using `Float64`.
- `sgp4c_wgs84_f32`: Constants based on WGS84 using `Float32`.
- `sgp4c_wgs72`: Constants based on WGS72 using `Float64`.
- `sgp4c_wgs72_f32`: Constants based on WGS72 using `Float32`.

> **Note**
> The propagator will use the same type of object `sgp4c` to propagate the orbit. Hence, if
> one selects `sgp4c_wgs84_f32`, the SGP4 will compute everything considering `Float32`
> numbers.

The SGP4 can also be initialized by passing the mean elements directly. For more
information, see the documentation of the function `sgp4_init`.

Afterward, we can propagate the orbit using the function `sgp4!(sgp4d, t)` that propagates
the mean elements defined in `sgp4d` by `t` minutes. This function returns the position [km]
and velocity [km/s] vectors represented in the True Equator, Mean Equinox (TEME) reference
frame.

```julia
# Propagate the orbit for 10 minutes.
julia> r_teme, v_teme = sgp4!(sgp4d, 10)
([-5300.1473032595195, 2356.780136349037, 4149.0611906521035], [4.464838382952148, -0.5106103512199875, 5.9760603775620815])
```

> **Warning**
> We do not use SI units here to keep consistency with the original SGP4/SDP4 algorithms.

The function `sgp4(t, args...)` creates the propagator and propagates the orbit defined in
`args...` by `t` minutes. It returns the same information as the function `sgp4!` and the
initialized propagator structure. `args...` must be the same arguments supported by `sgp4!`.

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

We also have the function `sgp4_init!` that initializes a SGP4 propagator structure
in-place, avoiding unnecessary allocations in some cases. For more information, see the
function documentation.

### TLE Fitting

This package also provides a way to fit a TLE for the SGP4 algorithm given a set of
osculating state vectors through the following function:

``` julia
function fit_sgp4_tle(vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector}
```

where the osculating elements are given by a set of position vectors `vr_teme` [km] and a
set of velocity vectors `vv_teme` [km / s] represented in the True-Equator, Mean-Equinox
reference frame (TEME) at instants in the array `vjd` [Julian Day].

The algorithm performs a least-square fitting to minimize the residue between the osculating
elements provided by the SGP4 propagator and the input data. It was based on **[4]**.

This function returns the fitted TLE and the last covariance matrix obtained from the
least-square algorithm.

> **Note**
> This algorithm version will allocate a new SGP4 propagator with the default constants
> `sgp4c_wgs84`. If another set of constants are required or if the user wants to reduce the
> allocations, use the function `fit_sgp4_tle!` instead.

The following keywords are avaible:

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops. (**Default** = 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in initial guess (see section  **Initial
    Guess**). (**Default** = true)
- `initial_guess::Union{Nothing, AbstractVector, TLE}`: Initial guess for the TLE fitting
    process. If it is `nothing`, the algorithm will obtain an initial estimate from the
    osculating elements in `vr_teme` and `vv_teme`. For more information, see the section
    **Initial Guess**. (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`. (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted TLE. (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal. (**Default** = `@SVector(ones(Bool, 6))`)
- `classification::Char`: Satellite classification character for the output TLE.
    (**Default** = 'U')
- `element_set_number::Int`: Element set number for the output TLE. (**Default** = 0)
- `international_designator::String`: International designator string for the output TLE.
    (**Default** = "999999")
- `name::String`: Satellite name for the output TLE. (**Default** = "UNDEFINED")
- `revolution_number::Int`: Revolution number for the output TLE. (**Default** = 0)
- `satellite_number::Int`: Satellite number for the output TLE. (**Default** = 9999)

#### Initial Guess

This algorithm uses a least-square algorithm to fit a TLE based on a set of osculating state
vectors. Since the system is chaotic, a good initial guess is paramount for algorithm
convergence. We can provide an initial guess using the keyword `initial_guess`.

If `initial_guess` is a `TLE`, we update the TLE epoch using the function
`update_sgp4_tle_epoch!` to the desired one in `mean_elements_epoch`. Afterward, we use this
new TLE as the initial guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

``` julia
┌                                    ┐
│ IDs 1 to 3: Mean position [km]     │
│ IDs 4 to 6: Mean velocity [km / s] │
│ ID  7:      Bstar         [1 / er] │
└                                    ┘
```

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted TLE is
obtained, the algorithm uses the function [`update_sgp4_tle_epoch!`](@ref) to change its
epoch to `mean_elements_epoch`.

> **Note**
> If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the TLE or
> the state vector. Hence, if `estimate_bstar` is `false`, it will be kept constant with
> this initial value.

#### Example

```julia
julia> vjd  = [2.46002818657856e6];

julia> vr_teme = [@SVector [-6792.402703741442, 2192.6458461287293, 0.18851758695295118]];

julia> vv_teme = [@SVector [0.3445760107690598, 1.0395135806993514, 7.393686131436984]];

julia> tle, P = fit_sgp4_tle(vjd, vr_teme, vv_teme; estimate_bstar = false)
ACTION:   Fitting the TLE.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:         47           2.7699e-07           2.5122e-10           2.7699e-07                 -100 %

(TLE: UNDEFINED (Epoch = 2023-03-24T16:28:40.388), [2.5682204663112826 0.5152694462560151 … -1.2385529801114805 0.0; 0.5152694462560317 4.6021551786709445 … -0.1449556257318217 0.0; … ; -1.2385529801114785 -0.14495562573181023 … 0.999604694355324 0.0; 0.0 0.0 … 0.0 0.0])

julia> tle
TLE:
                      Name : UNDEFINED
          Satellite number : 9999
  International designator : 999999
        Epoch (Year / Day) : 23 /  83.68657856 (2023-03-24T16:28:40.388)
        Element set number : 0
              Eccentricity :   0.00012470
               Inclination :  98.43040000 deg
                      RAAN : 162.10970000 deg
       Argument of perigee : 136.20170000 deg
              Mean anomaly : 223.92830000 deg
           Mean motion (n) :  14.40814394 revs / day
         Revolution number : 0
                        B* :            0 1 / er
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```

### TLE Epoch Update

We can also update SGP4 TLE epoch using the function:

```julia
function update_sgp4_tle_epoch(tle::TLE, new_epoch::Union{Number, DateTime}; kwargs...)
```

which returns a new TLE obtained by updating the epoch of `tle` to `new_epoch`.

> **Note**
> This algorithm version will allocate a new SGP4 propagator with the default constants
> `sgp4c_wgs84`. If another set of constants are required or if the user wants to reduce the
> allocations, use the function `update_sgp4_tle_epoch!` instead.

The following keywords are avaible:

- `atol::Number`: Tolerance for the residue absolute value. If, at any iteration, the
    residue is lower than `atol`, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If, at any
    iteration, the relative difference between the residues in two consecutive iterations is
    lower than `rtol`, the computation loop stops. (**Default** = 2e-4)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
    
#### Examples

``` julia
julia> tle = tle"""
           AMAZONIA 1
           1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
           2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""
TLE:
                      Name : AMAZONIA 1
          Satellite number : 47699
  International designator : 21015A
        Epoch (Year / Day) : 23 /  83.68657856 (2023-03-24T16:28:40.388)
        Element set number : 999
              Eccentricity :   0.00012470
               Inclination :  98.43040000 deg
                      RAAN : 162.10970000 deg
       Argument of perigee : 136.20170000 deg
              Mean anomaly : 223.92830000 deg
           Mean motion (n) :  14.40814394 revs / day
         Revolution number : 10865
                        B* :      4.3e-05 1 / er
                     ṅ / 2 :     -4.4e-07 rev / day²
                     n̈ / 6 :        1e-09 rev / day³

julia> update_sgp4_tle_epoch(tle, DateTime("2023-05-01"))
ACTION:   Fitting the TLE.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          2          6.78579e-06          2.29683e-08          6.78583e-06             -99.9999 %

TLE:
                      Name : AMAZONIA 1
          Satellite number : 47699
  International designator : 21015A
        Epoch (Year / Day) : 23 / 121.00000000 (2023-05-01T00:00:00)
        Element set number : 999
              Eccentricity :   0.00012481
               Inclination :  98.43040000 deg
                      RAAN : 198.88793445 deg
       Argument of perigee :  24.20140448 deg
              Mean anomaly :  86.69370896 deg
           Mean motion (n) :  14.40824649 revs / day
         Revolution number : 10865
                        B* :      4.3e-05 1 / er
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```

## References

The code in this package was built using the following references:

- **[1]** Hoots, F. R., Roehrich, R. L (1980). *Models for Propagation of NORAD Elements
  Set*. **Spacetrack Report No. 3**.
- **[2]** Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). *Revisiting
  Spacetrack Report #3: Rev1*. **AIAA**.
- **[3]** SGP4 Source code of [STRF](https://github.com/cbassa/strf), which the C code was
  converted by Paul. S. Crawford and Andrew R. Brooks.
- **[4]** Vallado, D. A., Crawford, P (2008). *SGP4 Orbit Determination*. **American Institute
  of Aeronautics ans Astronautics**.
