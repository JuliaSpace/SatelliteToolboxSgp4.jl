<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxTransformations.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

# SatelliteToolboxSgp4.jl

[![CI](https://img.shields.io/github/actions/workflow/status/JuliaSpace/SatelliteToolboxSgp4.jl/ci.yml?style=flat-square&logo=githubactions&logoColor=white&labelColor=475569&label=CI)](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/github/JuliaSpace/SatelliteToolboxSgp4.jl?token=480UYDX6H5&style=flat-square&logo=codecov&logoColor=white&labelColor=475569)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxSgp4.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495D1?style=flat-square&logo=julia&logoColor=white&labelColor=475569)](https://github.com/invenia/BlueStyle)
[![License](https://img.shields.io/github/license/JuliaSpace/SatelliteToolboxSgp4.jl?style=flat-square&logo=readme&logoColor=white&labelColor=475569&color=0284C7)](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/blob/main/LICENSE.txt)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11266935-DB2777?style=flat-square&logo=doi&logoColor=white&labelColor=475569)](https://zenodo.org/doi/10.5281/zenodo.11266935)

This package contains the implementation of the
[SGP4/SDP4](https://en.wikipedia.org/wiki/Simplified_perturbations_models) orbit propagator
for the Julia language.

## Installation

``` julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxSgp4")
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

- `SGP4C_WGS84`: (**DEFAULT**) Constants based on WGS84.
- `SGP4C_WGS72`: Constants based on WGS72.

> **Note**
> The propagator uses the number type of `sgp4c` in all computations. Hence, to propagate
> using another number type, convert the constants first. For example,
> `sgp4_init(tle; sgp4c = Sgp4Constants{Float32}(SGP4C_WGS84))` computes everything
> considering `Float32` numbers.

The SGP4 can also be initialized using an Orbit Mean-Elements Message (OMM) parsed by
[SatelliteToolboxOrbitDataMessages.jl](https://github.com/JuliaSpace/SatelliteToolboxOrbitDataMessages.jl),
which is re-exported by this package, provided that its mean element theory is SGP4:

```julia
julia> omm = read_omm("amazonia_1.xml")

julia> sgp4d = sgp4_init(omm)
```

Finally, the SGP4 can be initialized by passing the mean elements directly. For more
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
in-place, avoiding unnecessary allocations in some cases. The structure can be created with
`Sgp4Propagator(sgp4c)`, which only sets the constants. For more information, see the
function documentation.

### Mean Elements Fitting

This package also provides a way to fit a set of SGP4 mean elements, represented as a `TLE`
or as an `OrbitMeanElementsMessage` (OMM), given a set of osculating state vectors through
the following function:

``` julia
fit_sgp4_mean_elements(::Type{S}, vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) where {S <: Union{TLE, OrbitMeanElementsMessage}, Tjd <: Number, Tv <: AbstractVector} -> S, SMatrix{7, 7, T}, NamedTuple
```

where `S` selects the representation of the output (if it is omitted, the mean elements are
returned as an `OrbitMeanElementsMessage`), and the osculating elements are given
by a set of position vectors `vr_teme` [km] and a set of velocity vectors `vv_teme` [km / s]
represented in the True-Equator, Mean-Equinox reference frame (TEME) at instants in the
array `vjd` [Julian Day].

The algorithm performs a least-square fitting to minimize the residue between the osculating
elements provided by the SGP4 propagator and the input data. It was based on **[4]**.

This function returns the fitted mean elements, the last covariance matrix obtained from
the least-square algorithm, and a `NamedTuple` with the statistics of the fit: `converged`,
`iterations`, `position_rmse` [km], `velocity_rmse` [km / s], and `total_rmse`. When the
output is an OMM, the position and velocity block of the covariance is also stored in the
covariance matrix section of the message, unless the keyword `include_covariance` is
`false`.

> **Note**
> This algorithm version will allocate a new SGP4 propagator with the constants selected by
> the keyword `sgp4c`. If the user wants to reduce the allocations, use the function
> `fit_sgp4_mean_elements!` instead.

The following keywords are available:

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants, whose number type `T` is used in
    the fitting. Only available in `fit_sgp4_mean_elements`, since
    `fit_sgp4_mean_elements!` uses the constants of the propagator.
    (**Default**: `SGP4C_WGS84`)
- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in initial guess (see section **Initial
    Guess**).
    (**Default**: true)
- `include_covariance::Bool`: If `true`, the covariance of the mean position and velocity is
    stored in the output OMM. It is ignored when the output is a TLE.
    (**Default**: true)
- `initial_guess::Union{Nothing, AbstractVector, TLE, OrbitMeanElementsMessage}`: Initial
    guess for the fitting process. If it is `nothing`, the algorithm will obtain an initial
    estimate from the osculating elements in `vr_teme` and `vv_teme`. For more information,
    see the section **Initial Guess**.
    (**Default**: nothing)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix.
    Use `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    `ForwardDiff.jl` automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Union{Number, DateTime}`: Epoch for the fitted mean elements,
    represented by a Julian Day [UTC] or a `DateTime` [UTC].
    (**Default**: vjd[end])
- `template::Union{Nothing, S, NamedTuple}`: Source of the metadata of the output. If it
    is an object of type `S`, its metadata is copied, e.g. the satellite name and number of
    a `TLE` or the header, the metadata, and the TLE-related parameters of an
    `OrbitMeanElementsMessage`. If it is a `NamedTuple`, its entries are passed as keywords
    to the constructor of `S` on top of the default metadata, e.g.
    `(; object_name = "AMAZONIA 1", object_id = "2021-015A", norad_cat_id = 47699)`. If it
    is `nothing`, the metadata is filled with default values.
    (**Default**: nothing)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

#### Initial Guess

This algorithm uses a least-square algorithm to fit a set of mean elements based on a set of
osculating state vectors. Since the system is chaotic, a good initial guess is paramount for
algorithm convergence. We can provide an initial guess using the keyword `initial_guess`.

If `initial_guess` is a `TLE` or an `OrbitMeanElementsMessage`, we update its epoch to the
desired one in `mean_elements_epoch` using the same algorithm as
`update_sgp4_mean_elements_epoch!`. Afterward, we use the updated mean elements as the
initial guess.

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
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted mean elements
are obtained, the algorithm changes their epoch to `mean_elements_epoch`.

> **Note**
> If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the mean
> elements or the state vector. Hence, if `estimate_bstar` is `false`, it will be kept
> constant with this initial value.

#### Examples

```julia
julia> vr_teme = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118],
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736]
       ];

julia> vv_teme = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984],
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928]
       ];

julia> vjd = [2.46002818657856e6, 2.460028190050782e6];

julia> tle, P, stats = fit_sgp4_mean_elements(TLE, vjd, vr_teme, vv_teme; estimate_bstar = false);
ACTION:   Fitting the mean elements.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          5.79375e-09          4.38304e-07          4.38343e-07             -99.9514 %

julia> stats
(converged = true, iterations = 3, position_rmse = 5.793751158786877e-9, velocity_rmse = 4.383038942163017e-7, total_rmse = 4.383421785788155e-7)

julia> tle
TLE:
                      Name : UNDEFINED
          Satellite number : 9999
  International designator : 999999
        Epoch (Year / Day) : 23 /  83.69005079 (2023-03-24T16:33:40.388)
        Element set number : 0
              Eccentricity :   0.00012463
               Inclination :  98.43040000 deg
                      RAAN : 162.11312239 deg
       Argument of perigee : 136.15040637 deg
              Mean anomaly : 241.97934028 deg
           Mean motion (n) :  14.40814157 revs / day
         Revolution number : 0
                        B* :            0 1 / er
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³

julia> omm, P = fit_sgp4_mean_elements(
           vjd,
           vr_teme,
           vv_teme;
           estimate_bstar = false,
           template       = (; object_name = "AMAZONIA 1", norad_cat_id = 47699),
       );
```

### Mean Elements Epoch Update

We can also update the epoch of SGP4 mean elements, represented as a `TLE` or as an
`OrbitMeanElementsMessage`, using the function:

```julia
update_sgp4_mean_elements_epoch(me::S, new_epoch::Union{Number, DateTime}; kwargs...) where {S <: Union{TLE, OrbitMeanElementsMessage}} -> S
```

which returns a new object of the same type obtained by updating the epoch of `me` to
`new_epoch`.

> **Note**
> This algorithm version will allocate a new SGP4 propagator with the constants selected by
> the keyword `sgp4c`. If the user wants to reduce the allocations, use the function
> `update_sgp4_mean_elements_epoch!` instead.

The following keywords are available:

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants, whose number type is used in the
    fitting. Only available in `update_sgp4_mean_elements_epoch`, since
    `update_sgp4_mean_elements_epoch!` uses the constants of the propagator.
    (**Default**: `SGP4C_WGS84`)
- `atol::Number`: Tolerance for the residue absolute value. If, at any iteration, the
    residue is lower than `atol`, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If, at any
    iteration, the relative difference between the residues in two consecutive iterations is
    lower than `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix.
    Use `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    `ForwardDiff.jl` automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: true)

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
                     ṅ / 2 :     -4.4e-07 rev / day²
                     n̈ / 6 :        1e-09 rev / day³

julia> update_sgp4_mean_elements_epoch(tle, DateTime("2023-05-01"))
ACTION:   Updating the epoch of the mean elements.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          2          6.79508e-06          2.29731e-08          6.79512e-06             -99.9999 %

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
                     ṅ / 2 :            0 rev / day²
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
