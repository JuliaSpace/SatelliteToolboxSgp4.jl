## Description #############################################################################
#
# Least-square algorithm to fit SGP4 mean elements, represented as a TLE or an Orbit
# Mean-Elements Message (OMM), from a set of osculating state vectors.
#
## References ##############################################################################
#
# [1] Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute of
#     Aeronautics and Astronautics.
#
############################################################################################

export fit_sgp4_mean_elements, fit_sgp4_mean_elements!
export update_sgp4_mean_elements_epoch, update_sgp4_mean_elements_epoch!

# Types that can represent a set of SGP4 mean elements.
const _SGP4_MEAN_ELEMENTS = Union{TLE, OrbitMeanElementsMessage}

# Types accepted as the initial guess of the fitting algorithm.
const _INITIAL_GUESS_T = Union{Nothing, AbstractVector, TLE, OrbitMeanElementsMessage}

"""
    fit_sgp4_mean_elements(
        ::Type{S},
        vjd::AbstractVector{Tjd},
        vr_teme::AbstractVector{Tv},
        vv_teme::AbstractVector{Tv};
        kwargs...,
    ) where {
        S <: Union{TLE, OrbitMeanElementsMessage},
        Tjd <: Number,
        Tv <: AbstractVector
    } -> S, SMatrix{7, 7, T}, NamedTuple

    fit_sgp4_mean_elements(
        vjd::AbstractVector{Tjd},
        vr_teme::AbstractVector{Tv},
        vv_teme::AbstractVector{Tv};
        kwargs...,
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> OrbitMeanElementsMessage, SMatrix{7, 7, T}, NamedTuple

Fit a set of SGP4 mean elements represented as an object of type `S`, which can be `TLE` or
`OrbitMeanElementsMessage`, using the osculating elements represented by a set of position
vectors `vr_teme` [km] and a set of velocity vectors `vv_teme` [km / s] represented in the
True-Equator, Mean-Equinox reference frame (TEME) at instants in the array `vjd` [Julian
Day, UTC]. If `S` is omitted, the mean elements are returned as an
`OrbitMeanElementsMessage`.

This algorithm was based on **[1]**. It can fail if the least-square iterations diverge.

!!! note

    This algorithm version will allocate a new SGP4 propagator with the constants `sgp4c`.
    If the allocation must be avoided, use the function [`fit_sgp4_mean_elements!`](@ref)
    instead.

See also: [`fit_sgp4_mean_elements!`](@ref), [`update_sgp4_mean_elements_epoch`](@ref).

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)),
    whose number type `T` is used in the fitting.
    (**Default**: `SGP4C_WGS84`)

The other keywords are the same as in [`fit_sgp4_mean_elements!`](@ref).

# Returns

- `S`: The fitted mean elements.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm, whose state
    is the mean position [km], the mean velocity [km / s], and the drag term B* [1 / er].
- `NamedTuple`: Statistics of the least-square algorithm (see
    [`fit_sgp4_mean_elements!`](@ref)).

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics and Astronautics.

# Extended help

## Throws

- `ArgumentError`: If the lengths of `vjd`, `vr_teme`, and `vv_teme` differ, if the weight
    vector or the initial guess vector have the wrong size, if `max_iterations` is lower
    than 1, or if a `NamedTuple` template contains a field set by the fit.
- `Sgp4FitDivergenceError`: If the least-square iterations diverge.

## Examples

```julia-repl
julia> vr_teme = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118],
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736]
       ];

julia> vv_teme = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984],
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928]
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> tle, P, stats = fit_sgp4_mean_elements(TLE, vjd, vr_teme, vv_teme; estimate_bstar = false);
ACTION:   Fitting the mean elements.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          5.79375e-09          4.38304e-07          4.38343e-07             -99.9514 %

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

julia> stats
(converged = true, iterations = 3, position_rmse = 5.793751158786877e-9, velocity_rmse = 4.383038942163017e-7, total_rmse = 4.383421785788155e-7)

julia> omm, P, stats = fit_sgp4_mean_elements(
           OrbitMeanElementsMessage, vjd, vr_teme, vv_teme; estimate_bstar = false
       );
```
"""
function fit_sgp4_mean_elements(
    ::Type{S},
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    sgp4c::Sgp4Constants = SGP4C_WGS84,
    kwargs...,
) where {S <: _SGP4_MEAN_ELEMENTS, Tjd <: Number, Tv <: AbstractVector}
    sgp4d = Sgp4Propagator(sgp4c)
    return fit_sgp4_mean_elements!(sgp4d, S, vjd, vr_teme, vv_teme; kwargs...)
end

function fit_sgp4_mean_elements(
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_sgp4_mean_elements(
        OrbitMeanElementsMessage, vjd, vr_teme, vv_teme; kwargs...
    )
end

"""
    fit_sgp4_mean_elements!(
        sgp4d::Sgp4Propagator{Tepoch, T},
        ::Type{S},
        vjd::AbstractVector{Tjd},
        vr_teme::AbstractVector{Tv},
        vv_teme::AbstractVector{Tv};
        kwargs...,
    ) where {
        Tepoch <: Number,
        T <: Number,
        S <: Union{TLE, OrbitMeanElementsMessage},
        Tjd <: Number,
        Tv <: AbstractVector,
    } -> S, SMatrix{7, 7, T}, NamedTuple

    fit_sgp4_mean_elements!(
        sgp4d::Sgp4Propagator{Tepoch, T},
        vjd::AbstractVector{Tjd},
        vr_teme::AbstractVector{Tv},
        vv_teme::AbstractVector{Tv};
        kwargs...,
    ) where {
        Tepoch <: Number,
        T <: Number,
        Tjd <: Number,
        Tv <: AbstractVector,
    } -> OrbitMeanElementsMessage, SMatrix{7, 7, T}, NamedTuple

Fit a set of SGP4 mean elements for the propagator `sgp4d`, represented as an object of
type `S`, which can be `TLE` or `OrbitMeanElementsMessage`, using the osculating elements
represented by a set of position vectors `vr_teme` [km] and a set of velocity vectors
`vv_teme` [km / s] represented in the True-Equator, Mean-Equinox reference frame (TEME) at
instants in the array `vjd` [Julian Day, UTC]. If `S` is omitted, the mean elements are
returned as an `OrbitMeanElementsMessage`.

This algorithm was based on **[1]**. It can fail if the least-square iterations diverge.

!!! note

    The SGP4 orbit propagator `sgp4d` will be initialized with the mean elements returned
    by the function.

See also: [`fit_sgp4_mean_elements`](@ref), [`update_sgp4_mean_elements_epoch!`](@ref).

# Keywords

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
    (**Default**: `true`)
- `include_covariance::Bool`: If `true`, the covariance of the mean position and velocity
    obtained by the least-square algorithm is stored in the covariance matrix section of
    the output, represented in the TEME reference frame. It is only used when `S` is
    `OrbitMeanElementsMessage`.
    (**Default**: `true`)
- `initial_guess::Union{Nothing, AbstractVector, TLE, OrbitMeanElementsMessage}`: Initial
    guess for the fitting process. If it is `nothing`, the algorithm will obtain an initial
    estimate from the osculating elements in `vr_teme` and `vv_teme`. For more information,
    see the section **Initial Guess**.
    (**Default**: `nothing`)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix.
    Use `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    **ForwardDiff.jl** automatic differentiation.
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
- `mean_elements_epoch::Union{Number, DateTime}`: Epoch of the fitted mean elements,
    represented by a Julian Day [UTC] or a `DateTime` [UTC].
    (**Default**: `vjd[end]`)
- `template::Union{Nothing, S, NamedTuple}`: Source of the metadata of the output. If it
    is an object of type `S`, its metadata is copied, e.g. the satellite name and number of
    a `TLE` or the header, the metadata, and the TLE-related parameters of an
    `OrbitMeanElementsMessage`. If it is a `NamedTuple`, its entries are passed as keywords
    to the constructor of `S` on top of the default metadata, e.g.
    `(; name = "AMAZONIA 1", satellite_number = 47699)` for a `TLE` or
    `(; object_name = "AMAZONIA 1", object_id = "2021-015A", norad_cat_id = 47699)` for
    an `OrbitMeanElementsMessage`. In both cases, only the mean elements, the epoch, the
    drag term, and the covariance matrix are set by the fitted values, and a `NamedTuple`
    must not contain them. If it is `nothing`, the metadata is filled with default values.
    (**Default**: `nothing`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `SVector{6, Bool}(true, true, true, true, true, true)`)

# Returns

- `S`: The fitted mean elements.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm, whose state
    is the mean position [km], the mean velocity [km / s], and the drag term B* [1 / er].
- `NamedTuple`: Statistics of the least-square algorithm with the following fields:
    - `converged::Bool`: `true` if the iterations stopped because the residue was lower
        than `atol` or its relative variation was lower than `rtol`, or `false` if they
        stopped by reaching `max_iterations`.
    - `iterations::Int`: Number of iterations performed.
    - `position_rmse::T`: RMSE of the position residue in the last iteration [km].
    - `velocity_rmse::T`: RMSE of the velocity residue in the last iteration [km / s].
    - `total_rmse::T`: Weighted RMSE of the residue in the last iteration.

    The statistics refer to the fitting of the mean elements. If their epoch is updated
    afterward to match `mean_elements_epoch`, the statistics of that update are not
    returned.

# Initial Guess

This algorithm uses a least-square algorithm to fit a set of mean elements based on a set
of osculating state vectors. Since the system is chaotic, a good initial guess is paramount
for algorithm convergence. We can provide an initial guess using the keyword
`initial_guess`.

If `initial_guess` is a `TLE` or an `OrbitMeanElementsMessage`, we update its epoch to the
desired one in `mean_elements_epoch` using the same algorithm as
[`update_sgp4_mean_elements_epoch!`](@ref). Afterward, we use the updated mean elements as
the initial guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted mean
elements are obtained, the algorithm changes their epoch to `mean_elements_epoch` using the
same algorithm as [`update_sgp4_mean_elements_epoch!`](@ref).

!!! note

    If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the mean
    elements or the state vector. Hence, if `estimate_bstar` is `false`, it will be kept
    constant with this initial value.

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics and Astronautics.

# Extended help

## Throws

- `ArgumentError`: If the lengths of `vjd`, `vr_teme`, and `vv_teme` differ, if the weight
    vector or the initial guess vector have the wrong size, if `max_iterations` is lower
    than 1, or if a `NamedTuple` template contains a field set by the fit.
- `Sgp4FitDivergenceError`: If the least-square iterations diverge.

## Examples

```julia-repl
# Allocate a new SGP4 orbit propagator using a dummy TLE.
julia> sgp4d = sgp4_init(tle\"\"\"
           ISS (ZARYA)
           1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
           2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537\"\"\");

julia> vr_teme = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118],
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736]
       ];

julia> vv_teme = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984],
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928]
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> tle, P, stats = fit_sgp4_mean_elements!(
           sgp4d, TLE, vjd, vr_teme, vv_teme; estimate_bstar = false
       );
ACTION:   Fitting the mean elements.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          5.79375e-09          4.38304e-07          4.38343e-07             -99.9514 %

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
```
"""
function fit_sgp4_mean_elements!(
    sgp4d::Sgp4Propagator{Tepoch, T},
    ::Type{S},
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    atol::Number                                 = 2e-4,
    rtol::Number                                 = 2e-4,
    estimate_bstar::Bool                         = true,
    include_covariance::Bool                     = true,
    initial_guess::_INITIAL_GUESS_T              = nothing,
    jacobian_method::AbstractJacobianMethod      = FiniteDiffJacobian(),
    jacobian_perturbation::Number                = 1e-3,
    jacobian_perturbation_tol::Number            = 1e-7,
    max_iterations::Int                          = 50,
    mean_elements_epoch::Union{Number, DateTime} = vjd[end],
    template::Union{Nothing, S, NamedTuple}      = nothing,
    verbose::Bool                                = true,
    weight_vector::AbstractVector                = SVector{6, Bool}(true, true, true, true, true, true),
) where {
    Tepoch <: Number,
    T <: Number,
    S <: _SGP4_MEAN_ELEMENTS,
    Tjd <: Number,
    Tv <: AbstractVector,
}
    # Unpack.
    sgp4c = sgp4d.sgp4c

    # Desired epoch of the mean elements [Julian Day].
    desired_epoch = _julian_day(mean_elements_epoch)

    # Number of available measurements.
    num_measurements = length(vjd)

    # Check the inputs.
    length(vr_teme) != num_measurements && throw(
        ArgumentError("The number of elements in `vjd` and `vr_teme` must be the same.")
    )

    length(vv_teme) != num_measurements && throw(
        ArgumentError("The number of elements in `vjd` and `vv_teme` must be the same.")
    )

    if length(weight_vector) != 6
        throw(ArgumentError("The weight vector must have 6 elements."))
    end

    if (initial_guess isa AbstractVector) && (length(initial_guess) != 7)
        throw(ArgumentError("The initial guess state vector must have 7 elements."))
    end

    if max_iterations < 1
        throw(ArgumentError("The maximum number of iterations must be at least 1."))
    end

    # Check if `stdout` supports colors. This flag also selects whether the progress line
    # is updated in place using terminal escape sequences.
    has_color = get(stdout, :color, false)::Bool

    # Assemble the weight vector (diagonal of the weight matrix).
    W = @SVector T[
        weight_vector[1],
        weight_vector[2],
        weight_vector[3],
        weight_vector[4],
        weight_vector[5],
        weight_vector[6],
    ]

    # Convert the measurements to static vectors to avoid allocations inside the fitting
    # loop when the user provides dynamic vectors.
    vy = [
        SVector{6, T}(
            vr_teme[k][1],
            vr_teme[k][2],
            vr_teme[k][3],
            vv_teme[k][1],
            vv_teme[k][2],
            vv_teme[k][3],
        ) for k in 1:num_measurements
    ]

    # Keywords shared with the epoch update algorithm.
    update_kwargs = (;
        atol,
        rtol,
        jacobian_method,
        jacobian_perturbation,
        jacobian_perturbation_tol,
        max_iterations,
        verbose,
        has_color,
    )

    # == Initial Guess of the Mean Elements ================================================

    if initial_guess isa _SGP4_MEAN_ELEMENTS
        epoch = desired_epoch

        verbose && _fit_print_action(
            has_color,
            "Updating the epoch of the initial guess to match the desired one.",
        )

        # Convert the mean elements to the mean state vector at their own epoch and update
        # it to the desired epoch.
        x₁ = _update_sgp4_mean_state_vector!(
            sgp4d,
            _mean_state_vector(initial_guess; sgp4c = sgp4c),
            _mean_elements_epoch(initial_guess),
            epoch;
            update_kwargs...,
        )

    elseif initial_guess isa AbstractVector
        # In this case, the user must ensure that the provided mean elements are related to
        # the selected `mean_elements_epoch`.
        epoch = desired_epoch
        x₁    = SVector{7, T}(initial_guess...)

    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        id = _closest_measurement(vjd, desired_epoch)

        epoch = vjd[id]
        x₁    = SVector{7, T}(vr_teme[id]..., vv_teme[id]..., estimate_bstar ? T(0.00001) : T(0))
    end

    # == Least-Square Fitting ==============================================================

    verbose && _fit_print_action(has_color, "Fitting the mean elements.")

    x₂, P, stats = _fit_sgp4_mean_state_vector!(
        sgp4d,
        vjd,
        vy,
        x₁,
        epoch,
        W;
        atol,
        rtol,
        estimate_bstar,
        jacobian_method,
        jacobian_perturbation,
        jacobian_perturbation_tol,
        max_iterations,
        verbose,
        has_color,
    )

    # == Epoch Update ======================================================================

    # Update the epoch of the fitted mean elements to match the desired one.
    if abs(epoch - desired_epoch) > 0.001 / 86400
        verbose && _fit_print_action(
            has_color,
            "Updating the epoch of the fitted mean elements to match the desired one.",
        )

        x₂    = _update_sgp4_mean_state_vector!(sgp4d, x₂, epoch, desired_epoch; update_kwargs...)
        epoch = desired_epoch
    end

    # == Output ============================================================================

    # Build the mean elements with the requested representation.
    me = _build_mean_elements(
        S,
        x₂,
        epoch;
        sgp4c      = sgp4c,
        template   = template,
        covariance = include_covariance ? SMatrix{6, 6, T}(P[1:6, 1:6]) : nothing,
    )

    # Initialize the propagator with the fitted mean elements.
    sgp4_init!(sgp4d, me)

    return me, P, stats
end

function fit_sgp4_mean_elements!(
    sgp4d::Sgp4Propagator,
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_sgp4_mean_elements!(
        sgp4d, OrbitMeanElementsMessage, vjd, vr_teme, vv_teme; kwargs...
    )
end

"""
    update_sgp4_mean_elements_epoch(
        me::S,
        new_epoch::Union{Number, DateTime};
        kwargs...,
    ) where {S <: Union{TLE, OrbitMeanElementsMessage}} -> S

Update the epoch of the SGP4 mean elements `me`, which can be a `TLE` or an
`OrbitMeanElementsMessage`, to `new_epoch`, represented by a Julian Day [UTC] or a
`DateTime` [UTC].

!!! note

    This algorithm version will allocate a new SGP4 propagator with the constants `sgp4c`.
    If the allocation must be avoided, use the function
    [`update_sgp4_mean_elements_epoch!`](@ref) instead.

This function uses the following algorithm to update the epoch:

1. Initialize the SGP4 propagator with `me`;
2. Propagate the orbit to `new_epoch` and obtain the osculating state vector in TEME
    reference frame; and
3. Fit a new set of mean elements that provides the same osculating state vector but
    considering the new epoch.

The third step uses the same least-square algorithm as [`fit_sgp4_mean_elements!`](@ref).
Hence, some keywords are related to it and the function can fail if the iterations diverge.
The metadata of `me` is kept in the output, except for the creation date of an
`OrbitMeanElementsMessage`, which is set to the current time, and its covariance matrix,
which is removed.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)),
    whose number type is used in the fitting.
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
    **ForwardDiff.jl** automatic differentiation.
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
    (**Default**: `true`)

# Extended help

## Throws

- `ArgumentError`: If `max_iterations` is lower than 1.
- `Sgp4FitDivergenceError`: If the least-square iterations diverge.

## Examples

```julia-repl
julia> tle = tle\"\"\"
          AMAZONIA 1
          1 47699U 21015A   23083.68657856  .00000000  00000-8  43000-3 0  9999
          2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652
          \"\"\"
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
                        B* :      0.00043 1 / er
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³

julia> update_sgp4_mean_elements_epoch(tle, DateTime("2023-06-19"))
ACTION:   Updating the epoch of the mean elements.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          2.79315e-06          2.57774e-09          2.79315e-06             -99.9999 %

TLE:
                      Name : AMAZONIA 1
          Satellite number : 47699
  International designator : 21015A
        Epoch (Year / Day) : 23 / 170.00000000 (2023-06-19T00:00:00)
        Element set number : 999
              Eccentricity :   0.00012727
               Inclination :  98.43040000 deg
                      RAAN : 247.20081879 deg
       Argument of perigee : 237.78423781 deg
              Mean anomaly : 121.83893462 deg
           Mean motion (n) :  14.41052177 revs / day
         Revolution number : 10865
                        B* :      0.00043 1 / er
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```
"""
function update_sgp4_mean_elements_epoch(
    me::_SGP4_MEAN_ELEMENTS,
    new_epoch::Union{Number, DateTime};
    sgp4c::Sgp4Constants = SGP4C_WGS84,
    kwargs...,
)
    sgp4d = Sgp4Propagator(sgp4c)
    return update_sgp4_mean_elements_epoch!(sgp4d, me, new_epoch; kwargs...)
end

"""
    update_sgp4_mean_elements_epoch!(
        sgp4d::Sgp4Propagator,
        me::S,
        new_epoch::Union{Number, DateTime};
        kwargs...,
    ) where {S <: Union{TLE, OrbitMeanElementsMessage}} -> S

Update the epoch of the SGP4 mean elements `me`, which can be a `TLE` or an
`OrbitMeanElementsMessage`, to `new_epoch`, represented by a Julian Day [UTC] or a
`DateTime` [UTC], using the orbit propagator `sgp4d`.

!!! note

    The SGP4 orbit propagator `sgp4d` will be initialized with the mean elements returned
    by the function.

For more information about the algorithm and the keywords, see
[`update_sgp4_mean_elements_epoch`](@ref), except for `sgp4c`, since the constants are
those in `sgp4d`. The function can fail if the least-square iterations diverge.

# Extended help

## Throws

- `ArgumentError`: If `max_iterations` is lower than 1.
- `Sgp4FitDivergenceError`: If the least-square iterations diverge.
"""
function update_sgp4_mean_elements_epoch!(
    sgp4d::Sgp4Propagator{Tepoch, T},
    me::S,
    new_epoch::Union{Number, DateTime};
    atol::Number                            = 2e-4,
    rtol::Number                            = 2e-4,
    jacobian_method::AbstractJacobianMethod = FiniteDiffJacobian(),
    jacobian_perturbation::Number           = 1e-3,
    jacobian_perturbation_tol::Number       = 1e-7,
    max_iterations::Int                     = 50,
    verbose::Bool                           = true,
) where {Tepoch <: Number, T <: Number, S <: _SGP4_MEAN_ELEMENTS}
    if max_iterations < 1
        throw(ArgumentError("The maximum number of iterations must be at least 1."))
    end

    # Unpack.
    sgp4c = sgp4d.sgp4c

    # New epoch of the mean elements [Julian Day].
    epoch = _julian_day(new_epoch)

    # Check if `stdout` supports colors. This flag also selects whether the progress line
    # is updated in place using terminal escape sequences.
    has_color = get(stdout, :color, false)::Bool

    verbose && _fit_print_action(has_color, "Updating the epoch of the mean elements.")

    # Convert the mean elements to the mean state vector and update its epoch.
    x = _update_sgp4_mean_state_vector!(
        sgp4d,
        _mean_state_vector(me; sgp4c = sgp4c),
        _mean_elements_epoch(me),
        epoch;
        atol,
        rtol,
        jacobian_method,
        jacobian_perturbation,
        jacobian_perturbation_tol,
        max_iterations,
        verbose,
        has_color,
    )

    # Build the mean elements keeping the metadata of the input.
    new_me = _build_mean_elements(
        S, x, epoch; sgp4c = sgp4c, template = me, covariance = nothing
    )

    # Initialize the propagator with the updated mean elements.
    sgp4_init!(sgp4d, new_me)

    return new_me
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _fit_sgp4_mean_state_vector!(
        sgp4d::Sgp4Propagator{Tepoch, T},
        vjd::AbstractVector,
        vy::AbstractVector{SVector{6, T}},
        x₁::SVector{7, T},
        epoch::Number,
        W::SVector{6, T};
        kwargs...,
    ) where {Tepoch <: Number, T <: Number} -> SVector{7, T}, SMatrix{7, 7, T}, NamedTuple

Fit the SGP4 mean state vector at `epoch` [Julian Day, UTC] using the propagator `sgp4d`
and the least-square algorithm in **[1]**, starting from the initial guess `x₁` and using
the osculating state vectors `vy` (position [km] and velocity [km / s] in TEME) measured at
the instants `vjd` [Julian Day, UTC] with the weights `W`. The mean state vector has the
following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

The propagator `sgp4d` is used as a workspace and its final state is undefined. The
function can fail if the least-square iterations diverge.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value.
- `rtol::Number`: Tolerance for the relative difference between the residues.
- `estimate_bstar::Bool`: If `true`, the drag term B* is estimated. Otherwise, it is kept
    with the value in `x₁`.
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix.
- `jacobian_perturbation::Number`: Initial state perturbation used with
    `FiniteDiffJacobian()`.
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation used with
    `FiniteDiffJacobian()`.
- `max_iterations::Int`: Maximum number of iterations.
- `verbose::Bool`: If `true`, the progress is printed to `stdout`.
- `has_color::Bool`: If `true`, the progress is printed using colors and updated in place.

# Returns

- `SVector{7, T}`: The fitted mean state vector.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm.
- `NamedTuple`: Statistics of the least-square algorithm with the fields `converged`,
    `iterations`, `position_rmse`, `velocity_rmse`, and `total_rmse` (see
    [`fit_sgp4_mean_elements!`](@ref)).

# Extended help

## Throws

- `Sgp4FitDivergenceError`: If the least-square iterations diverge.
"""
function _fit_sgp4_mean_state_vector!(
    sgp4d::Sgp4Propagator{Tepoch, T},
    vjd::AbstractVector,
    vy::AbstractVector{SVector{6, T}},
    x₁::SVector{7, T},
    epoch::Number,
    W::SVector{6, T};
    atol::Number,
    rtol::Number,
    estimate_bstar::Bool,
    jacobian_method::AbstractJacobianMethod,
    jacobian_perturbation::Number,
    jacobian_perturbation_tol::Number,
    max_iterations::Int,
    verbose::Bool,
    has_color::Bool,
) where {Tepoch <: Number, T <: Number}
    num_measurements = length(vjd)

    # Number of states in the input vector.
    num_states = 7

    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    x₂ = x₁

    # Variable to store the last residue. It is only read from the second iteration on.
    σ_i_₁ = T(0)

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Statistics returned after the iterations.
    converged  = false
    iterations = 0
    σ_i        = T(0)
    σp_i       = T(0)
    σv_i       = T(0)

    # Header.
    verbose && _fit_print_header(has_color)

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)

    # == Workspaces ========================================================================
    #
    # The propagator is initialized once per mean state vector and propagated to all
    # measurements, since the initialization is more expensive than the propagation. Hence,
    # we need buffers to store the nominal propagated state vectors and the Jacobians of
    # all measurements. They are allocated once per fit.

    # Pre-allocate the Dual-typed propagator for ForwardDiff Jacobian computation so it is
    # reused across all iterations instead of being heap-allocated on every call.
    sgp4d_ad =
        jacobian_method isa ForwardDiffJacobian ? _create_ad_propagator(sgp4d) : nothing

    # Nominal propagated state vectors [km, km / s].
    vŷ = Vector{SVector{6, T}}(undef, num_measurements)

    # Jacobians of all measurements, in which `vJ[:, :, k]` is the Jacobian of the k-th one.
    vJ = Array{T, 3}(undef, 6, num_states, num_measurements)

    # Loop until the maximum allowed iteration.
    @inbounds for it in 1:max_iterations
        x₁ = x₂
        iterations = it

        # Variables to store the summations to compute the least square fitting algorithm.
        ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)
        ΣJ′Wb = @SVector zeros(T, num_states)

        # Variable to store the RMS errors in this iteration.
        σ_i  = T(0)
        σp_i = T(0)
        σv_i = T(0)

        # == Nominal Propagation ===========================================================

        # Initialize the SGP4 with the current estimated mean elements and propagate the
        # orbit to all measurements.
        _init_sgp4_with_state_vector!(sgp4d, x₁, epoch)

        for k in 1:num_measurements
            Δt = (vjd[k] - epoch) * 1440
            r̂_teme, v̂_teme = sgp4!(sgp4d, Δt)
            vŷ[k] = vcat(r̂_teme, v̂_teme)
        end

        # == Jacobian ======================================================================

        _sgp4_jacobian!(
            jacobian_method,
            vJ,
            sgp4d,
            sgp4d_ad,
            vjd,
            epoch,
            x₁,
            vŷ;
            perturbation     = jacobian_perturbation,
            perturbation_tol = jacobian_perturbation_tol,
        )

        # == Accumulation ==================================================================

        for k in 1:num_measurements
            # Compute the residue.
            b = vy[k] - vŷ[k]

            # Obtain the Jacobian of this measurement.
            J = SMatrix{6, 7, T}(@view vJ[:, :, k])

            ΣJ′WJ += J' * (W .* J)
            ΣJ′Wb += J' * (W .* b)
            σ_i   += dot(W .* b, b)
            σp_i  += b[1]^2 + b[2]^2 + b[3]^2
            σv_i  += b[4]^2 + b[5]^2 + b[6]^2
        end

        # Normalize and compute the RMS errors.
        σ_i  = √(σ_i / num_measurements)
        σp_i = √(σp_i / num_measurements)
        σv_i = √(σv_i / num_measurements)

        # == Estimate Update ===============================================================

        if estimate_bstar
            δx = ΣJ′WJ \ ΣJ′Wb
        else
            ΣJ′WJ_sub = SMatrix{6, 6, T}(@view ΣJ′WJ[1:6, 1:6])
            ΣJ′Wb_sub = SVector{6, T}(@view ΣJ′Wb[1:6])
            δx_sub    = ΣJ′WJ_sub \ ΣJ′Wb_sub

            δx = SVector{7, T}(
                δx_sub[1], δx_sub[2], δx_sub[3], δx_sub[4], δx_sub[5], δx_sub[6], 0
            )
        end

        # Limit the correction to avoid divergence, but it should not be applied to B*.
        for i in 1:6
            threshold = T(0.1)
            if !iszero(x₁[i]) && (abs(δx[i] / x₁[i]) > threshold)
                δx = setindex(δx, threshold * abs(x₁[i]) * sign(δx[i]), i)
            end
        end

        x₂ = x₁ + δx

        # == Convergence Check =============================================================

        # We cannot compute the RMSE variation in the first iteration.
        if it == 1
            verbose && _fit_print_progress(
                has_color,
                @sprintf("%10d %20g %20g %20g %20s", it, σp_i, σv_i, σ_i, "---")
            )

        else
            # Compute the RMSE variation.
            Δσ = (σ_i - σ_i_₁) / σ_i_₁

            verbose && _fit_print_progress(
                has_color,
                @sprintf("%10d %20g %20g %20g %20g %%", it, σp_i, σv_i, σ_i, 100 * Δσ)
            )

            # Check if the RMSE is increasing.
            if σ_i < σ_i_₁
                Δd = 0
            else
                Δd += 1
            end

            # If the RMSE increased by three iterations and its value is higher than 5e11,
            # we abort because the iterations are diverging.
            ((Δd ≥ 3) && (σ_i > 5e11)) && throw(Sgp4FitDivergenceError(it, σ_i))

            # Check if the condition to stop has been reached.
            if (abs(Δσ) < rtol) || (σ_i < atol)
                converged = true
                break
            end
        end

        σ_i_₁ = σ_i
    end

    verbose && println()

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Assemble the statistics.
    stats = (;
        converged,
        iterations,
        position_rmse = σp_i,
        velocity_rmse = σv_i,
        total_rmse    = σ_i,
    )

    return x₂, P, stats
end

"""
    _update_sgp4_mean_state_vector!(
        sgp4d::Sgp4Propagator{Tepoch, T},
        x::SVector{7, T},
        epoch::Number,
        new_epoch::Number;
        kwargs...,
    ) where {Tepoch <: Number, T <: Number} -> SVector{7, T}

Update the SGP4 mean state vector `x` from `epoch` [Julian Day, UTC] to `new_epoch` [Julian
Day, UTC] using the propagator `sgp4d` as a workspace. The mean state vector has the same
elements as described in [`_fit_sgp4_mean_state_vector!`](@ref).

The orbit is propagated to `new_epoch` and a new mean state vector is fitted so that it
provides the same osculating state vector at the new epoch. The drag term B* is kept
constant. If the epochs differ by less than 1 ms, `x` is returned unchanged. The function
can fail if the least-square iterations diverge.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value.
- `rtol::Number`: Tolerance for the relative difference between the residues.
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations.
- `verbose::Bool`: If `true`, the progress is printed to `stdout`.
- `has_color::Bool`: If `true`, the progress is printed using colors and updated in place.

# Extended help

## Throws

- `Sgp4FitDivergenceError`: If the least-square iterations diverge.
"""
function _update_sgp4_mean_state_vector!(
    sgp4d::Sgp4Propagator{Tepoch, T},
    x::SVector{7, T},
    epoch::Number,
    new_epoch::Number;
    atol::Number,
    rtol::Number,
    jacobian_method::AbstractJacobianMethod = FiniteDiffJacobian(),
    jacobian_perturbation::Number = 1e-3,
    jacobian_perturbation_tol::Number = 1e-7,
    max_iterations::Int,
    verbose::Bool,
    has_color::Bool,
) where {Tepoch <: Number, T <: Number}
    # Do not update the epoch if the new epoch is less than 1 ms from the current one.
    abs(epoch - new_epoch) < 0.001 / 86400 && return x

    # Propagate up to the desired epoch.
    _init_sgp4_with_state_vector!(sgp4d, x, epoch)
    r_teme, v_teme = sgp4!(sgp4d, 1440 * (new_epoch - epoch))

    # Assemble the initial guess using the osculating state vector together with the
    # current B*, which is kept constant.
    x₁ = SVector{7, T}(r_teme..., v_teme..., x[7])

    # Now, we want to fit a mean state vector at the new epoch that provides the same
    # position and velocity vectors computed previously.
    vjd = SVector{1, typeof(new_epoch)}(new_epoch)
    vy  = SVector{1, SVector{6, T}}(vcat(r_teme, v_teme))
    W   = SVector{6, T}(1, 1, 1, 1, 1, 1)

    x₂, ~, ~ = _fit_sgp4_mean_state_vector!(
        sgp4d,
        vjd,
        vy,
        x₁,
        new_epoch,
        W;
        atol,
        rtol,
        estimate_bstar = false,
        jacobian_method,
        jacobian_perturbation,
        jacobian_perturbation_tol,
        max_iterations,
        verbose,
        has_color,
    )

    return x₂
end

"""
    _check_template(template::NamedTuple, reserved::Tuple) -> Nothing

Throw an `ArgumentError` if the `NamedTuple` `template` contains any of the fields in
`reserved`, which are set by the fitting algorithm.
"""
function _check_template(template::NamedTuple, reserved::Tuple)
    invalid = filter(k -> k in reserved, keys(template))

    isempty(invalid) || throw(
        ArgumentError(
            "The template must not contain the fields set by the fit: " *
            join(string.(invalid), ", ") *
            ".",
        ),
    )

    return nothing
end

"""
    _julian_day(epoch::Union{Number, DateTime}) -> Number

Return the `epoch` as a Julian Day, converting it if it is a `DateTime` [UTC].
"""
_julian_day(epoch::Number) = epoch
_julian_day(epoch::DateTime) = datetime2julian(epoch)

"""
    _closest_measurement(vjd::AbstractVector, epoch::Number) -> Int

Return the index of the instant in `vjd` [Julian Day] closest to `epoch` [Julian Day].
"""
function _closest_measurement(vjd::AbstractVector, epoch::Number)
    id = firstindex(vjd)
    v  = abs(vjd[id] - epoch)

    for k in eachindex(vjd)
        vk = abs(vjd[k] - epoch)

        if vk < v
            id = k
            v  = vk
        end
    end

    return id
end

"""
    _init_sgp4_with_state_vector!(
        sgp4d::Sgp4Propagator,
        sv::SVector{7},
        epoch::Number,
    ) -> Nothing

Initialize the SGP4 orbit propagator `sgp4d` using the state vector `sv`, which must have
the following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

and be defined for the `epoch` [Julian Day].
"""
function _init_sgp4_with_state_vector!(sgp4d::Sgp4Propagator, sv::SVector{7}, epoch::Number)
    # Unpack.
    sgp4c = sgp4d.sgp4c

    # Obtain the initial mean Keplerian elements.
    r_teme   = @SVector [1000sv[1], 1000sv[2], 1000sv[3]]
    v_teme   = @SVector [1000sv[4], 1000sv[5], 1000sv[6]]
    bstar    = sv[7]
    orb_teme = rv_to_kepler(r_teme, v_teme, epoch)

    # Obtain the required mean elements to initialize the SGP4.
    a₀ = orb_teme.a / (1000 * sgp4c.R0) # ............................. Semi-major axis [er]
    e₀ = orb_teme.e                     # ................................. Eccentricity [ ]
    i₀ = orb_teme.i                     # ................................ Inclination [rad]
    Ω₀ = orb_teme.Ω                     # ....................................... RAAN [rad]
    ω₀ = orb_teme.ω                     # ............................ Arg. of perigee [rad]
    f₀ = orb_teme.f                     # ............................... True anomaly [rad]
    M₀ = true_to_mean_anomaly(e₀, f₀)   # ............................... Mean anomaly [rad]

    # Obtain the mean motion [rad / min].
    n₀ = sgp4c.XKE / √(a₀^3)

    # Initialize the orbit propagator.
    sgp4_init!(sgp4d, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)

    return nothing
end

"""
    _mean_state_vector_to_elements(sv::SVector{7}, sgp4c::Sgp4Constants) -> NTuple{7, Float64}

Convert the mean state vector `sv` to the SGP4 mean elements using the constants `sgp4c`.
The state vector must have the following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

# Returns

- `Float64`: Mean motion [rev/day].
- `Float64`: Eccentricity [-].
- `Float64`: Inclination [°].
- `Float64`: Right ascension of the ascending node [°].
- `Float64`: Argument of perigee [°].
- `Float64`: Mean anomaly [°].
- `Float64`: Drag term B* [1 / er].
"""
function _mean_state_vector_to_elements(sv::SVector{7}, sgp4c::Sgp4Constants)
    r_teme   = @SVector [1000sv[1], 1000sv[2], 1000sv[3]]
    v_teme   = @SVector [1000sv[4], 1000sv[5], 1000sv[6]]
    orb_teme = rv_to_kepler(r_teme, v_teme)

    # Obtain the Keplerian elements with the units used by the TLE and the OMM.
    a₀ = orb_teme.a / (1000 * sgp4c.R0)
    e₀ = orb_teme.e
    i₀ = rad2deg(orb_teme.i)
    Ω₀ = rad2deg(orb_teme.Ω)
    ω₀ = rad2deg(orb_teme.ω)
    M₀ = rad2deg(true_to_mean_anomaly(e₀, orb_teme.f))

    # Obtain the mean motion [rad/min] and convert it to [rev/day].
    n₀ = sgp4c.XKE / √(a₀^3)

    # The TLE and the OMM store the elements as `Float64`. Hence, we convert them here so
    # that the fitting also works with other number types, e.g. `Float32`.
    return (
        Float64(720n₀ / π),
        Float64(e₀),
        Float64(i₀),
        Float64(Ω₀),
        Float64(ω₀),
        Float64(M₀),
        Float64(sv[7]),
    )
end

"""
    _elements_to_mean_state_vector(
        n₀::Number,
        e₀::Number,
        i₀::Number,
        Ω₀::Number,
        ω₀::Number,
        M₀::Number,
        bstar::Number,
        sgp4c::Sgp4Constants{T}
    ) where {T <: Number} -> SVector{7, T}

Convert the SGP4 mean elements to the mean state vector using the constants `sgp4c`. The
elements are the mean motion `n₀` [rev/day], the eccentricity `e₀` [-], the inclination
`i₀` [°], the right ascension of the ascending node `Ω₀` [°], the argument of perigee `ω₀`
[°], the mean anomaly `M₀` [°], and the drag term `bstar` [1 / er]. The state vector has
the following structure:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘
"""
function _elements_to_mean_state_vector(
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number,
    sgp4c::Sgp4Constants{T},
) where {T <: Number}
    # Convert the mean motion to [rad/min] and obtain the semi-major axis [m].
    n = n₀ * π / 720
    a = (1000 * sgp4c.R0) * (sgp4c.XKE / n)^(2 / 3)
    f = mean_to_true_anomaly(e₀, deg2rad(M₀))

    # Convert to state vector.
    orb = KeplerianElements(0, a, e₀, deg2rad(i₀), deg2rad(Ω₀), deg2rad(ω₀), f)
    r_i, v_i = kepler_to_rv(orb)

    return SVector{7, T}(
        r_i[1] / 1000,
        r_i[2] / 1000,
        r_i[3] / 1000,
        v_i[1] / 1000,
        v_i[2] / 1000,
        v_i[3] / 1000,
        bstar,
    )
end

# == Jacobian ==============================================================================

"""
    _create_ad_propagator(
        sgp4d::Sgp4Propagator{Tepoch, T},
    ) where {Tepoch <: Number, T <: Number} -> Sgp4Propagator

Create a Dual-typed SGP4 propagator for use with `ForwardDiffJacobian`. The returned
propagator is passed to [`_sgp4_jacobian!`](@ref) so that it is allocated once per fit.
"""
function _create_ad_propagator(sgp4d::Sgp4Propagator{Tepoch, T}) where {Tepoch, T}
    tag   = ForwardDiff.Tag{Nothing, T}
    D     = ForwardDiff.Dual{tag, T, 7}
    sgp4c = sgp4d.sgp4c

    return Sgp4Propagator{Tepoch}(Sgp4Constants{D}(sgp4c))
end

"""
    _sgp4_jacobian!(
        ::FiniteDiffJacobian,
        vJ::AbstractArray{T, 3},
        sgp4d::Sgp4Propagator{Tepoch, T},
        sgp4d_ad::Nothing,
        vjd::AbstractVector,
        epoch::Number,
        x₁::SVector{7, T},
        vŷ::AbstractVector{SVector{6, T}};
        kwargs...,
    ) where {Tepoch <: Number, T <: Number} -> Nothing

Compute by finite differences the SGP4 Jacobians with respect to the mean state vector `x₁`
at `epoch` [Julian Day] for all the measurement instants `vjd` [Julian Day], storing the
Jacobian of the k-th measurement in `vJ[:, :, k]`. The vector `vŷ` must contain the state
vectors propagated with `x₁` to the instants `vjd`. Hence:

                 ∂sgp4(x, Δt) │
    vJ[:, :, k] = ──────────── │
                      ∂x      │ x = x₁, Δt = vjd[k] - epoch

The propagator `sgp4d` is used as a workspace and its final state is undefined. The
propagator `sgp4d_ad` is not used.

# Keywords

- `perturbation::Number`: Initial state perturbation to compute the finite-difference:
    `Δx = x * perturbation`.
    (**Default**: 1e-3)
- `perturbation_tol::Number`: Tolerance to accept the perturbation. If the computed
    perturbation is lower than `perturbation_tol`, we increase it until its absolute value
    is higher than `perturbation_tol`.
    (**Default**: 1e-7)
"""
function _sgp4_jacobian!(
    ::FiniteDiffJacobian,
    vJ::AbstractArray{T, 3},
    sgp4d::Sgp4Propagator{Tepoch, T},
    sgp4d_ad::Nothing,
    vjd::AbstractVector,
    epoch::Number,
    x₁::SVector{7, T},
    vŷ::AbstractVector{SVector{6, T}};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
) where {Tepoch <: Number, T <: Number}
    num_measurements = length(vjd)

    @inbounds for j in 1:7
        # State that will be perturbed.
        α = x₁[j]

        # Obtain the perturbation, taking care to avoid small values.
        ϵ = α * T(perturbation)

        for _ in 1:5
            abs(ϵ) > perturbation_tol && break
            ϵ *= T(1.4)
        end

        # Avoid division by zero in cases that α is very small. In this situation, we force
        # `|ϵ| = perturbation_tol`.
        if abs(ϵ) < perturbation_tol
            ϵ = signbit(α) ? -perturbation_tol : perturbation_tol
        end

        # Initialize the propagator with the perturbed state and propagate it to all the
        # measurements to obtain the j-th column of every Jacobian.
        _init_sgp4_with_state_vector!(sgp4d, setindex(x₁, α + ϵ, j), epoch)

        for k in 1:num_measurements
            Δt = (vjd[k] - epoch) * 1440
            r_teme, v_teme = sgp4!(sgp4d, Δt)
            y₂ = vcat(r_teme, v_teme)
            ∂y = (y₂ - vŷ[k]) / ϵ

            for i in 1:6
                vJ[i, j, k] = ∂y[i]
            end
        end
    end

    return nothing
end

"""
    _sgp4_jacobian!(
        ::ForwardDiffJacobian,
        vJ::AbstractArray{T, 3},
        sgp4d::Sgp4Propagator{Tepoch, T},
        sgp4d_ad::Sgp4Propagator{Tepoch, D},
        vjd::AbstractVector,
        epoch::Number,
        x₁::SVector{7, T},
        vŷ::AbstractVector{SVector{6, T}};
        kwargs...,
    ) where {Tepoch <: Number, T <: Number, D <: ForwardDiff.Dual} -> Nothing

Compute by forward-mode automatic differentiation the SGP4 Jacobians with respect to the
mean state vector `x₁` at `epoch` [Julian Day] for all the measurement instants `vjd`
[Julian Day], storing the Jacobian of the k-th measurement in `vJ[:, :, k]`. Hence:

                 ∂sgp4(x, Δt) │
    vJ[:, :, k] = ──────────── │
                      ∂x      │ x = x₁, Δt = vjd[k] - epoch

The Dual-typed propagator `sgp4d_ad`, created by [`_create_ad_propagator`](@ref), is used
as a workspace and its final state is undefined. The propagator `sgp4d` and the propagated
state vectors `vŷ` are not used. The keywords are accepted for compatibility with the
finite-difference method and are ignored.
"""
function _sgp4_jacobian!(
    ::ForwardDiffJacobian,
    vJ::AbstractArray{T, 3},
    sgp4d::Sgp4Propagator{Tepoch, T},
    sgp4d_ad::Sgp4Propagator{Tepoch, D},
    vjd::AbstractVector,
    epoch::Number,
    x₁::SVector{7, T},
    vŷ::AbstractVector{SVector{6, T}};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
) where {Tepoch <: Number, T <: Number, D <: ForwardDiff.Dual}
    num_measurements = length(vjd)

    # Seed the dual numbers so that the j-th partial is the derivative with respect to the
    # j-th element of the mean state vector.
    seeds  = ntuple(i -> ForwardDiff.Partials(ntuple(j -> T(i == j), Val(7))), Val(7))
    x_dual = SVector{7, D}(ntuple(i -> D(x₁[i], seeds[i]), Val(7)))

    # Initialize the propagator once and propagate it to all the measurements.
    _init_sgp4_with_state_vector!(sgp4d_ad, x_dual, epoch)

    @inbounds for k in 1:num_measurements
        Δt     = (vjd[k] - epoch) * 1440
        r, v   = sgp4!(sgp4d_ad, Δt)
        y_dual = vcat(r, v)

        for j in 1:7, i in 1:6
            vJ[i, j, k] = ForwardDiff.partials(y_dual[i], j)
        end
    end

    return nothing
end

# == Printing Helpers ======================================================================

"""
    _fit_print_action(has_color::Bool, msg::AbstractString) -> Nothing

Print to `stdout` the action message `msg` of the fitting algorithm, prefixed by an
`ACTION:` tag, which is highlighted if `has_color` is `true`.
"""
# The helper is not inlined so that its allocation sites, which are only reachable when the
# algorithm is verbose, are counted once regardless of the number of call sites.
@noinline function _fit_print_action(has_color::Bool, msg::AbstractString)
    println(_fit_decorated(_FIT_ACTION_TAG, has_color), "   ", msg)
    return nothing
end

"""
    _fit_print_header(has_color::Bool) -> Nothing

Print to `stdout` the header of the progress table of the fitting algorithm. The header is
decorated if `has_color` is `true`. In this case, it is followed by an empty line that the
first progress line overwrites.
"""
# The helper is not inlined so that its allocation sites, which are only reachable when the
# algorithm is verbose, are counted once regardless of the number of call sites.
@noinline function _fit_print_header(has_color::Bool)
    print(
        "          ",
        _fit_decorated(_FIT_HEADER, has_color),
        "\n          ",
        _fit_decorated(_FIT_UNITS, has_color),
        "\n",
    )

    # The empty line is only required if the progress line is updated in place.
    has_color && println()

    return nothing
end

"""
    _fit_print_progress(has_color::Bool, msg::AbstractString) -> Nothing

Print to `stdout` the progress line `msg` of the fitting algorithm, prefixed by a
`PROGRESS:` tag. If `has_color` is `true`, the tag is highlighted and the previous line is
erased first using terminal escape sequences, so consecutive calls update the same terminal
line. Otherwise, each call prints a new line, keeping the output readable when it is
redirected to a file.
"""
# The helper is not inlined so that its allocation sites, which are only reachable when the
# algorithm is verbose, are counted once regardless of the number of call sites.
@noinline function _fit_print_progress(has_color::Bool, msg::AbstractString)
    has_color && print("\x1b[A\x1b[2K\r")
    println(_fit_decorated(_FIT_PROGRESS_TAG, has_color), " ", msg)
    return nothing
end

"""
    _fit_decorated(versions::Tuple{String, String}, has_color::Bool) -> String

Return the colored version of a string printed by the fitting algorithm, stored as the
second element of `versions`, if `has_color` is `true`. Otherwise, return the plain version
stored as its first element.
"""
function _fit_decorated(versions::Tuple{String, String}, has_color::Bool)
    return versions[has_color ? 2 : 1]
end
