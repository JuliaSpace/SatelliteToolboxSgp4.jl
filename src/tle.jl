# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions related to SGP4 TLEs.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute of
#       Aeronautics ans Astronautics.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export fit_sgp4_tle, fit_sgp4_tle!, update_sgp4_tle_epoch, update_sgp4_tle_epoch!

const _INITIAL_GUESS_T = Union{Nothing, AbstractVector, TLE}

"""
    fit_sgp4_tle(vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> TLE, SMatrix{7, 7, T}

Fit a Two-Line Element set (`TLE`) for the SGP4 orbit propagator using the osculating
elements represented by a set of position vectors `vr_teme` [km] and a set of velocity
vectors `vv_teme` [km / s] represented in the True-Equator, Mean-Equinox reference frame
(TEME) at instants in the array `vjd` [Julian Day].

This algorithm was based on **[1]**.

!!! note
    This algorithm version will allocate a new SGP4 propagator with the default constants
    `sgp4c_wgs84`. If another set of constants are required, use the function
    [`fit_sgp4_tle!`](@ref) instead.

# Keywords

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

# Returns

- `TLE`: The fitted TLE.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm.

# Initial Guess

This algorithm uses a least-square algorithm to fit a TLE based on a set of osculating state
vectors. Since the system is chaotic, a good initial guess is paramount for algorithm
convergence. We can provide an initial guess using the keyword `initial_guess`.

If `initial_guess` is a `TLE`, we update the TLE epoch using the function
[`update_sgp4_tle_epoch!`](@ref) to the desired one in `mean_elements_epoch`. Afterward, we
use this new TLE as the initial guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted TLE is
obtained, the algorithm uses the function [`update_sgp4_tle_epoch!`](@ref) to change its
epoch to `mean_elements_epoch`.

!!! note
    If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the TLE or
    the state vector. Hence, if `estimate_bstar` is `false`, it will be kept constant with
    this initial value.

# Examples

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

julia> tle, P = fit_sgp4_tle(vjd, vr_teme, vv_teme; estimate_bstar = false)
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          1               13.589           0.00473939               13.589                  ---
PROGRESS:          2          0.000902827          4.59551e-06          0.000902839             -99.9934 %
PROGRESS:          3          5.80079e-09          4.38304e-07          4.38343e-07             -99.9514 %

(TLE: UNDEFINED (Epoch = 2023-03-24T16:33:40.388), [1.0855427857632278 -0.024689551085771335 … -0.07505082051235405 -5691.217934764823; -0.024689551096333886 1.0180284581267744 … 0.023015858293549355 1745.6307022631847; … ; -0.07505082051312215 0.023015858293646975 … 0.06515845129134726 4946.150014962679; -5691.217934824044 1745.6307022705219 … 4946.150014962695 3.7558624425028527e8])

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
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics ans Astronautics.
"""
function fit_sgp4_tle(
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    # We must initialize the SGP4 propagator structure together with any mutable fields.
    sgp4d = Sgp4Propagator{Float64, Float64}()
    sgp4d.sgp4c = sgp4c_wgs84
    sgp4d.sgp4ds = Sgp4DeepSpace{Float64}()

    return fit_sgp4_tle!(sgp4d, vjd, vr_teme, vv_teme; kwargs...)
end

"""
    fit_sgp4_tle!(sgp4d::Sgp4Propagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> TLE, SMatrix{7, 7, T}

Fit a Two-Line Element set (`TLE`) for the SGP4 orbit propagator `sgp4d` using the
osculating elements represented by a set of position vectors `vr_teme` [km] and a set of
velocity vectors `vv_teme` [km / s] represented in the True-Equator, Mean-Equinox reference
frame (TEME) at instants in the array `vjd` [Julian Day].

This algorithm was based on **[1]**.

!!! notes
    The SGP4 orbit propagator `sgp4d` will be initialized with the TLE returned by the
    function.

# Keywords

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

# Returns

- `TLE`: The fitted TLE.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm.

# Initial Guess

This algorithm uses a least-square algorithm to fit a TLE based on a set of osculating state
vectors. Since the system is chaotic, a good initial guess is paramount for algorithm
convergence. We can provide an initial guess using the keyword `initial_guess`.

If `initial_guess` is a `TLE`, we update the TLE epoch using the function
[`update_sgp4_tle_epoch!`](@ref) to the desired one in `mean_elements_epoch`. Afterward, we
use this new TLE as the initial guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted TLE is
obtained, the algorithm uses the function [`update_sgp4_tle_epoch!`](@ref) to change its
epoch to `mean_elements_epoch`.

!!! note
    If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the TLE or
    the state vector. Hence, if `estimate_bstar` is `false`, it will be kept constant with
    this initial value.

# Examples

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

julia> tle, P = fit_sgp4_tle!(sgp4d, vjd, vr_teme, vv_teme; estimate_bstar = false)
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          1               13.589           0.00473939               13.589                  ---
PROGRESS:          2          0.000902827          4.59551e-06          0.000902839             -99.9934 %
PROGRESS:          3          5.80079e-09          4.38304e-07          4.38343e-07             -99.9514 %

(TLE: UNDEFINED (Epoch = 2023-03-24T16:33:40.388), [1.0855427857632278 -0.024689551085771335 … -0.07505082051235405 -5691.217934764823; -0.024689551096333886 1.0180284581267744 … 0.023015858293549355 1745.6307022631847; … ; -0.07505082051312215 0.023015858293646975 … 0.06515845129134726 4946.150014962679; -5691.217934824044 1745.6307022705219 … 4946.150014962695 3.7558624425028527e8])

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
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics ans Astronautics.
"""
function fit_sgp4_tle!(
    sgp4d::Sgp4Propagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    atol::Number                      = 2e-4,
    rtol::Number                      = 2e-4,
    estimate_bstar::Bool              = true,
    initial_guess::_INITIAL_GUESS_T   = nothing,
    jacobian_perturbation::Number     = 1e-3,
    jacobian_perturbation_tol::Number = 1e-7,
    max_iterations::Int               = 50,
    mean_elements_epoch::Number       = vjd[end],
    verbose::Bool                     = true,
    weight_vector::AbstractVector     = @SVector(ones(Bool, 6)),
    # Keywords with the fields of the output TLE.
    classification::Char              = 'U',
    element_set_number::Int           = 0,
    international_designator::String  = "999999",
    name::String                      = "UNDEFINED",
    revolution_number::Int            = 0,
    satellite_number::Int             = 9999,
) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector}
    # Unpack.
    sgp4c = sgp4d.sgp4c

    # Number of available measurements.
    num_measurements = length(vjd)

    # Check the inputs.
    length(vr_teme) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vr_teme` must be the same."))

    length(vv_teme) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vv_teme` must be the same."))

    if length(weight_vector) != 6
        throw(ArgumentError("The weight vector must have 6 elements."))
    end

    if (initial_guess isa AbstractVector) && (length(initial_guess) != 7)
        throw(ArgumentError("The initial guess state vector must have 7 elements."))
    end

    # Check if stdout supports colors.
    has_color = get(stdout, :color, false)::Bool
    cd = has_color ? string(_D) : ""
    cb = has_color ? string(_B) : ""
    cy = has_color ? string(_Y) : ""

    # Assemble the weight matrix.
    W = Diagonal(
        @SVector T[
            weight_vector[1],
            weight_vector[2],
            weight_vector[3],
            weight_vector[4],
            weight_vector[5],
            weight_vector[6]
        ]
    )

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    if !isnothing(initial_guess)
        epoch = mean_elements_epoch

        # The user can provide a TLE or the initial mean state vector.
        if initial_guess isa TLE
            verbose && println("$(cy)ACTION:$(cd)   Updating the epoch of the initial TLE guess to match the desired one.")

            # If a TLE is provided, we need to update its epoch to match the desired one.
            tle = update_sgp4_tle_epoch!(
                sgp4d,
                initial_guess,
                epoch;
                verbose = verbose
            )

            # Now, we can convert the TLE to the mean state vector.
            x₁ = _tle_to_mean_state_vector(tle; sgp4c = sgp4c)

        else
            # In this case, the user must ensure that the provided mean elements are related
            # to the selected `mean_elements_epoch`.
            x₁ = SVector{7, T}(initial_guess...)

        end
    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        ~, id = findmin(abs.(vjd .- mean_elements_epoch))
        epoch = vjd[id]
        x₁ = SVector{7, T}(vr_teme[id]..., vv_teme[id]..., estimate_bstar ? T(0.00001) : T(0))
    end

    x₂ = x₁

    # Number of states in the input vector.
    num_states = 7

    # Number of observations in each instant.
    num_observations = 6

    # Covariance matrix.
    P = SMatrix{num_states, num_states, T}(I)

    # Variable to store the last residue.
    σ_i_₁ = nothing

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Allocate the Jacobian matrix.
    J = zeros(T, num_observations, num_states)

    # Header.
    verbose && println("$(cy)ACTION:$(cd)   Fitting the TLE.")
    verbose && @printf("          %s%10s %20s %20s %20s %20s%s\n", cy, "Iteration", "Position RMSE", "Velocity RMSE", "Total RMSE", "RMSE Variation", cd)
    verbose && @printf("          %s%10s %20s %20s %20s %20s%s\n", cb, "", "[km]", "[km / s]", "[ ]", "", cd)

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    ΣJ′WJ = nothing

    # Loop until the maximum allowed iteration.
    @inbounds for it in 1:max_iterations
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting algorithm.
        ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)
        ΣJ′Wb = @SVector zeros(T, num_states)

        # Variable to store the RMS errors in this iteration.
        σ_i  = T(0)
        σp_i = T(0)
        σv_i = T(0)

        for k in 1:num_measurements
            # Obtain the measured ephemerides.
            y = vcat(vr_teme[k], vv_teme[k])

            # Initialize the SGP4 with the current estimated mean elements.
            _init_sgp4_with_state_vector!(sgp4d, x₁, epoch)

            # Obtain the propagation time for this measurement.
            Δt = (vjd[k] - epoch) * 1440

            # Propagate the orbit.
            r̂_teme, v̂_teme = sgp4!(sgp4d, Δt)
            ŷ = vcat(r̂_teme, v̂_teme)

            # Compute the residue.
            b = y - ŷ

            # Compute the Jacobian inplace.
            _sgp4_jacobian!(
                J,
                sgp4d,
                Δt,
                x₁,
                ŷ;
                perturbation     = jacobian_perturbation,
                perturbation_tol = jacobian_perturbation_tol
            )

            # Convert the Jacobian matrix to a static matrix, leading to substantial
            # performance gains in the following computation.
            Js = SMatrix{num_observations, num_states, T}(J)

            # Accumulation.
            ΣJ′WJ += Js' * W * Js
            ΣJ′Wb += Js' * W * b
            σ_i   += b'  * W * b
            σp_i  += @views b[1:3]' * b[1:3]
            σv_i  += @views b[4:6]' * b[4:6]
        end

        # Normalize and compute the RMS errors.
        σ_i  = √(σ_i  / num_measurements)
        σp_i = √(σp_i / num_measurements)
        σv_i = √(σv_i / num_measurements)

        # Update the estimate.
        @views if estimate_bstar
            δx = ΣJ′WJ \ ΣJ′Wb
        else
            ΣJ′WJ_sub = SMatrix{6, 6, T}(ΣJ′WJ[1:6, 1:6])
            ΣJ′Wb_sub = SVector{6, T}(ΣJ′Wb[1:6])
            δx_sub    = ΣJ′WJ_sub \ ΣJ′Wb_sub

            δx  = SVector{7, T}(
                δx_sub[1],
                δx_sub[2],
                δx_sub[3],
                δx_sub[4],
                δx_sub[5],
                δx_sub[6],
                0
            )
        end

        # Limit the correction to avoid divergence, but it should not be applied to B*.
        for i in 1:6
            threshold = T(0.1)
            if abs(δx[i] / x₁[i]) > threshold
                δx = setindex(δx, threshold * abs(x₁[i]) * sign(δx[i]), i)
            end
        end

        x₂ = x₁ + δx

        # We cannot compute the RMSE variation in the first iteration.
        if it == 1
            verbose && @printf("%sPROGRESS:%s %10d %20g %20g %20g %20s\n", cb, cd, it, σp_i, σv_i, σ_i, "---")

        else
            # Compute the RMSE variation.
            Δσ = (σ_i - σ_i_₁) / σ_i_₁

            verbose && @printf("%sPROGRESS:%s %10d %20g %20g %20g %20g %%\n", cb, cd, it, σp_i, σv_i, σ_i, 100 * Δσ)

            # Check if the RMSE is increasing.
            if σ_i < σ_i_₁
                Δd = 0
            else
                Δd += 1
            end

            # If the RMSE increased by three iterations and its value is higher than 5e11,
            # we abort because the iterations are diverging.
            ((Δd ≥ 3) && (σ_i > 5e11)) && error("The iterations diverged!")

            # Check if the condition to stop has been reached.
            ((abs(Δσ) < rtol) || (σ_i < atol) || (it ≥ max_iterations)) && break
        end

        σ_i_₁ = σ_i
    end

    verbose && println()

    # Convert the state vector to TLE.
    tle = _mean_state_vector_to_tle(
        x₂,
        epoch;
        sgp4c                    = sgp4c,
        classification           = classification,
        element_set_number       = element_set_number,
        international_designator = international_designator,
        name                     = name,
        revolution_number        = revolution_number,
        satellite_number         = satellite_number,
    )

    # Update the epoch of the fitted TLE to match the desired one.
    if abs(epoch - mean_elements_epoch) > 0.001 / 86400
        verbose && println("$(cy)ACTION:$(cd)   Updating the epoch of the fitted TLE to match the desired one.")
        tle = update_sgp4_tle_epoch!(sgp4d, tle, mean_elements_epoch; verbose = verbose)
    end

    # Initialize the propagator with the TLE.
    sgp4_init!(sgp4d, tle)

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Return the TLE and the covariance.
    return tle, P
end

"""
    update_sgp4_tle_epoch(tle::TLE, new_epoch::Union{Number, DateTime}; kwargs...) -> TLE

Update the `tle` epoch with SGP4 mean elements to `new_epoch`. `new_epoch` can be
represented by a Julian Day or a `DateTime`.

!!! notes
    This algorithm version will allocate a new SGP4 propagator with the default constants
    `sgp4c_wgs84`. If another set of constants are required, use the function
    [`update_sgp4_tle_epoch!`](@ref) instead.

This function uses the following algorithm to update the TLE epoch:

1. Initialize the SGP4 propagator with `tle`;
2. Propagate the orbit to `new_epoch` and obtain the osculating state vector in TEME
    reference frame; and
3. Fit a new TLE that provides the same osculating state vector but considering the new
    epoch.

The third step uses the function [`fit_sgp4_tle!`](@ref). Hence, some keywords are related
to it. For more information, see the documentation of [`fit_sgp4_tle!`](@ref).

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If, at any iteration, the
    residue is lower than `atol`, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If, at any
    iteration, the relative difference between the residues in two consecutive iterations is
    lower than `rtol`, the computation loop stops. (**Default** = 2e-4)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)

# Examples

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
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³

julia> update_sgp4_tle_epoch(tle, DateTime("2023-06-19"))
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          1               15.085           0.00784051               15.085                  ---
PROGRESS:          2              8.69233           1.2628e-05              8.69233             -42.3778 %
PROGRESS:          3              3.33201           4.8632e-06              3.33201             -61.6673 %
PROGRESS:          4          2.79394e-06          2.57745e-09          2.79395e-06             -99.9999 %

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
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```
"""
function update_sgp4_tle_epoch(tle::TLE, new_epoch::Union{Number, DateTime}; kwargs...)
    # We must initialize the SGP4 propagator structure together with any mutable fields.
    sgp4d = Sgp4Propagator{Float64, Float64}()
    sgp4d.sgp4c = sgp4c_wgs84
    sgp4d.sgp4ds = Sgp4DeepSpace{Float64}()

    return update_sgp4_tle_epoch!(sgp4d, tle, new_epoch; kwargs...)
end

"""
    update_sgp4_tle_epoch!(sgp4d::Sgp4Propagator, tle::TLE, new_epoch::Union{Number, DateTime}; kwargs...) -> TLE

Update the `tle` epoch with SGP4 mean elements to `new_epoch` using the orbit propagator
`sgp4d`. `new_epoch` can be represented by a Julian Day or a `DateTime`.

!!! notes
    The SGP4 orbit propagator `sgp4d` will be initialized with the TLE returned by the
    function.

This function uses the following algorithm to update the TLE epoch:

1. Initialize the SGP4 propagator with `tle`;
2. Propagate the orbit to `new_epoch` and obtain the osculating state vector in TEME
    reference frame; and
3. Fit a new TLE that provides the same osculating state vector but considering the new
    epoch.

The third step uses the function [`fit_sgp4_tle!`](@ref). Hence, some keywords are related
to it. For more information, see the documentation of [`fit_sgp4_tle!`](@ref).

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If, at any iteration, the
    residue is lower than `atol`, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If, at any
    iteration, the relative difference between the residues in two consecutive iterations is
    lower than `rtol`, the computation loop stops. (**Default** = 2e-4)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)

# Examples

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
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³

# Allocate a new SGP4 orbit propagator using the created TLE. Notice that any TLE can be
# used here.
julia> sgp4d = sgp4_init(tle)

julia> update_sgp4_tle_epoch!(sgp4d, tle, DateTime("2023-06-19"))
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          1               15.085           0.00784051               15.085                  ---
PROGRESS:          2              8.69233           1.2628e-05              8.69233             -42.3778 %
PROGRESS:          3              3.33201           4.8632e-06              3.33201             -61.6673 %
PROGRESS:          4          2.79394e-06          2.57745e-09          2.79395e-06             -99.9999 %

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
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```
"""
function update_sgp4_tle_epoch!(sgp4d::Sgp4Propagator, tle::TLE, new_epoch::DateTime; kwargs...)
    dt = datetime2julian(new_epoch)
    return update_sgp4_tle_epoch!(sgp4d, tle, dt; kwargs...)
end

function update_sgp4_tle_epoch!(
    sgp4d::Sgp4Propagator,
    tle::TLE,
    new_epoch::Number;
    atol::Number           = 2e-4,
    rtol::Number           = 2e-4,
    max_iterations::Number = 50,
    verbose::Bool          = true
)
    # First, we need to initialize the SGP4 propagator with the initial TLE.
    sgp4_init!(sgp4d, tle)

    # Do not update the epoch if the new epoch is less than 1ms from the desired one.
    if abs(sgp4d.epoch - new_epoch) < 0.001 / 86400
        return tle
    end

    # Propagate up to the desired epoch.
    r_teme, v_teme = sgp4!(sgp4d, 1440 * (new_epoch - sgp4d.epoch))

    # Assemble the initial guess vector using the osculating information together with the
    # B* in the input TLE.
    initial_guess = SVector{7, eltype(r_teme)}(r_teme..., v_teme..., tle.bstar)

    # Now, we want to fit a TLE in the new epoch that provides the same position and
    # velocity vector computed previously.

    # NOTE: We are asking to not estimate B* since we are providing the same value as
    # obtained in the input TLE through the initial guess. Here we assume that B* does not
    # change.
    vjd     = @SVector [new_epoch]
    vr_teme = @SVector [r_teme]
    vv_teme = @SVector [v_teme]

    tle, ~ = fit_sgp4_tle!(
        sgp4d,
        vjd,
        vr_teme,
        vv_teme;
        atol                     = atol,
        rtol                     = rtol,
        max_iterations           = max_iterations,
        mean_elements_epoch      = new_epoch,
        verbose                  = verbose,
        classification           = tle.classification,
        element_set_number       = tle.element_set_number,
        estimate_bstar           = false,
        initial_guess            = initial_guess,
        international_designator = tle.international_designator,
        name                     = tle.name,
        revolution_number        = tle.revolution_number,
        satellite_number         = tle.satellite_number
    )

    return tle
end

############################################################################################
#                                    Private Functions
############################################################################################

"""
    _init_sgp4_with_state_vector!(sgp4d::Sgp4Propagator, sv::SVector{8}, epoch::Number) -> Nothing

Initialize the SGP4 orbit propagator `sgp4d` using the state vector `sv`, which must have
the following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

and be defined for the `epoch` [Julian Day].
"""
function _init_sgp4_with_state_vector!(
    sgp4d::Sgp4Propagator,
    sv::SVector{7},
    epoch::Number,
)
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
    _mean_state_vector_to_tle(sv::SVector{8}; kwargs...) -> TLE

Create a TLE given the mean state vector `sv`, which must have the following elements:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 propagator constants. (**Default** = sgp4c_wgs84)
- `classification::Char`: Satellite classification for the TLE. (**Default** = 'U')
- `element_set_number::Int`: Satellite element set number for the TLE. (**Default** = 0)
- `international_designator::String`: Satellite international designator for the TLE.
    (**Default** = "999999")
- `name::String`: Satellite name for the TLE. (**Default** = "UNDEFINED")
- `revolution_number::Int`: Satellite revolution number for the TLE. (**Default** = 0)
- `satellite_number::Int`: Satellite number for the TLE. (**Default** = 9999)
"""
function _mean_state_vector_to_tle(
    sv::SVector{7},
    epoch::Number;
    sgp4c::Sgp4Constants = sgp4c_wgs84,
    classification::Char = 'U',
    element_set_number::Int = 0,
    international_designator::String = "999999",
    name::String = "UNDEFINED",
    revolution_number::Int = 0,
    satellite_number::Int = 9999,
)

    r_teme   = @SVector [1000sv[1], 1000sv[2], 1000sv[3]]
    v_teme   = @SVector [1000sv[4], 1000sv[5], 1000sv[6]]
    bstar    = sv[7]
    orb_teme = rv_to_kepler(r_teme, v_teme)

    # Compute the data as required by the TLE format.
    dt  = julian2datetime(epoch)
    dt₀ = DateTime(Year(dt))

    dt_year    = year(dt)
    epoch_year = dt_year < 1980 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt - dt₀).value / 1000 / 86400 + 1

    # Obtain the Keplerian elements with the correct units for the TLE.
    a₀ = orb_teme.a / (1000 * sgp4c.R0)
    e₀ = orb_teme.e
    i₀ = rad2deg(orb_teme.i)
    Ω₀ = rad2deg(orb_teme.Ω)
    ω₀ = rad2deg(orb_teme.ω)
    M₀ = rad2deg(true_to_mean_anomaly(e₀, orb_teme.f))

    # Obtain the mean motion [rad/min].
    n₀ = sgp4c.XKE / √(a₀^3)

    # Construct the TLE.
    tle = TLE(
        name,
        satellite_number,
        classification,
        international_designator,
        epoch_year,
        epoch_day,
        0,
        0,
        bstar,
        element_set_number,
        i₀,
        Ω₀,
        e₀,
        ω₀,
        M₀,
        720n₀ / π,
        revolution_number
    )

    return tle
end

"""
    _tle_to_mean_state_vector(tle::TLE; sgp4c::Sgp4Constants = sgp4c_wgs84) -> SVector{7, T}

Convert the `tle` to a mean state vector using the SGP4 constants `sgp4c`. The state vector
has the following structure:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘
"""
function _tle_to_mean_state_vector(
    tle::TLE;
    sgp4c::Sgp4Constants{T} = sgp4c_wgs84
) where T<:Number
    # Unpack information.
    n = tle.mean_motion * π / 720
    a = (1000 * sgp4c.R0) * (sgp4c.XKE / n)^(2 / 3)
    e = tle.eccentricity
    i = deg2rad(tle.inclination)
    Ω = deg2rad(tle.raan)
    ω = deg2rad(tle.argument_of_perigee)
    M = deg2rad(tle.mean_anomaly)
    f = mean_to_true_anomaly(e, M)

    # Convert to state vector.
    orb = KeplerianElements(0, a, e, i, Ω, ω, f)
    r_i, v_i = kepler_to_rv(orb)

    sv = SVector{7, T}(
        r_i[1] / 1000,
        r_i[2] / 1000,
        r_i[3] / 1000,
        v_i[1] / 1000,
        v_i[2] / 1000,
        v_i[3] / 1000,
        tle.bstar
    )

    return sv
end

"""
    _sgp4_jacobian!(J::AbstractMatrix{T}, sgp4d::Sgp4Propagator{Tepoch, T}, Δt::Number, x₁::SVector{8, T}, y₁::SVector{7, T}; kwargs...) where {T<:Number, Tepoch<:Number} -> Nothing

Compute the SGP4 Jacobian by finite-differences using the propagator `sgp4d` at instant `Δt`
considering the input mean elements `x₁` that must provide the output vector `y₁`. The
result is written to the matrix `J`. Hence:

        ∂sgp4(x, Δt) │
    J = ──────────── │
             ∂x      │ x = x₁

# Keywords

- `perturbation::T`: Initial state perturbation to compute the finite-difference:
    `Δx = x * perturbation`. (**Default** = 1e-3)
- `perturbation_tol::T`: Tolerance to accept the perturbation. If the computed perturbation
    is lower than `perturbation_tol`, we increase it until it absolute value is higher than
    `perturbation_tol`. (**Default** = 1e-7)
"""
function _sgp4_jacobian!(
    J::AbstractMatrix{T},
    sgp4d::Sgp4Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{7, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7)
) where {T<:Number, Tepoch<:Number}

    num_states = 7
    dim_obs    = 6

    # Auxiliary variables.
    x₂ = x₁

    @inbounds for j in 1:num_states
        # State that will be perturbed.
        α = x₂[j]

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

        α += ϵ

        # Modify the perturbed state.
        x₂ = setindex(x₂, α, j)

        # Obtain the Jacobian by finite differentiation.
        _init_sgp4_with_state_vector!(sgp4d, x₂, sgp4d.epoch)
        r_teme, v_teme = sgp4!(sgp4d, Δt)
        y₂ = @SVector [
            r_teme[1],
            r_teme[2],
            r_teme[3],
            v_teme[1],
            v_teme[2],
            v_teme[3],
        ]

        J[:, j] .= (y₂ .- y₁) ./ ϵ

        # Restore the value of the perturbed state for the next iteration.
        x₂ = setindex(x₂, x₁[j], j)
    end

    return nothing
end
