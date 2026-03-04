## Description #############################################################################
#
# Types and stub functions for ML-∂SGP4 (neural-network-corrected SGP4).
# Actual implementations are provided by the SatelliteToolboxSgp4LuxExt extension,
# which is loaded automatically when Lux.jl, Optimisers.jl, and Zygote.jl are available.
#
# Reference:
#   Acciarini, G., Baydin, A. G., & Izzo, D. (2025). "Closing the gap between SGP4 and
#   high-precision propagation via differentiable programming." Acta Astronautica, 226,
#   694-701.  https://doi.org/10.1016/j.actaastro.2024.10.063
#
############################################################################################

export MLdSGP4Config
export ml_dsgp4_init, ml_dsgp4, ml_dsgp4!, ml_dsgp4_train, ml_dsgp4_save, ml_dsgp4_load
export freeze

# ------------------------------------------------------------------------------------------
#                                     MLdSGP4Config
# ------------------------------------------------------------------------------------------

"""
    MLdSGP4Config

Architecture and normalization hyperparameters for ML-∂SGP4.

Default values match the ESA dSGP4 reference implementation
([mldsgp4.py](https://github.com/esa/dSGP4/blob/master/dsgp4/mldsgp4.py)).

# Fields

- `hidden_layers::Vector{Int}`: Neuron counts for each hidden layer.  Length determines the
  number of hidden layers; each element sets the width of that layer.
- `normalization_R::Float64`: Position normalization constant [km].
- `normalization_V::Float64`: Velocity normalization constant [km/s].
- `input_correction::Float64`: Initial scale for input correction parameters (α).
- `output_correction::Float64`: Initial scale for output correction parameters (β).

# Examples

```julia
MLdSGP4Config()                              # paper defaults: [100, 100]
MLdSGP4Config(hidden_layers = [64, 64])      # smaller for fast experiments
MLdSGP4Config(hidden_layers = [128, 64, 32]) # 3 hidden layers, tapering
```
"""
struct MLdSGP4Config{RT<:Number, VT<:Number, IT<:Number, OT<:Number}
    hidden_layers::Vector{Int}
    normalization_R::RT
    normalization_V::VT
    input_correction::IT
    output_correction::OT
end

function MLdSGP4Config(;
    hidden_layers     = [100, 100],
    normalization_R   = 6958.137,
    normalization_V   = 7.947155867983262,
    input_correction  = 1e-2,
    output_correction = 0.8,
)
    MLdSGP4Config(
        hidden_layers, normalization_R, normalization_V,
        input_correction, output_correction,
    )
end

# ------------------------------------------------------------------------------------------
#                                   Stub Functions
# ------------------------------------------------------------------------------------------

const _LUX_EXT_ERR = "This function requires Lux.jl, Optimisers.jl, and Zygote.jl. " *
                      "Load them with `using Lux, Optimisers, Zygote`."

function ml_dsgp4_init end
function ml_dsgp4 end
function ml_dsgp4! end
function ml_dsgp4_train end
function ml_dsgp4_save end
function ml_dsgp4_load end
function freeze end
