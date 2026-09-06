## Description #############################################################################
#
# Implement the copy for the structures related to the SGP4 propagator.
#
############################################################################################

# Define the copy function for the structure `Sgp4Propagator`. All fields are immutable,
# hence copying them is enough to obtain an independent structure.
let
    fields = fieldnames(Sgp4Propagator)
    expressions = [:(new_sgp4d.$f = sgp4d.$f) for f in fields]

    @eval begin
        function Base.copy(
            sgp4d::Sgp4Propagator{Tepoch, T}
        ) where {Tepoch <: Number, T <: Number}
            new_sgp4d = Sgp4Propagator{Tepoch, T}()
            $(expressions...)
            return new_sgp4d
        end
    end
end
