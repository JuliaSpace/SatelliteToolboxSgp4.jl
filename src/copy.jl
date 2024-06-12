## Description #############################################################################
#
# Implement the copy for the structures related to the SGP4 propagator.
#
############################################################################################

# Define the copy function for the structure `Sgp4DeepSpace`.
let
    fields = fieldnames(Sgp4DeepSpace)
    expressions = [
        :(new_sgp4ds.$f = sgp4ds.$f)
        for f in fields
    ]

    @eval begin
        function Base.copy(sgp4ds::Sgp4DeepSpace{T}) where T<:Number
            new_sgp4ds = Sgp4DeepSpace{T}()
            $(expressions...)
            return new_sgp4ds
        end
    end
end

# Define the copy function for the structure `Sgp4Propagator`.
let
    fields = fieldnames(Sgp4Propagator)
    expressions = [
        :(new_sgp4d.$f = sgp4d.$f)
        for f in fields if f != :sgp4ds
    ]

    @eval begin
        function Base.copy(sgp4d::Sgp4Propagator{Tepoch, T}) where {Tepoch<:Number, T<:Number}
            new_sgp4d = Sgp4Propagator{Tepoch, T}()
            $(expressions...)
            # `sgp4ds` is the only mutable field in the structure.
            new_sgp4d.sgp4ds = copy(sgp4d.sgp4ds)
            return new_sgp4d
        end
    end
end
