## Description #############################################################################
#
# Implement the copy for the structures related to the SGP4 propagator.
#
############################################################################################

# Define the copy function for the structure `Sgp4DeepSpace`.
@generated function Base.copy(sgp4ds::Sgp4DeepSpace{T}) where T<:Number
    fields = fieldnames(sgp4ds)

    expressions = Expr[:(new_sgp4ds = Sgp4DeepSpace{$T}())]
    sizehint!(expressions, length(fields) + 1)

    @inbounds for f in fields
        push!(expressions, :(new_sgp4ds.$f = sgp4ds.$f))
    end

    return :(
        $(expressions...);
        return new_sgp4ds
    )
end

# Define the copy function for the structure `Sgp4Propagator`.
@generated function Base.copy(
    sgp4d::Sgp4Propagator{Tepoch, T}
) where {Tepoch<:Number, T<:Number}
    fields = fieldnames(sgp4d)

    expressions = Expr[:(new_sgp4d = Sgp4Propagator{$Tepoch, $T}())]
    sizehint!(expressions, length(fields) + 1)

    @inbounds for f in fields
        if f != :sgp4ds
            push!(expressions, :(new_sgp4d.$f = sgp4d.$f))
        else
            # `sgp4ds` is the only mutable element inside this structure.
            push!(expressions, :(new_sgp4d.$f = copy(sgp4d.$f)))
        end
    end

    return :(
        $(expressions...);
        return new_sgp4d
    )
end
