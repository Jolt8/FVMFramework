using ComponentArrays

@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    local properties
    if G1 <: VirtualFVMArray
        properties = fieldnames(G1.parameters[2])
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        properties = keys(Ax.parameters[1])
    else
        properties = ()
    end

    exprs = []
    for n in properties
        view_exprs = []
        for i in 1:N
            if groups[i] <: VirtualFVMArray
                push!(view_exprs, :(_resolve_unified_val(groups[$i], Val{$(QuoteNode(n))}())))
            else
                push!(view_exprs, :(getproperty(groups[$i], $(QuoteNode(n)))))
            end
        end
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

@generated function for_fields!(f, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    exprs = []
    if G1 <: VirtualFVMArray
        # For VirtualFVMArray, axes are in the NamedTuple 'axes'
        AxType = G1.parameters[2] # NamedTuple type of axes
        for name in fieldnames(AxType)
            # Use getproperty at runtime to get the axis
            push!(exprs, :(f(getfield(getfield(groups[1], :axes), $(QuoteNode(name))), groups...)))
        end
    elseif G1 <: ComponentArray
        # For ComponentArray, the axes are in the 4th type parameter
        # AxTupleType: Tuple{Axis{...}}
        AxType = G1.parameters[4].parameters[1]
        # AxType: Axis{layout}
        layout = AxType.parameters[1] # This is the NamedTuple instance
        for i in eachindex(layout)
            axis_obj = layout[i]
            # Try to get the index/range from the type of the axis object
            # ComponentArrays uses the first type parameter for the index/range in most axes
            T = typeof(axis_obj)
            if hasproperty(T, :parameters) && length(T.parameters) >= 1
                idx = T.parameters[1]
                push!(exprs, :(f($idx, groups...)))
            else
                # Fallback
                push!(exprs, :(f($axis_obj, groups...)))
            end
        end
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

const foreach_field_at! = for_fields!