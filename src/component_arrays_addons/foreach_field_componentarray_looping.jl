#=@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{Any, N}) where {N}
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
end=#

#= Example:
    foreach_field_at!(cell_id, du.mass_fractions) do species, mass_fractions
        mass_fractions[species] = 1.0
    end
=#

@generated function foreach_field_at!(f, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    local properties
    if G1 <: VirtualFVMArray
        properties = fieldnames(G1.parameters[2])
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        properties = keys(Ax.parameters[1])
    else
        return :(nothing)
    end
    
    exprs = []
    if G1 <: VirtualFVMArray
        for n in properties
            push!(exprs, :(f(getfield(getfield(groups[1], :axes), $(QuoteNode(n))), groups...)))
        end
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        layout = Ax.parameters[1]
        offset = 0
        for (i, n) in enumerate(properties)
            layout_val = layout[i]
            if layout_val isa ViewAxis
                arr_length = length(layout_val.ax)
                args = [:(view(groups[$j], :)) for j in 1:N]
                push!(exprs, :(f($offset + 1:$(offset + arr_length), $(args...))))
                offset += arr_length
            else
                push!(exprs, :(f($layout_val, groups...)))
            end
        end
    end
    
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end
