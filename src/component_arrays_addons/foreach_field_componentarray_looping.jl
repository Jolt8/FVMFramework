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

@generated function foreach_field_at!(f, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    local layout_keys
    if G1 <: VirtualFVMArray
        layout_keys = fieldnames(G1.parameters[2])
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        layout_keys = keys(Ax.parameters[1])
    else
        layout_keys = ()
    end

    exprs = []
    for n in layout_keys
        push!(exprs, :(f($(QuoteNode(n)), groups...)))
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end