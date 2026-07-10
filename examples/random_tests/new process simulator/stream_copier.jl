
function _get_leaf_paths(layout::NamedTuple, current_path = ())
    paths = []
    for (k, v) in pairs(layout)
        if hasproperty(v, :ax) && typeof(v.ax) <: ComponentArrays.Axis
            inner_layout = typeof(v.ax).parameters[1]
            append!(paths, _get_leaf_paths(inner_layout, (current_path..., k)))
        else
            push!(paths, (current_path..., k))
        end
    end
    return paths
end

@generated function copy_stream!(streams::ComponentArray, src_id::Int, dest_id::Int)
    AxType = streams.parameters[4].parameters[1]
    layout = AxType.parameters[1]
    
    leaf_paths = _get_leaf_paths(layout)
    
    exprs = []
    for path in leaf_paths
        prop_expr = :streams
        for p in path
            prop_expr = :(getproperty($prop_expr, $(QuoteNode(p))))
        end
        push!(exprs, :($prop_expr[dest_id] = $prop_expr[src_id]))
    end
    
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end