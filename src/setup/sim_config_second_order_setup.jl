function common_merging_of_second_order_vars!(array_found_in, merged_caches, field, n_cells)
    grad_sym = Symbol(String(field) * "_grad")
    min_sym = Symbol(String(field) * "_min")
    max_sym = Symbol(String(field) * "_max")

    if getproperty(array_found_in, field) isa ComponentVector
        names = keys(getproperty(array_found_in, field))
        grad_vec = ComponentVector(NamedTuple{names}(
            Tuple(zeros(n_cells, 3) .* unit(array_found_in[field][1]) .* 0.0 ./ 1.0u"m" for _ in 1:length(names))
        ))
        max_vec = ComponentVector(NamedTuple{names}(
            Tuple(zeros(n_cells) .* unit(array_found_in[field][1]) .* 0.0 for _ in 1:length(names))
        ))
        min_vec = ComponentVector(NamedTuple{names}( 
            Tuple(zeros(n_cells) .* unit(array_found_in[field][1]) .* 0.0 for _ in 1:length(names))
        ))
        merged_caches = merge_properties(merged_caches, ComponentVector(
            NamedTuple{(grad_sym, min_sym, max_sym)}((grad_vec, min_vec, max_vec)))
        )
    else
        merged_caches = merge_properties(merged_caches, ComponentVector(
            NamedTuple{(grad_sym, min_sym, max_sym)}((zeros(n_cells, 3), zeros(n_cells), zeros(n_cells))))
        )
    end

    return merged_caches
end

function merge_in_second_order_caches(config, merged_caches, merged_properties)
    n_cells = length(config.geo.cell_volumes)
    second_order_syms = config.second_order_syms

    for field in propertynames(merged_properties)
        if field in second_order_syms
            merged_caches = common_merging_of_second_order_vars!(merged_properties, merged_caches, field, n_cells)
        end
    end

    for field in propertynames(merged_caches)
        if field in second_order_syms
            merged_caches = common_merging_of_second_order_vars!(merged_caches, merged_caches, field, n_cells)
        end
    end

    for field in propertynames(config.u_proto)
        if field in second_order_syms
            merged_caches = common_merging_of_second_order_vars!(config.u_proto, merged_caches, field, n_cells)
        end
    end

    return merged_caches
end