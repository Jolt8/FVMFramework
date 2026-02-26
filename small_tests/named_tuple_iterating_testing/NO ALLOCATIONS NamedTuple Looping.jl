using BenchmarkTools

# An iterative recursive @inline approach without @generated
@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

function test_map_keys_inline(v1, v2, v3, axes)
    named_tuple_1 = create_views_inline(v1, axes)
    named_tuple_2 = create_views_inline(v2, axes)
    named_tuple_3 = create_views_inline(v3, axes)

    map(keys(named_tuple_1), values(named_tuple_1)) do name_sym, species_arr
        other_arr = getproperty(named_tuple_2, name_sym)
        third_arr = getproperty(named_tuple_3, name_sym)

        for i in eachindex(species_arr)
            species_arr[i] += other_arr[i] + third_arr[i]
        end
    end

    return named_tuple_1, named_tuple_2, named_tuple_3
end

function test_map_keys_inline_simple(v1, v2, v3, axes, cell_id)
    named_tuple_1 = create_views_inline(v1, axes)
    named_tuple_2 = create_views_inline(v2, axes)
    named_tuple_3 = create_views_inline(v3, axes)

    map(keys(named_tuple_1)) do name_sym
        named_tuple_1[name_sym][cell_id] += 1.0
        named_tuple_2[name_sym][cell_id] += 1.0
        named_tuple_3[name_sym][cell_id] += 1.0
    end

    return named_tuple_1, named_tuple_2, named_tuple_3
end

v1 = rand(100);
v2 = rand(100);
v3 = rand(100);
axes = (c=1:50, d=51:100)

function test_looped_inline(v1, v2, v3, axes, cell_ids)
    named_tuple_1 = create_views_inline(v1, axes)
    named_tuple_2 = create_views_inline(v2, axes)
    named_tuple_3 = create_views_inline(v3, axes)

    for cell_id in cell_ids
        for name_sym in keys(named_tuple_1)
            named_tuple_1[name_sym][cell_id] += 1.0
            named_tuple_2[name_sym][cell_id] += 1.0
            named_tuple_3[name_sym][cell_id] += 1.0
        end
    end

    return named_tuple_1, named_tuple_2, named_tuple_3
end

@benchmark test_looped_inline(v1, v2, v3, axes, 1)

@generated function create_views_gen(v::AbstractVector, axes::NamedTuple{Names,Types}) where {Names,Types}
    exprs = []
    for (name, T) in zip(Names, Types.parameters)
        if T <: NamedTuple
            push!(exprs, :($(name) = create_views_gen(v, axes.$(name))))
        else
            push!(exprs, :($(name) = view(v, axes.$(name))))
        end
    end
    return :((; $(exprs...)))
end

function test_map_keys_gen(v1, v2, v3, axes)
    named_tuple_1 = create_views_gen(v1, axes)
    named_tuple_2 = create_views_gen(v2, axes)
    named_tuple_3 = create_views_gen(v3, axes)

    map(keys(named_tuple_1), values(named_tuple_1)) do name_sym, species_arr
        other_arr = getproperty(named_tuple_2, name_sym)
        third_arr = getproperty(named_tuple_3, name_sym)

        for i in eachindex(species_arr)
            species_arr[i] += other_arr[i] + third_arr[i]
        end
    end

    return named_tuple_1, named_tuple_2, named_tuple_3
end

@benchmark test_map_keys_gen(v1, v2, v3, axes)

