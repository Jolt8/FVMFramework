using BenchmarkTools

function test_foreach_tuple(named_tuple)
    # The Julia compiler handles `map` and `foreach` over tuples very well
    # BUT, Keys itself is an iterator. We need to unroll using a Tuple of Symbols
    foreach(keys(named_tuple)) do property_name
        v = getproperty(named_tuple, property_name)
        for i in eachindex(v)
            v[i] += 1.0
        end
    end
end

function test_val_tuple(named_tuple)
    # Using Val to make the compiler unroll
    foreach(Val.(keys(named_tuple))) do val_name
        property_name = typeof(val_name).parameters[1]
        v = getproperty(named_tuple, property_name)
        for i in eachindex(v)
            v[i] += 1.0
        end
    end
end

function test_generated(named_tuple)
    _test_generated(named_tuple, Val(keys(named_tuple)))
end

@generated function _test_generated(named_tuple, ::Val{Names}) where Names
    exprs = []
    for name in Names
        push!(exprs, quote
            v = named_tuple.$name
            for i in eachindex(v)
                v[i] += 1.0
            end
        end)
    end
    return quote
        Base.@_inline_meta
        $(exprs...)
    end
end

v = rand(300)
axes = (a=1:100, b=(c=101:200, d=201:300))

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

named_tuple = create_views_inline(v, axes)

@btime test_foreach_tuple($named_tuple.b)
@btime test_val_tuple($named_tuple.b)
@btime test_generated($named_tuple.b)
