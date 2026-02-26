using BenchmarkTools

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

# An iterative recursive @inline approach without @generated
@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

function main()
    axes = (a=1:100, b=(c=101:200, d=201:300))
    v = rand(300)

    # Needs to be passed through a function barrier or const
    b1 = @benchmark create_views_gen($v, $axes)
    b2 = @benchmark create_views_inline($v, $axes)

    println("Using Generated Function:")
    display(b1)
    println()
    println("Using Inline Map:")
    display(b2)
    println()
end

main()

axes = (a=1:100, b=(c=101:200, d=201:300))
v = rand(300)

@btime create_views_gen(v, axes)
@btime create_views_inline(v, axes)

function test_views_gen(v, axes)
    named_tuple = create_views_gen(v, axes)

    named_tuple.a[1] = 1.0
    named_tuple.a[1] += 1.0
    for property_name in keys(named_tuple.b)
        for i in eachindex(named_tuple.b[property_name])
            named_tuple.b[property_name][i] += 1.0
        end
    end
end

function test_views_inline(v, axes)
    named_tuple = create_views_inline(v, axes)

    named_tuple.a[1] = 1.0
    named_tuple.a[1] += 1.0

    for property_name in keys(named_tuple.b)
        for i in eachindex(named_tuple.b[property_name])
            named_tuple.b[property_name][i] += 1.0
        end
    end
end

@btime test_views_gen($v, $axes)
@btime test_views_inline($v, $axes)
