using BenchmarkTools
using Polyester

@inline create_views_inline(v, ax::Tuple{UnitRange{Int64},NTuple{N,Int64}}) where N = eachcol(reshape(view(v, ax[1]), ax[2]))

function test_alloc(v, ax)
    m = create_views_inline(v, ax)
    s = 0.0
    for i in 1:100
        for j in 1:4
            s += m[i][j]
        end
    end
    s
end

function test_alloc_batch(v, ax)
    m = create_views_inline(v, ax)
    s = zeros(100)
    @batch for i in 1:100
        for j in 1:4
            s[i] += m[i][j]
        end
    end
    sum(s)
end

v = ones(400)
ax = (1:400, (4, 100))
println(test_alloc(v, ax))
println(test_alloc_batch(v, ax))

@btime test_alloc($v, $ax)
@btime test_alloc_batch($v, $ax)
