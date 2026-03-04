using BenchmarkTools
using Polyester

function test_alloc(v)
    m = eachcol(reshape(v[1], 4, 100))
    @batch for i in 1:100
        for j in 1:4
            m[i][j] += 1.0
        end
    end
    
end

v = ones(400)

v_test = (v,)

@btime test_alloc($v)
@btime test_alloc($v_test)
