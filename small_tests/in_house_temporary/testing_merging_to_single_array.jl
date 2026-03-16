using PreallocationTools
using ForwardDiff
using BenchmarkTools

#float only for comparison
test1 = zeros(100)

properties1 = ones(80)

cached1 = ones(20) 

function test_allocations(test, properties, cache)
    @views test[1:80] = properties

    @views test[81:100] = cache

    test[50] = test[1] + test[100]

    return nothing
end

@btime test_allocations($test1, $properties1, $cached1) #25.677 ns (0 allocations: 0 bytes)

#all duals
test = zeros(typeof(ForwardDiff.Dual(1.0, 1.0)), 100)

properties = ones(typeof(ForwardDiff.Dual(1.0, 1.0)), 80)

cached = ones(20) 
diff_cache = DiffCache(cached, 10)

function test_allocations(test, properties, cache)
    @views test[1:80] = properties

    @views test[81:100] = get_tmp(cache, ForwardDiff.Dual(1.0, 1.0))

    test[50] = test[1] + test[100]

    return nothing
end

@btime test_allocations($test, $properties, $diff_cache) #28.715 ns, 0 allocations


#half duals
test2 = zeros(typeof(ForwardDiff.Dual(1.0, 1.0)), 100)

properties2 = ones(80)

cached2 = ones(20) 
diff_cache2 = DiffCache(cached, 10)

function test_allocations2(test, properties, cache)
    @views test[1:80] = properties

    @views test[81:100] = get_tmp(cache, ForwardDiff.Dual(1.0, 1.0))

    test[50] = test[1] + test[100]

    return nothing
end

@btime test_allocations2($test2, $properties2, $diff_cache2) #19.559 ns, 0 allocations


#double cache
test3 = zeros(100)
test3_diff_cache = DiffCache(test3, 2)

properties3 = ones(80)

cached3 = ones(20) 
diff_cache3 = DiffCache(cached, 2)

function test_allocations4(test, properties, cache)
    test = get_tmp(test, ForwardDiff.Dual(1.0, 1.0))

    @views test[1:80] = properties

    @views test[81:100] = get_tmp(cache, ForwardDiff.Dual(1.0, 1.0))

    test[50] = test[1] + test[100]

    return nothing
end

@btime test_allocations4($test3_diff_cache, $properties3, $diff_cache3) #76.687 ns, 0 allocations


#triple cache
test4 = zeros(100)
test4_diff_cache = DiffCache(test3, 2)

properties4 = ones(80)
properties4_diff_cache = DiffCache(properties3, 1)

cached4 = ones(20) 
diff_cache4 = DiffCache(cached, 2)

function test_allocations4(test, properties, cache)
    test = get_tmp(test, ForwardDiff.Dual(1.0, 1.0))

    @views test[1:80] = get_tmp(properties, ForwardDiff.Dual(1.0, 1.0))

    @views test[81:100] = get_tmp(cache, ForwardDiff.Dual(1.0, 1.0))

    test[50] = test[1] + test[100]

    return nothing
end

@btime test_allocations4($test4_diff_cache, $properties4_diff_cache, $diff_cache4) #77.572 ns (0 allocations: 0 bytes)


function test_allocations4_floats_only(test, properties, cache)
    test = get_tmp(test, 1.0)

    @views test[1:80] = get_tmp(properties, 1.0)

    @views test[81:100] = get_tmp(cache, 1.0)

    test[50] = test[1] + test[100]

    return nothing
end

@btime test_allocations4_floats_only($test4_diff_cache, $properties4_diff_cache, $diff_cache4) #21.565 ns (0 allocations: 0 bytes)

#1. That's tricky, I think we're going to have to make the pre-allocated buffer also a DiffCache 
#(even for properties) which isn't ideal, but I don't think we should focus on that right 
#now because we can always make a separate type that holds differently sized Duals 
#(i.e same length for state and cached variables but properties would have a dual with only 1 number in it. 


