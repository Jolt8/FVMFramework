using Polyester
using BenchmarkTools
using ComponentArrays
using LoopVectorization
using FVMFramework

test = ComponentVector(pressure = zeros(10000000))
cache = ComponentVector(temp = zeros(10000000), zoop = (zoop_1 = zeros(10000000), zoop_2 = zeros(10000000)))

axes = virtual_merge_axes((test, cache))

function test_function(test, cache, axes)
    @batch for i in eachindex(test)
        test.pressure[i] += 1.0
        cache.temp[i] += 1.0
    end
end

@btime test_function($test, $cache, $axes) #14.754 ms

test.pressure .= 0.0
cache.temp .= 0.0
function test_function_2(test, cache, axes)
    du = VirtualFVMArray((test, cache), axes)
    @batch for i in eachindex(test)
        test.pressure[i] += 1.0
        cache.temp[i] += 1.0
    end
end

@btime test_function_2($test, $cache, $axes) #14.933 ms

test .= 0.0
cache .= 0.0
function test_function_3(test, cache, axes)
    du = VirtualFVMArray((test, cache), axes)
    @batch for i in eachindex(test)
        du.pressure[i] += 1.0
        du.temp[i] += 1.0
    end
end

@btime test_function_3($test, $cache, $axes) 
VSCodeServer.@profview test_function_3(test, cache, axes)
#12.146 ms (1 allocation, 64 bytes) with @batch if the StrideArray has no .ptr field is not triggered 
    #- the way that error is triggred is by accessing a ComponentVector inside the @batch loop with something like du.pressure
    #- when this is triggered, it means that all other @batch loops are slowed down for the rest of the session
    #- until it's closed and reopened and the bug is not triggered
    #- I wonder if the small amount of allocations would be worth it for the speed boost
    #- it seems like it allocates once per @batch call with ((32 * n_top_level_fields_present_in_du) bytes)
    #- woah, for some reason setindexing a nested field twice brings the allocations up to 29.603 ms (7771 allocations: 419.17 KiB)
    #- therefore, it probably isn't worth it, but it's interesting to know that it's not just a simple 

#14.364 ms (0 allocations, 0 bytes) with @batch if the StrideArray has a .ptr field error is triggered
#12.031 ms (0 allocations, 0 bytes) with @tturbo (even if the StrideArray array error has been triggered)

test.pressure .= 0.0
cache.temp .= 0.0
function test_function_non_batched(test, cache, axes)
    for i in eachindex(test)
        test.pressure[i] += 1.0
        cache.temp[i] += 1.0
    end
end

@btime test_function_non_batched($test, $cache, $axes) #16.275 ms

test.pressure .= 0.0
cache.temp .= 0.0
function test_function_non_batched_2(test, cache, axes)
    du = VirtualFVMArray((test, cache), axes)
    for i in eachindex(test)
        test.pressure[i] += 1.0
        cache.temp[i] += 1.0
    end
end

@btime test_function_non_batched_2($test, $cache, $axes) #16.159 ms

test.pressure .= 0.0
cache.temp .= 0.0
function test_function_non_batched_3(test, cache, axes)
    du = VirtualFVMArray((test, cache), axes)
    for i in eachindex(test)
        du.pressure[i] += 1.0
        du.temp[i] += 1.0
    end
end

@btime test_function_non_batched_3($test, $cache, $axes) #18.780 ms

test_raw = zeros(10000000)
cache_raw = zeros(10000000)

function test_raw_array(test, cache)
    for i in eachindex(test)
        test[i] += 1.0
        cache[i] += 1.0
    end
end

@btime test_raw_array($test_raw, $cache_raw) #12.508 ms

test_raw .= 0.0
cache_raw .= 0.0
function test_raw_array_batched(test, cache)
    @batch for i in eachindex(test)
        test[i] += 1.0
        cache[i] += 1.0
    end
end

@btime test_raw_array_batched($test_raw, $cache_raw) 
#14.819 ms with @batch 
#11.988 ms with @tturbo