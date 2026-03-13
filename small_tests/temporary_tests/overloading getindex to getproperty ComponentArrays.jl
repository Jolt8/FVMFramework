using ComponentArrays
using BenchmarkTools

u = ComponentVector(a = [0.0, 0.0], b = [0.0, 0.0])
Base.getindex(u::ComponentArray, i::Symbol) = Base._views(getproperty(u, i))
#@inline Base.setindex(x::ComponentArray, f::Symbol, v::Float64)::test = setproperty!(x::ComponentArray, f::Symbol, v::Float64)

getaxes(u)[1]

function test_get_index(x::ComponentVector, i::Symbol)
    @views getproperty(x, i)
end

function test_getindex_performance(ca) 
    for i in eachindex(ca.a)
        
        #ca[:a][i] #6 allocations
        ca.a[i] #0 allocations
        getproperty(ca, :a)[i] #0 allocations
        view(ca, :a)[i] #0 allocations
        #getindex(ca, :a)[i] #6 allocations
        test_get_index(ca, :a)[i] #6 allocations, 128 bytes

        #ca[:a][i] += 1.0 #8 allocations 
        @views ca[:a][i] += 1.0 #0 allocations
        ca.a[i] += 1.0 #0 allocations
        getproperty(ca, :a)[1] += 1.0 #0 allocations
        view(ca, :a)[i] += 1.0 #0 allocations
        #getindex(ca, :a)[i] += 1.0 #8 allocations
    end
    return 
end

test_getindex_performance(u)

@code_warntype test_getindex_performance(u)
@benchmark test_getindex_performance($u)
