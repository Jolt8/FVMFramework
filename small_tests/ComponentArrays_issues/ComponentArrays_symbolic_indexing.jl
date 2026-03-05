using ComponentArrays
using BenchmarkTools
using PreallocationTools
using Polyester

u = ComponentArray(a = [0.0, 0.0], b = [0.0, 0.0])

getindex(u, :a)
getproperty(u, :a)

u.a.= 0.0
u[:a][1] += 1.0
u.a[1]
#0.0 (should be 1.0)

u.a.= 0.0
u.a[1] += 1.0
u.a[1]
#1.0 (works fine)

#Strangely, symbolic indexing works when a single float is contained within the field

u = ComponentArray(a = 0.0, b = 0.0)

u.a = 0.0
u[:a] += 1.0
u.a 
#1.0 (works fine)


n_cells = 1000

u = ComponentArray(mass_fractions = (methanol = [0.0, 0.0], water = [0.0, 0.0]))

#symbolic indexing 
u.mass_fractions.methanol .= 0.0
u.mass_fractions[:methanol][1] += 1.0
u.mass_fractions.methanol[1]

function test_symbolic_indexing!(u)
    u.mass_fractions.methanol[1] = 0.0 
    u.mass_fractions[:methanol][1] += 1.0
    return 
end

@btime test_symbolic_indexing!($u) #476.842 ns (2 allocations: 7.87 KiB)
u.mass_fractions.methanol[1] #returns 0.0 (didn't work)

#dot indexing
u.mass_fractions.methanol .= 0.0
u.mass_fractions.methanol[1] += 1.0
u.mass_fractions.methanol[1]

function test_dot_indexing!(u)
    u.mass_fractions.methanol[1] = 0.0
    u.mass_fractions.methanol[1] += 1.0
    return 
end

@btime test_dot_indexing!($u) #2.900 ns, (0 allocations, 0 bytes) 
u.mass_fractions.methanol[1] #returns 1.0 (works)

#views 
u.mass_fractions.methanol .= 0.0
view(u.mass_fractions, :methanol)[1] += 1.0
u.mass_fractions.methanol[1]

function test_view_indexing!(u)
    u.mass_fractions.methanol[1] = 0.0 
    view(u.mass_fractions, :methanol)[1] += 1.0
end

@btime test_view_indexing!($u) #2.900 ns, (0 allocations, 0 bytes) 
u.mass_fractions.methanol[1] #returns 1.0 (works)

#getproperty indexing
u.mass_fractions.methanol .= 0.0
getproperty(u.mass_fractions, :methanol)[1] += 1.0
u.mass_fractions.methanol[1]

function test_getproperty_indexing!(u)
    u.mass_fractions.methanol[1] = 0.0 
    getproperty(u.mass_fractions, :methanol)[1] += 1.0
end

@btime test_getproperty_indexing!($u) #2.900 ns (0 allocations: 0 bytes)
u.mass_fractions.methanol[1] #returns 1.0 (works)

#Example without vectors contained in methanol and water
u = ComponentArray(mass_fractions = ComponentArray(methanol = 0.0, water = 0.0))

#symbolic indexing 
u.mass_fractions.methanol = 0.0
u.mass_fractions[:methanol] += 1.0
u.mass_fractions.methanol

function test_symbolic_indexing!(u)
    u.mass_fractions.methanol = 0.0 
    u.mass_fractions[:methanol] += 1.0
    return 
end

@btime test_symbolic_indexing!($u) #2.900 ns (0 allocations: 0 bytes)
u.mass_fractions.methanol #returns 1.0 (works)

#dot indexing
u.mass_fractions.methanol = 0.0
u.mass_fractions.methanol += 1.0
u.mass_fractions.methanol

function test_dot_indexing!(u)
    u.mass_fractions.methanol = 0.0
    u.mass_fractions.methanol += 1.0
    return 
end

@btime test_dot_indexing!($u) #2.900 ns, (0 allocations, 0 bytes) 
u.mass_fractions.methanol #returns 1.0 (works)


