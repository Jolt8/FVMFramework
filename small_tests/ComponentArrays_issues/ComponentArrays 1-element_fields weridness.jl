using ComponentArrays
using BenchmarkTools
using PreallocationTools
using Polyester

u = ComponentArray(mass_fractions = ComponentArray(methanol = 0.0, water = 0.0))

#symbolic indexing (doesn't work in @batch)
u.mass_fractions.methanol = 0.0
u.mass_fractions[:methanol] += 1.0
u.mass_fractions.methanol[1] #returns 1.0, works

function test_symbolic_indexing!(u)
    u.mass_fractions.methanol = 0.0 
    u.mass_fractions[:methanol] += 1.0
    return 
end

@btime test_symbolic_indexing!($u) #2.600 ns (0 allocations: 0 bytes) #faster than indexing the returned vector
u.mass_fractions.methanol[1] #returns 1.0 (works)

#dot indexing
u.mass_fractions.methanol = 0.0
u.mass_fractions.methanol += 1.0
u.mass_fractions.methanol #works

function test_dot_indexing!(u)
    u.mass_fractions.methanol = 0.0
    u.mass_fractions.methanol += 1.0
    return 
end

@btime test_dot_indexing!($u) #2.600 ns (0 allocations: 0 bytes)
u.mass_fractions.methanol[1] #returns 1.0 (works)


