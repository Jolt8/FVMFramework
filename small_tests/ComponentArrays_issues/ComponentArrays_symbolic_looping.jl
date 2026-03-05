using ComponentArrays
using BenchmarkTools
using PreallocationTools
using Polyester

n_cells = 100

#ComponentArrays
u_ca = ComponentArray(mass_fractions = (methanol = 0.0, water = 0.0))

u_ca.mass_fractions.methanol = 0.0
u_ca.mass_fractions.methanol += 1.0
u_ca.mass_fractions.methanol #returns 1.0 (works)

function test_symbolic_looping!(u)
    map(keys(u.mass_fractions)) do species_name
        view(u.mass_fractions, species_name)[1] = 0.0
        view(u.mass_fractions, species_name)[1] += 1.0
    end
    return 
end

@btime test_symbolic_looping!($u_ca) #886.364 ns (12 allocations: 448 bytes)
u_ca.mass_fractions #returns (methanol = 1.0, water = 1.0) (works)


#ComponentArrays
u_ca_nested_vector = ComponentArray(mass_fractions = ComponentArray(methanol = zeros(n_cells), water = zeros(n_cells)))

function test_symbolic_looping!(u, n_cells)
    #adding batch here also makes this much slower (1.769 ms (26936 allocations: 858.38 KiB) vs 801.200 μs (14934 allocations: 483.34 KiB))
    for cell_id in 1:n_cells
        map(keys(u.mass_fractions)) do species_name
            #doing u.mass_fractions[species_name][cell_id] = 0.0 doesn't updated the array stored inside the component array
            #it's also significantly slower
            u.mass_fractions[species_name][cell_id] = 0.0 
            u.mass_fractions[species_name][cell_id] += 1.0
        end
    end
    return 
end
#this may fail the first time with:
#ERROR: FieldError: type StrideArraysCore.AbstractPtrArray has no field `mass_fractions`, available fields: `ptr`, `sizes`, `strides`, `offsets`
#but running it again will work for some reason
@btime test_symbolic_looping!($u_ca_nested_vector, $n_cells) #801.200 μs (14934 allocations: 483.34 KiB)
u_ca_nested_vector.mass_fractions #returns (methanol = [0.0, 0.0...], water = [0.0, 0.0...]) (doesn't work)



#NamedTuples
u_nt = (mass_fractions = (methanol = [0.0], water = [0.0]), )

u_nt.mass_fractions.methanol[1] = 0.0
u_nt.mass_fractions.methanol[1] += 1.0
u_nt.mass_fractions.methanol[1] #returns 1.0 (works)

function test_symbolic_looping!(u)
    map(keys(u.mass_fractions)) do species_name
        u.mass_fractions[species_name][1] = 0.0
        u.mass_fractions[species_name][1] += 1.0
    end
    return 
end

@btime test_symbolic_looping!($u_nt) #2.900 ns (0 allocations: 0 bytes)
u_nt.mass_fractions #returns (methanol = [1.0], water = [1.0]) (works)


u_nt_nested_vector = (mass_fractions = (methanol = zeros(n_cells), water = zeros(n_cells)), )

function test_symbolic_looping!(u, n_cells)
    #adding @batch makes this much slower 685.500 μs (13936 allocations: 327.14 KiB) vs 613.376 ns (0 allocations: 0 bytes)
    for cell_id in 1:n_cells
        map(keys(u.mass_fractions)) do species_name
            u.mass_fractions[species_name][cell_id] = 0.0
            u.mass_fractions[species_name][cell_id] += 1.0
        end
    end
    return 
end

@btime test_symbolic_looping!($u_nt_nested_vector, $n_cells) #63.376 ns (0 allocations: 0 bytes)
u_nt_nested_vector.mass_fractions #returns (methanol = [1.0, 1.0...], water = [1.0, 1.0...]) (works)
