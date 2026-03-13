using ComponentArrays
using BenchmarkTools
using Polyester

#the reason we can't use component arrays is because looping through nested fields doesn't work because it returns a copy instead of a view into the original vector
du = ComponentArray(a = [0.0, 0.0], b = [0.0, 0.0])

getindex(du, :a) #returns a copy of the original vector
getproperty(du, :a) #returns a view into the original vector

du.a.= 0.0
du[:a][1] += 1.0
du.a[1] #0.0 (should be 1.0)

du.a.= 0.0
du.a[1] += 1.0
du.a[1] #1.0 (works fine)

n_cells = 100

du_for_looping = ComponentArray(mass_fractions = (methanol = zeros(n_cells), water = zeros(n_cells)))

function test_symbolic_looping!(du, n_cells)
    #adding batch here also makes this much slower 184.900 μs (2802 allocations: 425.03 KiB) vs #98.700 μs (1600 allocations: 387.50 KiB)
    @batch for cell_id in 1:n_cells
        map(keys(du.mass_fractions)) do species_name
            #doing du.mass_fractions[species_name][cell_id] = 0.0 doesn't updated the array stored inside the component array
            #it's also significantly slower
            du.mass_fractions[species_name][cell_id] = 0.0 
            du.mass_fractions[species_name][cell_id] += 1.0
        end
    end
    return 
end
#this may fail the first time with:
#ERROR: FieldError: type StrideArraysCore.AbstractPtrArray has no field `mass_fractions`, available fields: `ptr`, `sizes`, `strides`, `offsets`
#but running it again will work for some reason
@btime test_symbolic_looping!($du_for_looping, $n_cells) #98.700 μs (1600 allocations: 387.50 KiB)
du_for_looping.mass_fractions #returns (methanol = [0.0, 0.0...], water = [0.0, 0.0...]) (doesn't work)