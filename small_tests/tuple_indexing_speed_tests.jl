initial_mass_fractions = (
    methanol=[1.0],
    water=[1.3],
    carbon_monoxide=[0.0001],
    hydrogen=[0.02],
    carbon_dioxide=[0.0001]
)

total_mass_fractions = 0.0

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    total_mass_fractions += initial_mass_fractions[species_name][1]
end

total_mass_fractions

function test_indexing_1(initial_mass_fractions)
    for (species_name, mass_fraction) in pairs(initial_mass_fractions)
        initial_mass_fractions[species_name][1] /= total_mass_fractions
    end
end

function test_indexing_2(initial_mass_fractions)
    for species_name in keys(initial_mass_fractions)
        initial_mass_fractions[species_name][1] /= total_mass_fractions
    end
end

function test_indexing_3(initial_mass_fractions)
    for (species_name, mass_fraction) in pairs(initial_mass_fractions)
        mass_fraction[1] /= total_mass_fractions
    end
end

function test_indexing_4(initial_mass_fractions)
    foreach(keys(initial_mass_fractions)) do key
        initial_mass_fractions[key][1] /= total_mass_fractions #10 allocations
        initial_mass_fractions[key][1] /= total_mass_fractions #10 allocations
    end
    initial_mass_fractions[:methanol] ./= 1.0 #0 allocation
    initial_mass_fractions[:methanol][1] /= 1.0 #0 allocation
    initial_mass_fractions[:methanol] ./= total_mass_fractions #1 allocation
    initial_mass_fractions[:methanol][1] /= total_mass_fractions #2 allocation
    #wtf, why does simply accessing total_mass_fractions cause allocations?
    initial_mass_fractions[:methanol] #0 allocations
    initial_mass_fractions[:methanol][1] #0 allocations
    getproperty(initial_mass_fractions, :methanol)[1] = total_mass_fractions #0 allocations
    initial_mass_fractions.methanol[1] = total_mass_fractions #0 allocations
    initial_mass_fractions[:methanol][1] = total_mass_fractions #0 allocations
end
#this is not good, every single access of the NamedTuple allocates 
#for example, the first two in the loop cause 20 allocations and the single property access after that causes 2 allocations

using BenchmarkTools

@btime test_indexing_1($initial_mass_fractions)
@btime test_indexing_2($initial_mass_fractions)
@btime test_indexing_3($initial_mass_fractions)
@btime test_indexing_4($initial_mass_fractions)

