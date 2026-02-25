using ComponentArrays
using BenchmarkTools

n_cells = 1000

initial_mass_fractions = ComponentVector(
    methanol = (1.0 * ones(n_cells)), 
    water = (1.3 * ones(n_cells)), 
    carbon_monoxide = (0.0001 * ones(n_cells)), 
    hydrogen = (0.02 * ones(n_cells)), 
    carbon_dioxide = (0.0001 * ones(n_cells))
)

total_mass_fractions = 2.3

for species_name in propertynames(initial_mass_fractions)
    for cell_id in 1:n_cells
        #this is another thing we have to keep in mind, we can't use the dot syntax on component when using propertynames
        getproperty(initial_mass_fractions, species_name)[cell_id] /= total_mass_fractions
    end
end

using BenchmarkTools

function no_pair_getproperty(initial_mass_fractions, total_mass_fractions, species)
    for species_name in species
        for cell_id in 1:n_cells
            getproperty(initial_mass_fractions, species_name)[cell_id] /= total_mass_fractions
        end
    end
end

function no_pair_view(initial_mass_fractions, total_mass_fractions, species)
    for species_name in species
        for cell_id in 1:n_cells
            view(initial_mass_fractions, species_name)[cell_id] /= total_mass_fractions
        end
    end
end 

function pair_getproperty(initial_mass_fractions, total_mass_fractions, species_tuple)
    for (species_name, value) in pairs(species_tuple)
        for cell_id in 1:n_cells
            getproperty(initial_mass_fractions, species_name)[cell_id] /= total_mass_fractions
        end
    end
end

function pair_view(initial_mass_fractions, total_mass_fractions, species_tuple)
    for (species_name, species_id) in pairs(species_tuple)
        for cell_id in 1:n_cells
            view(initial_mass_fractions, species_name)[cell_id] /= total_mass_fractions
        end
    end
end 

function pair_matrix(initial_mass_fractions, total_mass_fractions, species_tuple)
    for (species_name, species_id) in pairs(species_tuple)
        for cell_id in 1:n_cells
            initial_mass_fractions[species_id, cell_id] /= total_mass_fractions
        end
    end
end

function raw_matrix(initial_mass_fractions, total_mass_fractions)
    for cell_id in 1:n_cells
        initial_mass_fractions[1, cell_id] /= total_mass_fractions
        initial_mass_fractions[2, cell_id] /= total_mass_fractions
        initial_mass_fractions[3, cell_id] /= total_mass_fractions
        initial_mass_fractions[4, cell_id] /= total_mass_fractions
        initial_mass_fractions[5, cell_id] /= total_mass_fractions
    end
end

function raw_matrix_flipped(initial_mass_fractions, total_mass_fractions)
    for cell_id in 1:n_cells
        initial_mass_fractions[cell_id, 1] /= total_mass_fractions
        initial_mass_fractions[cell_id, 2] /= total_mass_fractions
        initial_mass_fractions[cell_id, 3] /= total_mass_fractions
        initial_mass_fractions[cell_id, 4] /= total_mass_fractions
        initial_mass_fractions[cell_id, 5] /= total_mass_fractions
    end
end

species_tuple = (methanol = 1, water = 2, carbon_monoxide = 3, hydrogen = 4, carbon_dioxide = 5)

initial_mass_fractions_component_vector = ComponentVector(
    methanol = (1.0 * ones(n_cells)), 
    water = (1.3 * ones(n_cells)), 
    carbon_monoxide = (0.0001 * ones(n_cells)), 
    hydrogen = (0.02 * ones(n_cells)), 
    carbon_dioxide = (0.0001 * ones(n_cells))
)

initial_mass_fractions_matrix = stack([
    1.0 * ones(n_cells), 
    1.3 * ones(n_cells), 
    0.0001 * ones(n_cells), 
    0.02 * ones(n_cells), 
    0.0001 * ones(n_cells)
], dims = 1)

initial_mass_fractions_matrix_2_dim = stack([
    1.0 * ones(n_cells), 
    1.3 * ones(n_cells), 
    0.0001 * ones(n_cells), 
    0.02 * ones(n_cells), 
    0.0001 * ones(n_cells)
], dims = 2)

total_mass_fractions = 2.3

@btime no_pair_getproperty($initial_mass_fractions_component_vector, $total_mass_fractions, $species) #42.8 us 305 alloc
@btime no_pair_view($initial_mass_fractions_component_vector, $total_mass_fractions, $species) #16.0 us 305 alloc

@btime pair_getproperty($initial_mass_fractions_component_vector, $total_mass_fractions, $species_tuple) #42.8 us 305 alloc
@btime pair_view($initial_mass_fractions_component_vector, $total_mass_fractions, $species_tuple) #16.0 us 305 alloc

@btime pair_matrix($initial_mass_fractions_matrix, $total_mass_fractions, $species_tuple) #7.225 us 205 alloc
@btime raw_matrix($initial_mass_fractions_matrix, $total_mass_fractions) #5.043 us 161 alloc

@btime raw_matrix_flipped($initial_mass_fractions_matrix_2_dim, $total_mass_fractions) #5.043 us 161 alloc

