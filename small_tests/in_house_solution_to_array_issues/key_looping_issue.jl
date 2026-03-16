using Polyester
include("FVMArray.jl")

n_cells = 1000
proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
        ethanol = zeros(n_cells),
    ),
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
        ethanol = zeros(n_cells),
    ),
    pressure = zeros(n_cells),
)

ax, len = create_axes(proto)
data = zeros(len)
u = FVMArray(data, ax)

function loop_map!(u, n_cells)
    @batch for cell_id in 1:n_cells
        map(keys(u.mass_fractions)) do species
                u.mass_fractions[species][cell_id] += 1.0
                u.molar_concentrations[species][cell_id] += 1.0
        end
        nothing
    end
end

@btime loop_map!($u, $n_cells) #1.435ms, 46995 allocations


function loop_map_non_batch!(u, n_cells)
    for cell_id in 1:n_cells
        map(keys(u.mass_fractions)) do species
                u.mass_fractions[species][cell_id] += 1.0
                u.molar_concentrations[species][cell_id] += 1.0
        end
        nothing
    end
end

@btime loop_map_non_batch!($u, $n_cells) #6.520 μs (0 allocations: 0 bytes)


function test_symbolic_indexing_with_batch!(u, n_cells)
    @batch for cell_id in 1:n_cells
        #map(keys(u.mass_fractions)) do species
                u.mass_fractions[:methylene_blue][cell_id] += 1.0
                u.mass_fractions[:water][cell_id] += 1.0
                u.mass_fractions[:ethanol][cell_id] += 1.0
                u.molar_concentrations[:methylene_blue][cell_id] += 1.0
                u.molar_concentrations[:water][cell_id] += 1.0
                u.molar_concentrations[:ethanol][cell_id] += 1.0
        #end
        nothing
    end
end

@btime test_symbolic_indexing_with_batch!($u, $n_cells) #1.700 μs (1 allocation: 128 bytes)

#shrunk down
shrunk_n_cells = 100
shrunk_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
        ethanol = zeros(n_cells),
    ),
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
        ethanol = zeros(n_cells),
    ),
    pressure = zeros(n_cells),
)

shrunk_ax, shrunk_len = create_axes(proto)
shrunk_data = zeros(shrunk_len)
shrunk_u = FVMArray(shrunk_data, shrunk_ax)

function test_symbolic_indexing_with_batch_shrunk!(u, n_cells)
    @batch for cell_id in 1:n_cells
        #map(keys(u.mass_fractions)) do species
                u.mass_fractions[:methylene_blue][cell_id] += 1.0
                u.molar_concentrations[:methylene_blue][cell_id] += 1.0
        #end
        nothing
    end
end

@btime test_symbolic_indexing_with_batch_shrunk!($shrunk_u, $shrunk_n_cells) #587.222 ns (1 allocation: 128 bytes)