pipe_inside_diameter = 16.0u"mm"
pipe_length = 12.1u"inch" |> u"m"

stripped_pipe_length = ustrip(pipe_length |> u"m")
pipe_width = ustrip(pipe_inside_diameter |> u"m")

n_cells_axial = 100
n_layers = 5

grid_dimensions = (1, n_layers, n_cells_axial)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((pipe_width, pipe_width * n_layers, stripped_pipe_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

evaporator_endpoint = 3u"inch" |> u"m"
evaporator_endpoint_cell = Int(round(evaporator_endpoint / (pipe_length / n_cells_axial)))

addcellset!(grid, "pipe_inlet", xyz -> xyz[3] <= (1 * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "pipe_inlet")

addcellset!(grid, "silicon_carbide_preheater", xyz -> xyz[3] >= (1 * (stripped_pipe_length / n_cells_axial)) && xyz[3] <= (evaporator_endpoint_cell * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "silicon_carbide_preheater")

addcellset!(grid, "copper_mesh_reformer", xyz -> xyz[3] >= (evaporator_endpoint_cell * (stripped_pipe_length / n_cells_axial)) && xyz[3] <= ((n_cells_axial - 1) * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "copper_mesh_reformer")

addcellset!(grid, "pipe_outlet", xyz -> xyz[3] >= ((n_cells_axial - 1) * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "pipe_outlet")

addcellset!(grid, "steel_pipe_wall", xyz -> xyz[2] >= pipe_width && xyz[2] <= 2 * pipe_width)
getcellset(grid, "steel_pipe_wall")

addcellset!(grid, "thermocouple_and_jacket", xyz -> xyz[2] >= 2 * pipe_width && xyz[2] <= 3 * pipe_width)
getcellset(grid, "thermocouple_and_jacket")

addfacetset!(grid, "thermocouple_to_heating_wire", xyz -> xyz[2] == 3 * pipe_width)
getfacetset(grid, "thermocouple_to_heating_wire")

addcellset!(grid, "heating_wire", xyz -> xyz[2] >= 3 * pipe_width && xyz[2] <= 4 * pipe_width)
getcellset(grid, "heating_wire")

addcellset!(grid, "insulation", xyz -> xyz[2] >= 4 * pipe_width)
getcellset(grid, "insulation")

addfacetset!(grid, "insulation_to_air", xyz -> xyz[2] >= 5 * pipe_width)
getfacetset(grid, "insulation_to_air")

addfacetset!(grid, "pipe_endcaps_to_air", xyz -> xyz[2] >= pipe_width && xyz[2] <= 2 * pipe_width && (xyz[3] == 0.0 || xyz[3] == stripped_pipe_length))
getfacetset(grid, "pipe_endcaps_to_air")

n_cells = length(grid.cells)

u_proto = ComponentVector(
    mass_fractions = (
        water = zeros(n_cells)u"kg/kg",
        air = zeros(n_cells)u"kg/kg"
    ),
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)]

common_properties = get_5_layer_common_properties(pipe_length, cell_lengths_along_pipe, n_cells_axial)

pre_copper_mesh_reformer_properties = get_copper_mesh_reformer_properties(pipe_length, n_cells_axial, grid.cellsets["copper_mesh_reformer"], cell_lengths_along_pipe)
copper_mesh_reformer_properties = merge_properties(pre_copper_mesh_reformer_properties, common_properties)

pre_silicon_carbide_sand_properties = get_silicon_carbide_sand_properties()
silicon_carbide_sand_properties = merge_properties(pre_silicon_carbide_sand_properties, common_properties)

pre_steel_pipe_wall_properties = get_steel_pipe_wall_properties()
steel_pipe_wall_properties = merge_properties(pre_steel_pipe_wall_properties, common_properties)

pre_thermocouple_and_jacket_properties = get_thermocouple_and_jacket_properties(grid, pipe_length, n_cells_axial, cell_lengths_along_pipe)
thermocouple_and_jacket_properties = merge_properties(pre_thermocouple_and_jacket_properties, common_properties)

pre_heating_wire_properties = get_heating_wire_properties(grid, pipe_length, common_properties)
heating_wire_properties = merge_properties(pre_heating_wire_properties, common_properties)

pre_insulation_properties = get_insulation_properties()
insulation_properties = merge_properties(pre_insulation_properties, common_properties)

edit_pseudo_2D_geometry_5_layers!(config.geo, grid, common_properties)

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
reaction_names = keys(copper_mesh_reformer_properties.reactions.reforming_reactions)
species_names = keys(copper_mesh_reformer_properties.molecular_weights)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)",
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        insulation_to_air_overall_heat_transfer_coefficient_to_environment = u"W/(m^2*K)",
        pipe_endcaps_to_air_thermal_conductance = u"W/K",
        heater_weight_1 = u"1",
        heater_weight_2 = u"1",
        heater_weight_3 = u"1",
        heater_weight_4 = u"1",
        fluid_to_steel_pipe_convective_heat_transfer_coefficient = u"W/(m^2*K)",
        steel_thermal_mass_multiplier = u"1",
        thermocouple_to_heating_wire_thermal_resistance = u"K/W",
        wattage_received_per_m = u"W/m",
        molar_concentrations = u"mol/m^3",
        species_mass_flows = u"kg/s",
        net_rates = u"mol/s",
        mass = u"kg",
        mass_face = u"kg",
        mass_evaporated = u"kg",
        superficial_velocity = u"m/s",
        measured_heater_wattage_per_cell = u"W", 
        pipe_mass_flow = u"kg/s"
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",
        net_rates = (
            reforming_reactions = NamedTuple{reaction_names}(
                Tuple(zeros(n_cells)u"mol/s" for _ in 1:length(reaction_names))
            ),
        ), 
        molar_concentrations = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"mol/m^3" for _ in 1:length(species_names))
        ),
        species_mass_flows = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"kg" for _ in 1:length(species_names))
        )
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(
        insulation_to_air_overall_heat_transfer_coefficient_to_environment = 0.0u"W/(m^2*K)",
        pipe_endcaps_to_air_thermal_conductance = 0.0u"W/K",
        heater_weight_1 = 0.0u"1",
        heater_weight_2 = 0.0u"1",
        heater_weight_3 = 0.0u"1",
        heater_weight_4 = 0.0u"1",
        fluid_to_steel_pipe_convective_heat_transfer_coefficient = 0.0u"W/(m^2*K)",
        steel_thermal_mass_multiplier = 0.0,
        TC1_thermal_resistance = 0.0u"K/W",
        TC2_thermal_resistance = 0.0u"K/W",
        TC3_thermal_resistance = 0.0u"K/W",
        TC4_thermal_resistance = 0.0u"K/W",
        TC5_thermal_resistance = 0.0u"K/W"
    )
)

function fluid_physics_functions!(du, u, cell_id, vol)
end

function solid_physics_functions!(du, u, cell_id, vol)
end

placeholder_temperature = 99999.9u"K"
placeholder_mass_fractions = ComponentVector(
    water = 1e-99u"kg/kg",
    air = 1e-99u"kg/kg"
)

add_region!(
    config, "pipe_inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = silicon_carbide_sand_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]
        
        du.heat[cell_id] *= 0.0

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end
)

add_region!(
    config, "silicon_carbide_preheater";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = silicon_carbide_sand_properties,
    region_function =
    function evaporator!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "copper_mesh_reformer";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = copper_mesh_reformer_properties,
    region_function =
    function reactor!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "pipe_outlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = copper_mesh_reformer_properties,
    region_function =
    function outlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)
        
        du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]

        for_fields!(u.mass_fractions, du.species_mass_flows) do species, u_mass_fractions, du_species_mass_flows
            du_species_mass_flows[species[cell_id]] -= u.pipe_mass_flow[cell_id] * u_mass_fractions[species[cell_id]]
        end

        du.heat[cell_id] -= u.pipe_mass_flow[cell_id] * u.fluid_cp[cell_id] * u.temp[cell_id] 

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "steel_pipe_wall";
    type = Solid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = steel_pipe_wall_properties,
    region_function =
    function steel_pipe_wall!(du, u, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "thermocouple_and_jacket";
    type = ThermocoupleAndJacket(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = thermocouple_and_jacket_properties,
    region_function =
    function thermocouple_and_jacket!(du, u, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_patch!(
    config, "thermocouple_to_heating_wire";
    properties = ComponentVector(),
    patch_function =
    function thermocouple_to_heating_wire!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.temp[idx_b] - u.temp[idx_a]) / u.thermocouple_to_heating_wire_thermal_resistance[idx_a]
    end
)

add_region!(
    config, "heating_wire";
    type = HeatingWire(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = heating_wire_properties,
    region_function =
    function heating_wire!(du, u, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "insulation";
    type = Solid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = insulation_properties,
    region_function =
    function insulation!(du, u, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_patch!(
    config, "insulation_to_air";
    properties = ComponentVector(),
    patch_function =
    function insulation_to_air!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.insulation_to_air_overall_heat_transfer_coefficient_to_environment[idx_a] * u.insulation_to_environment_cell_areas[idx_a]) * (u.room_temperature[idx_a] - u.temp[idx_a])
    end
)

add_patch!(
    config, "pipe_endcaps_to_air";
    properties = ComponentVector(),
    patch_function =
    function pipe_endcaps_to_air!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.pipe_endcaps_to_air_thermal_conductance[idx_a]) * (u.room_temperature[idx_a] - u.temp[idx_a])
    end
)

function get_heated_cells(heater_start, heater_end)
    heated_cells = Int64[]
    for cell_id in collect(grid.cellsets["heating_wire"])
        distance_along_pipe = heating_wire_properties.cell_lengths_along_pipe[cell_id]
        if upreferred(heater_start) <= distance_along_pipe <= upreferred(heater_end)
            push!(heated_cells, cell_id)
        end
    end
    return heated_cells
end

total_heater_length = 10u"inch"
heater_start = 1u"inch"
n_heater_regions = 5

heater_1_start = heater_start + (total_heater_length / n_heater_regions) * 0
heater_1_end = heater_start + (total_heater_length / n_heater_regions) * 1

heater_2_start = heater_start + (total_heater_length / n_heater_regions) * 1
heater_2_end = heater_start + (total_heater_length / n_heater_regions) * 2

heater_3_start = heater_start + (total_heater_length / n_heater_regions) * 2
heater_3_end = heater_start + (total_heater_length / n_heater_regions) * 3

heater_4_start = heater_start + (total_heater_length / n_heater_regions) * 3
heater_4_end = heater_start + (total_heater_length / n_heater_regions) * 4

heater_5_start = heater_start + (total_heater_length / n_heater_regions) * 4
heater_5_end = heater_start + (total_heater_length / n_heater_regions) * 5

heater_1_cells = get_heated_cells(heater_1_start, heater_1_end)
n_heater_1_cells = length(heater_1_cells)

heater_2_cells = get_heated_cells(heater_2_start, heater_2_end)
n_heater_2_cells = length(heater_2_cells)

heater_3_cells = get_heated_cells(heater_3_start, heater_3_end)
n_heater_3_cells = length(heater_3_cells)

heater_4_cells = get_heated_cells(heater_4_start, heater_4_end)
n_heater_4_cells = length(heater_4_cells)

heater_5_cells = get_heated_cells(heater_5_start, heater_5_end)
n_heater_5_cells = length(heater_5_cells)

n_total_heated_cells = n_heater_1_cells + n_heater_2_cells + n_heater_3_cells + n_heater_4_cells + n_heater_5_cells

TC1_TC2_midpoint = (thermocouple_and_jacket_properties.TC1_length_along_pipe + thermocouple_and_jacket_properties.TC2_length_along_pipe) / 2
TC2_TC3_midpoint = (thermocouple_and_jacket_properties.TC2_length_along_pipe + thermocouple_and_jacket_properties.TC3_length_along_pipe) / 2
TC3_TC4_midpoint = (thermocouple_and_jacket_properties.TC3_length_along_pipe + thermocouple_and_jacket_properties.TC4_length_along_pipe) / 2
TC4_TC5_midpoint = (thermocouple_and_jacket_properties.TC4_length_along_pipe + thermocouple_and_jacket_properties.TC5_length_along_pipe) / 2

TC1_area_start = 0.0u"inch"
TC1_area_end = TC1_TC2_midpoint

TC2_area_start = TC1_TC2_midpoint
TC2_area_end = TC2_TC3_midpoint

TC3_area_start = TC2_TC3_midpoint
TC3_area_end = TC3_TC4_midpoint

TC4_area_start = TC3_TC4_midpoint
TC4_area_end = TC4_TC5_midpoint

TC5_area_start = TC4_TC5_midpoint
TC5_area_end = common_properties.pipe_length

function get_cells_along_pipe_range(common_properties, region_names, start_point, end_point)
    cells_along_pipe_range = Int64[]

    for region_name in region_names
        region_cell_ids = grid.cellsets[region_name]
        for cell_id in region_cell_ids
            if start_point <= common_properties.cell_lengths_along_pipe[cell_id] < end_point
                push!(cells_along_pipe_range, cell_id)
            end
        end
    end

    return cells_along_pipe_range
end

TC1_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC1_area_start, TC1_area_end)
TC2_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC2_area_start, TC2_area_end)
TC3_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC3_area_start, TC3_area_end)
TC4_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC4_area_start, TC4_area_end)
TC5_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC5_area_start, TC5_area_end)
