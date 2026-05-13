using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak

using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ForwardDiff
using DataFrames
using CSV
using DataInterpolations

using FVMFramework

#to see an explanation of the physical build of this reactor check out:
joinpath(@__DIR__, "..", "..", "notes", "reactor_design.md")

#plan for mapping layering to FVM groups:
    #cell group 1: fluid 
    #cell group 2: pipe wall (steel)
    #cell group 3: thermocouple with the two layers of mica tape and fiberglass insulation 
    #cell group 4: heating wire (this one should probably be lumped with the thermocouple group)
    #cell group 5: ceramic blanket + aluminum foil insulation

pipe_inside_diameter = 16.0u"mm"
pipe_length = 12.1u"inch" |> u"m"

stripped_pipe_length = ustrip(pipe_length |> u"m")
pipe_width = ustrip(pipe_inside_diameter |> u"m")

n_cells_axial = 100
n_layers = 4

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

addcellset!(grid, "thermocouple_and_heating_wire", xyz -> xyz[2] >= 2 * pipe_width && xyz[2] <= 3 * pipe_width)
getcellset(grid, "thermocouple_and_heating_wire")

addcellset!(grid, "insulation", xyz -> xyz[2] >= 3 * pipe_width)
getcellset(grid, "insulation")

addfacetset!(grid, "insulation_to_air", xyz -> xyz[2] >= 4 * pipe_width)
getfacetset(grid, "insulation_to_air")

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "..", "..", "physics", "multiphase.jl")) #the .. makes it go up one directory
Revise.includet(joinpath(@__DIR__, "..", "..", "physics", "energy.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "physics", "momentum.jl"))

#we use a polynomial due to the fact that the endcaps allow for a lot more heat to escape the ends of the reactor 
#than the cells closer to the center of the reactor 
#we don't need p1

Revise.includet(joinpath(@__DIR__, "packing_and_fluid_property_merging.jl"))

function update_copper_mesh_properties!(du, u, cell_id, vol, system)
    mw_avg!(u, cell_id)
    get_fluid_and_copper_mesh_packing_rho!(du, u, cell_id, vol)
    get_fluid_and_copper_mesh_packing_k!(du, u, cell_id, vol)
    get_fluid_and_copper_mesh_packing_cp!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
end

function update_silicon_carbide_bed_packing_properties!(du, u, cell_id, vol, system)
    mw_avg!(u, cell_id)
    get_fluid_and_silicon_carbide_bed_packing_rho!(du, u, cell_id, vol)
    get_fluid_and_silicon_carbide_bed_packing_k!(du, u, cell_id, vol)
    get_fluid_and_silicon_carbide_bed_packing_cp!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
end

function update_fluid_properties!(du, u, cell_id, vol, system)
    mw_avg!(u, cell_id)
    #rho_ideal!(u, cell_id)
    #rho_multiphase!(du, u, cell_id, vol)

    properties = ComponentVector(system.properties_vec, system.properties_axes)
    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]

    molar_concentrations!(u, cell_id)
    #update_velocity!(du, u, cell_id, vol)
end

function update_solid_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]
    #we still have to come up with a way to automatically do this
    #we could probably make a generated function that automatically applies fixed values to cached values 
    #I checked and this is working for now
end

#you could also add any of these functions individually

function fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    #it's definitely one of these
    sum_mass_flux_face_to_cell!(du, u, cell_id) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)

    #cap_evaporation_rate_to_phase_holdup!(du, u, cell_id, vol)
    #HAHA, I foud it, this function above causes the solver to dt_nan after the first time step!
end

function solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

u_proto = ComponentVector(
    mass_fractions = (
        methanol = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
        carbon_monoxide = zeros(n_cells)u"kg/kg",
        hydrogen = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        air = zeros(n_cells)u"kg/kg"
    ),
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
    #liquid_holdup = zeros(n_cells)u"m^3/m^3",
    #gas_holdup = zeros(n_cells)u"m^3/m^3"
)

config = create_fvm_config(grid, u_proto);

cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)]

Revise.includet(joinpath(@__DIR__, "..", "..", "properties", "common_properties.jl"))
common_properties = get_common_properties(pipe_length, cell_lengths_along_pipe, n_cells_axial)

Revise.includet(joinpath(@__DIR__, "..", "..", "properties", "packed_media_properties", "copper_mesh_reformer_properties.jl"))
pre_copper_mesh_reformer_properties = get_copper_mesh_reformer_properties(pipe_length, n_cells_axial, grid.cellsets["copper_mesh_reformer"], cell_lengths_along_pipe)
copper_mesh_reformer_properties = merge_properties(pre_copper_mesh_reformer_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "properties", "packed_media_properties", "silicon_carbide_sand_properties.jl"))
pre_silicon_carbide_sand_properties = get_silicon_carbide_sand_properties()
silicon_carbide_sand_properties = merge_properties(pre_silicon_carbide_sand_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "properties", "steel_pipe_wall_properties", "steel_pipe_wall_properties.jl"))
pre_steel_pipe_wall_properties = get_steel_pipe_wall_properties()
steel_pipe_wall_properties = merge_properties(pre_steel_pipe_wall_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "properties", "thermocouple_and_heating_wire_properties", "thermocouple_and_heating_wire_properties.jl"))
pre_thermocouple_and_heating_wire_properties = get_heating_wire_and_thermocouple_properties(grid, pipe_length, n_cells_axial, cell_lengths_along_pipe)
thermocouple_and_heating_wire_properties = merge_properties(pre_thermocouple_and_heating_wire_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "properties", "insulation_properties", "insulation_properties.jl"))
pre_insulation_properties = get_insulation_properties()
insulation_properties = merge_properties(pre_insulation_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "pseudo_2D_geometry_editing.jl"))
edit_pseudo_2D_geometry!(config.geo, grid, common_properties)

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
reaction_names = keys(copper_mesh_reformer_properties.reactions.reforming_reactions)
species_names = keys(copper_mesh_reformer_properties.molecular_weights)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)", #k and cp cause a dimsnion error for some reason
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        overall_heat_transfer_coefficient_to_environment = u"W/(m^2*K)",
        wattage_received_per_m = u"W/m",
        molar_concentrations = u"mol/m^3",
        species_mass_flows = u"kg/s",
        net_rates = u"mol/s",
        mass = u"kg",
        mass_face = u"kg",
        mass_evaporated = u"kg",
        superficial_velocity = u"m/s",
        #a correlation will be developed to account for the fact that the heating wire inside the reactor is unevenly spaced
        measured_heater_wattage_per_cell = u"W", 
        pipe_mass_flow = u"kg/s"
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",
        net_rates = (
            reforming_reactions = NamedTuple{reaction_names}(
                Tuple(zeros(n_cells)u"mol/s" for _ in 1:length(reaction_names))
            ), #don't forget the comma!
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
        overall_heat_transfer_coefficient_to_environment = 0.0u"W/(m^2*K)",
    )
)

function fluid_physics_functions!(du, u, cell_id, vol)
    #vaporization_model!(du, u, cell_id, vol)
    #ergun_momentum_friction!(du, u, cell_id, vol)
end

function solid_physics_functions!(du, u, cell_id, vol)
end

#TODO: remember to turn off the heater and to set 
#the inlet mass flow to what was observed in the experiment

placeholder_temperature = 99999.9u"K" #we want this to error hard if it doesn't get update by a trial's measured room temperature
placeholder_mass_fractions = ComponentVector(
    methanol = 1e-99u"kg/kg",
    water = 1e-99u"kg/kg",
    carbon_monoxide = 1e-99u"kg/kg",
    hydrogen = 1e-99u"kg/kg",
    carbon_dioxide = 1e-99u"kg/kg",
    air = 1e-99u"kg/kg"
)

add_region!(
    config, "pipe_inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
        #liquid_holdup = 1.0,
        #gas_holdup = 0.0
    ),
    properties = silicon_carbide_sand_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]
        
        du.heat[cell_id] *= 0.0
        #du.heat[cell_id] += u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id]

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
        #liquid_holdup = 0.0,
        #gas_holdup = 1.0
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
        #liquid_holdup = 0.0,
        #gas_holdup = 1.0
    ),
    properties = copper_mesh_reformer_properties,
    region_function =
    function reactor!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

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
        #liquid_holdup = 0.0,
        #gas_holdup = 1.0
    ),
    properties = copper_mesh_reformer_properties,
    region_function =
    function outlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)
        
        du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]

        for_fields!(u.mass_fractions, du.species_mass_flows) do species, u_mass_fractions, du_species_mass_flows
            du_species_mass_flows[species[cell_id]] -= u.pipe_mass_flow[cell_id] * u_mass_fractions[species[cell_id]]
        end
        #this is to prevent the concentration of all species from building up at the outlet

        #du.heat[cell_id] *= 0.0
        du.heat[cell_id] -= u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id] 
        #this is to prevent the temperature from building up at the outlet

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
        #liquid_holdup = 0.0,
        #gas_holdup = 1.0,
    ),
    properties = steel_pipe_wall_properties,
    region_function =
    function steel_pipe_wall!(du, u, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "thermocouple_and_heating_wire";
    type = Solid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
        #liquid_holdup = 0.0,
        #gas_holdup = 1.0,
    ),
    properties = thermocouple_and_heating_wire_properties,
    region_function =
    function thermocouple_and_heating_wire!(du, u, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        du.heat[cell_id] += u.measured_heater_wattage_per_cell[cell_id]

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
        #liquid_holdup = 0.0,
        #gas_holdup = 1.0,
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
    properties = ComponentVector(), #no new properties here
    patch_function =
    function insulation_to_air!(
        du, u,
        idx_a, idx_b, face_idx, #idx_b is not applicable here because it connects to nothing
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.overall_heat_transfer_coefficient_to_environment[idx_a] * u.insulation_to_environment_cell_areas[idx_a]) * (u.room_temperature[idx_a] - u.temp[idx_a])
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    #=new_pressure_driven_mass_flux!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    all_species_advection!(
        du, u, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
    
    enthalpy_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
    

    #=mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#
end

#we don't need a fluid_solid_flux! because temperature diffusion is the only thing happening here

function solid_solid_flux!(du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return solid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!
end

Revise.includet(joinpath(@__DIR__, "..", "packed_bed_reactor_dry", "empty_trial_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "packed_bed_reactor_dry", "air_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "thermocouple_data_processing", "thermocouple_data_with_wattage.jl"))


function dry_run_trial()
    dry_run_config = deepcopy(config)
    dry_run_properties = deepcopy(common_properties)

    empty_trial_properties = get_empty_trial_properties() 
    dry_run_properties = merge_properties(dry_run_properties, empty_trial_properties)

    air_properties = get_air_properties()
    dry_run_properties = merge_properties(dry_run_properties, air_properties)

    region_names = [dry_run_config.regions[i].name for i in eachindex(dry_run_config.regions)]
    regions = Dict(region_names .=> dry_run_config.regions)

    regions["pipe_inlet"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["pipe_inlet"].initial_conditions.mass_fractions = dry_run_properties.inlet_mass_fractions

    function dry_run_inlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end

    regions["pipe_inlet"].region_function = dry_run_inlet!

    regions["silicon_carbide_preheater"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["silicon_carbide_preheater"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions

    regions["copper_mesh_reformer"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["copper_mesh_reformer"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions

    regions["pipe_outlet"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["pipe_outlet"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions

    regions["steel_pipe_wall"].initial_conditions.temp = dry_run_properties.room_temperature

    regions["thermocouple_and_heating_wire"].initial_conditions.temp = dry_run_properties.room_temperature

    regions["insulation"].initial_conditions.temp = dry_run_properties.room_temperature

    thermocouple_data_path = joinpath(@__DIR__, "..", "thermocouple_data_processing", "heated_trial_tc_temps.csv")
    thermocouple_data = get_thermocouple_data(thermocouple_data_path)

    function measured_heater_wattage_per_cell(t)
        return thermocouple_data.heater_power_interp(ForwardDiff.value(t)) / n_heated_cells
    end

    return dry_run_properties, dry_run_config, thermocouple_data, measured_heater_wattage_per_cell
end

Revise.includet(joinpath(@__DIR__, "..", "packed_bed_reactor_water_flow_trial", "hot_water_trial_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "packed_bed_reactor_water_flow_trial", "water_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "packed_bed_reactor_water_flow_trial", "recorded_data_processing", "inlet_and_outlet_temperatures.jl"))
Revise.includet(joinpath(@__DIR__, "..", "packed_bed_reactor_water_flow_trial", "recorded_data_processing", "values_of_note.jl"))
Revise.includet(joinpath(@__DIR__, "..", "thermocouple_data_processing", "thermocouple_data.jl"))

function hot_water_trial()
    hot_water_config = deepcopy(config)
    hot_water_trial_properties = deepcopy(common_properties)

    values_of_note = get_packed_bed_reactor_water_flow_trial_values_of_note()
    pump_shutoff_timestamp = ustrip(upreferred(values_of_note.pump_shut_off_time))

    general_properties = get_hot_water_trial_properties(values_of_note)
    hot_water_trial_properties = merge_properties(hot_water_trial_properties, general_properties)

    water_properties = get_water_properties()
    hot_water_trial_properties = merge_properties(hot_water_trial_properties, water_properties)

    inlet_and_outlet_temperatures = get_inlet_and_outlet_temperature_correlations()
    inlet_temp_interp = inlet_and_outlet_temperatures.inlet_temp_interp

    region_names = [hot_water_config.regions[i].name for i in eachindex(hot_water_config.regions)]
    regions = Dict(region_names .=> hot_water_config.regions)

    regions["pipe_inlet"].initial_conditions.temperature = inlet_temp_interp(0.0u"s")
    regions["pipe_inlet"].initial_conditions.mass_fractions = hot_water_trial_properties.inlet_mass_fractions

    function hot_water_inlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]
        
        du.heat[cell_id] *= 0.0
        #du.heat[cell_id] += u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id]
        #the only reason we don't have to do this for the outlet is because the outlet actually does the du.heat[cell_id] += etc..

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end

    regions["pipe_inlet"].region_function = hot_water_inlet!

    regions["silicon_carbide_preheater"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["silicon_carbide_preheater"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions

    regions["copper_mesh_reformer"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["copper_mesh_reformer"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions

    regions["pipe_outlet"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["pipe_outlet"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions

    regions["steel_pipe_wall"].initial_conditions.temp = hot_water_trial_properties.room_temperature

    regions["thermocouple_and_heating_wire"].initial_conditions.temp = hot_water_trial_properties.room_temperature

    regions["insulation"].initial_conditions.temp = hot_water_trial_properties.room_temperature

    inlet_cell_id = collect(grid.cellsets["pipe_inlet"])[1]

    function pump_shut_off(du, u, cell_id, t)
        if t <= pump_shutoff_timestamp #pump on
            u.temp[inlet_cell_id] = inlet_temp_interp(ForwardDiff.value(t))

            u.pipe_mass_flow[cell_id] = ustrip(upreferred(hot_water_trial_properties.pipe_mass_flow))
        else #pump shut off
            u.pipe_mass_flow[cell_id] = 0.0
        end
    end

    thermocouple_data_path = joinpath(@__DIR__, "..", "thermocouple_data_processing", "hot_water_flow_tc_temps.csv")
    thermocouple_data = get_thermocouple_data(thermocouple_data_path)

    function measured_heater_wattage_per_cell(t)
        return thermocouple_data.heater_power_interp(ForwardDiff.value(t)) / n_heated_cells
    end
    
    return hot_water_config, hot_water_trial_properties, pump_shut_off, thermocouple_data, measured_heater_wattage_per_cell
end

#TODO: I think we could have an integrated way that allows for the properties of many trials to be sequentially added by doing something like in the initial conditions:
#=
    initial_conditions = (
        mass_fractions = TrialSelector(
            dry_run = empty_mass_fractions,
            hot_water = common_properties.inlet_mass_fractions,
        ),
        temp = TrialSelector(
            dry_run = common_properties.room_temperature,
            hot_water = inlet_temp_interp(0.0u"s")
        ),
        pressure = 1.0u"bar", #this is automatically applied to both trials
    )
=#
#we'll do that next time if it's necessary

fluid_regions = ["pipe_inlet", "silicon_carbide_preheater", "copper_mesh_reformer", "pipe_outlet"]
advecting_fluid_cells = vcat(collect(grid.cellsets["pipe_inlet"]), collect(grid.cellsets["silicon_carbide_preheater"]), collect(grid.cellsets["copper_mesh_reformer"]), collect(grid.cellsets["pipe_outlet"]))

heater_start = 1.0u"inch"
heater_end = 11.0u"inch"

heated_cells = findall(x -> heater_start <= x <= heater_end, [common_properties.cell_lengths_along_pipe[x] for x in collect(grid.cellsets["thermocouple_and_heating_wire"])])
n_heated_cells = length(heated_cells)

function trial_independent_solve_system!(du, u, p_vec, t, geo, system)
    #VERY IMPORTANT: since most software uses 0-based indexing, you need to adjust the cell id by +1
    #for example, if you mouse over cell_id 5161 in ParaView, you need to use 5162 in the code because julia uses 1-based indexing 

    #TODO: we probably should do something about the issue below
    #we could probably just add a pre_calculations function for this in each region group, but this works for now
    #we would have to decide whether or not it would be based on type like Fluid() vs Solid() or if each region would have its own function
    #I think region based is better just because being explicit is always great and a lot more flexible. 
    for reg in system.region_groups
        if reg.name == "pipe_inlet" || reg.name == "silicon_carbide_preheater"
            for cell_id in reg.region_cells
                update_silicon_carbide_bed_packing_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        elseif reg.name == "copper_mesh_reformer" || reg.name == "pipe_outlet"
            for cell_id in reg.region_cells
                update_copper_mesh_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        elseif reg.name in fluid_regions
            for cell_id in reg.region_cells
                update_fluid_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        else
            for cell_id in reg.region_cells
                update_solid_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        end
    end

    p = ComponentVector(p_vec, p_axes)

    for cell_id in eachindex(geo.cell_volumes)
        u.overall_heat_transfer_coefficient_to_environment[cell_id] = p.overall_heat_transfer_coefficient_to_environment[1]
    end

    for i in 1:length(advecting_fluid_cells) - 1 #we don't take mass out of the outlet
        idx_a = advecting_fluid_cells[i]
        idx_b = advecting_fluid_cells[i + 1]
        
        du.mass_face[idx_a, 6] -= u.pipe_mass_flow[idx_a]
        du.mass_face[idx_b, 1] += u.pipe_mass_flow[idx_a]
    end

    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system) #this seems to be the culprit
end

dry_run_properties, dry_run_config, dry_run_thermocouple_data, dry_run_measured_heater_wattage_per_cell = dry_run_trial();
dry_run_du0_vec, dry_run_u0_vec, dry_run_state_axes, dry_run_geo, dry_run_system = finish_fvm_config(dry_run_config, connection_map_function, check_units = false);

hot_water_properties, hot_water_config, hot_water_pump_shut_off, hot_water_thermocouple_data, hot_water_measured_heater_wattage_per_cell = hot_water_trial();
hot_water_du0_vec, hot_water_u0_vec, hot_water_state_axes, hot_water_geo, hot_water_system = finish_fvm_config(hot_water_config, connection_map_function, check_units = false);

p_axes = hot_water_system.p_axes

function dry_run_solve_system!(du, u, p_vec, t, geo, system)
    for cell_id in heated_cells
        u.measured_heater_wattage_per_cell[cell_id] = dry_run_measured_heater_wattage_per_cell(t) 
    end

    trial_independent_solve_system!(du, u, p_vec, t, geo, system)
end

function hot_water_solve_system!(du, u, p_vec, t, geo, system)
    for cell_id in advecting_fluid_cells
        hot_water_pump_shut_off(du, u, cell_id, t)
    end

    for cell_id in heated_cells
        u.measured_heater_wattage_per_cell[cell_id] = hot_water_measured_heater_wattage_per_cell(t) 
    end
    
    trial_independent_solve_system!(du, u, p_vec, t, geo, system)
end

f_closure_dry_run = (du, u, p, t) -> fvm_operator!(du, u, p, t, dry_run_solve_system!, dry_run_geo, dry_run_system)

f_closure_hot_water = (du, u, p, t) -> fvm_operator!(du, u, p, t, hot_water_solve_system!, hot_water_geo, hot_water_system)

#NOTE: p is not actually being used in the dry run or hot water functions right now because I haven't set them up yet
#I'm just using them as a placeholder for now
p_guess = ustrip.(Vector(ComponentVector(
    overall_heat_transfer_coefficient_to_environment = 9.612656359973201u"W/(m^2*K)", 
)))

function built_trial_implicit_prob(f_closure, du0_vec, u0_vec, thermocouple_data)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tMax = ustrip(upreferred(thermocouple_data.timestamps[end]))
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end

dry_run_implicit_prob = built_trial_implicit_prob(f_closure_dry_run, dry_run_du0_vec, dry_run_u0_vec, dry_run_thermocouple_data)
hot_water_implicit_prob = built_trial_implicit_prob(f_closure_hot_water, hot_water_du0_vec, hot_water_u0_vec, hot_water_thermocouple_data)

@time sol = solve(dry_run_implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)
@time sol = solve(dry_run_implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

#u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)]

#sim_file = @__FILE__

#root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\Julia_cfd_output_files"

#sol_to_vtk(sol, u_named, grid, sim_file, root_dir)

dry_run_timestamps = ustrip.(dry_run_thermocouple_data.timestamps)
hot_water_timestamps = ustrip.(hot_water_thermocouple_data.timestamps)

TC1_cell_id = Int(ustrip(thermocouple_and_heating_wire_properties.TC1_closest_cell_id))
TC2_cell_id = Int(ustrip(thermocouple_and_heating_wire_properties.TC2_closest_cell_id))
TC3_cell_id = Int(ustrip(thermocouple_and_heating_wire_properties.TC3_closest_cell_id))
TC4_cell_id = Int(ustrip(thermocouple_and_heating_wire_properties.TC4_closest_cell_id))
TC5_cell_id = Int(ustrip(thermocouple_and_heating_wire_properties.TC5_closest_cell_id))

n_trials = 2
n_saves = 100
dry_run_save_interval = ustrip(upreferred(dry_run_thermocouple_data.timestamps[end])) / n_saves
hot_water_save_interval = ustrip(upreferred(hot_water_thermocouple_data.timestamps[end])) / n_saves

function loss(θ)
    #Dry Run
    dry_run_loss_prob = remake(implicit_prob, p = θ)

    sol = solve(
        dry_run_loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        #InterpolatingAdjoint(autodiff = AutoMooncake()),
        #callback = approximate_time_to_finish_cb
        saveat = dry_run_save_interval
    )

    u_named = [ComponentVector(sol.u[i], state_axes) for i in eachindex(sol.u)]

    mean_squared_error = 0.0

    for i in eachindex(sol.t)
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC1_cell_id]) - ustrip(thermocouple_data.TC1_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC2_cell_id]) - ustrip(thermocouple_data.TC2_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC3_cell_id]) - ustrip(thermocouple_data.TC3_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC4_cell_id]) - ustrip(thermocouple_data.TC4_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC5_cell_id]) - ustrip(thermocouple_data.TC5_temps_interp(sol.t[i])))
        #mean_squared_error += abs2(ustrip(u_named[i].temp[end]) - ustrip(outlet_temp_interp(sol.t[i])))
        #I don't think using the measured outlet temperature is actually useful
    end

    #Hot Water 
    hot_water_loss_prob = remake(implicit_prob, p = θ)

    hot_water_sol = solve(
        hot_water_loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        #InterpolatingAdjoint(autodiff = AutoMooncake()),
        #callback = approximate_time_to_finish_cb
        saveat = hot_water_save_interval
    )

    #this is to prevent the optimizer from crashing the simulation to get a lower mean_squared_error
    if length(sol.t) + length(hot_water_sol.t) < n_trials * n_saves
        return 1e10
    end

    hot_water_u_named = [ComponentVector(hot_water_sol.u[i], state_axes) for i in eachindex(hot_water_sol.u)]

    for i in eachindex(hot_water_sol.t)
        mean_squared_error += abs2(ustrip(hot_water_u_named[i].temp[TC1_cell_id]) - ustrip(thermocouple_data.TC1_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(hot_water_u_named[i].temp[TC2_cell_id]) - ustrip(thermocouple_data.TC2_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(hot_water_u_named[i].temp[TC3_cell_id]) - ustrip(thermocouple_data.TC3_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(hot_water_u_named[i].temp[TC4_cell_id]) - ustrip(thermocouple_data.TC4_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(hot_water_u_named[i].temp[TC5_cell_id]) - ustrip(thermocouple_data.TC5_temps_interp(hot_water_sol.t[i])))
        #mean_squared_error += abs2(ustrip(u_named[i].temp[end]) - ustrip(outlet_temp_interp(sol.t[i])))
        #I don't think using the measured outlet temperature is actually useful
    end

    return mean_squared_error / (n_trials * n_saves)
end

p_guess_init = ComponentVector(
    overall_heat_transfer_coefficient_to_environment = 10.0u"W/(m^2*K)",
)

p_axes = getaxes(p_guess_init)
p_guess = ustrip.(upreferred.(Vector(p_guess_init)))

loss(p_guess)
#why is getting the gradient take so much longer than just running the simulation
#usually finding the gradient takes around 2x-5x longer than running the simulation once
#while the Tsit5() solver does obey this rule of thumb, any kind of FBDF method takes forever

@time grad = ForwardDiff.gradient(loss, p_guess)
#1.967 s (1368576 allocations: 369.83 MiB)
#that's great, it means that these types of 1D problems are not very computationally intensive for finding gradients

#jac = ForwardDiff.jacobian(loss, p_guess)


#OPTIMIZATION

#SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true #turn to true to debug EnzymeJVP

#Logging.disable_logging(Logging.Warn)  # Disable all warnings
#Logging.disable_logging(Logging.Warn - 1)  # enable all warnings

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

p_lower_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    overall_heat_transfer_coefficient_to_environment = 0.001u"W/(m^2*K)",
))))

p_upper_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    overall_heat_transfer_coefficient_to_environment = 1000.0u"W/(m^2*K)", 
))))

optprob = Optimization.OptimizationProblem(optf, p_guess, lb=p_lower_bounds, ub=p_upper_bounds)

function randomize(lower, upper)
    return lower + (upper - lower) * rand()
end

#p_ensemble = [[randomize(p_lower_bounds[i], p_upper_bounds[i]) for i in eachindex(p_lower_bounds)] for _ in 1:Sys.CPU_THREADS]

p_ensemble = exp.(range(log(p_lower_bounds[1]), log(p_upper_bounds[1]), length = 100))
#length=Sys.CPU_THREADS))

for i in eachindex(p_ensemble)
    println(i)
    loss_i = loss([p_ensemble[i]])
    println(loss_i)
    push!(resulting_losses, loss_i)
end

#resulting_losses = [loss([p_ensemble[i]]) for i in eachindex(p_ensemble)]

plot(p_ensemble[1:80], resulting_losses[101:180])

p_ensemble[67]
resulting_losses[167]

CSV.write("optimization_results_2026-05-11_09-01-04.csv", DataFrame(loss = resulting_losses[101:end], overall_heat_transfer_coefficient_to_environment = p_ensemble))

#a heat transfer coefficient of 10.0  seems to give the lowest loss

function prob_func(prob, i, repeat)
    return remake(prob, p = p_ensemble[i])  
end

ensembleprob = EnsembleProblem(optprob; prob_func)

LOSS = Float64[]
PARS = []

using Dates

# Ensure the directory exists and use a filename-safe date format (colons are invalid on Windows)
mkpath(joinpath(@__DIR__, "optimization_results"))
results_path = joinpath(@__DIR__, "optimization_results", "optimization_results_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).csv")

# Create the file and manually write the header string using propertynames
open(results_path, "w") do io
    header_str = "loss," * join(string.(propertynames(p_guess_init)), ",")
    println(io, header_str)
end

const cb_lock = ReentrantLock()

cb = function (state, l)
    display(l)
    display(state.u)
    
    lock(cb_lock) do
        push!(LOSS, l)
        push!(PARS, state.u)
        
        # Convert state.u to a named tuple using your p_axes so the CSV has nice column headers
        p_named = NamedTuple(ComponentVector(state.u, p_axes))
        row = merge((loss = l, ), p_named)
        
        # CSV.write with append=true automatically opens, appends, and closes (flushes) the file
        CSV.write(results_path, DataFrame([row]), append=true)
    end
    
    false
end

#=
@time res = Optimization.solve(
    ensembleprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    EnsembleThreads(),
    trajectories = 50,
    #Sys.CPU_THREADS,
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    #f_abstol=1e-8,
    #g_abstol=1e-8,
)
=#

@time res = Optimization.solve(
    optprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    #EnsembleThreads(),
    #trajectories = 50,
    #Sys.CPU_THREADS,
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    #f_abstol=1e-8,
    #g_abstol=1e-8,
)

#=
using OptimizationBBO #for BlackBoxOptim
@time res = solve(
    ensembleprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    EnsembleThreads(),
    trajectories = Sys.CPU_THREADS,
    callback = cb,
)
=#

res.u0_vec

loss(ustrip.(upreferred.(Vector(p_vec_fitted))))

p_fitted = ComponentVector(res.u, p_axes)

#0.41741
#0.045
#7050
#8.71

results_path = joinpath(@__DIR__, "optimization_results", "optimization_results_2026-05-10_14-09-03.csv")

results_data = CSV.read(results_path, DataFrame)

results_data.overall_heat_transfer_coefficient_to_environment

middle_thermocouple_temps = []

for i in eachindex(results_data.overall_heat_transfer_coefficient_to_environment[1:20])
    p_test = [results_data.overall_heat_transfer_coefficient_to_environment[i]]

    test_prob = remake(implicit_prob, p = p_test)
    
    sol = solve(
        test_prob,
        FBDF(linsolve = SparspakFactorization(),), 
    )

    u_named = [ComponentVector(sol.u[i], state_axes) for i in eachindex(sol.u)]

    push!(middle_thermocouple_temps, [u_named[j].temp[TC3_cell_id] for j in eachindex(u_named)])
end


