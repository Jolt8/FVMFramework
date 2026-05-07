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
using Mooncake

using FVMFramework

pipe_inside_diameter = 0.5u"inch" |> u"m"
pipe_length = 12.1u"inch" |> u"m"

stripped_pipe_length = ustrip(pipe_length |> u"m")
pipe_width = ustrip(pipe_inside_diameter |> u"m")

n_cells = 100

grid_dimensions = (1, 1, n_cells)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((pipe_width, pipe_width, stripped_pipe_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "inlet", xyz -> xyz[3] <= (1 * (stripped_pipe_length / n_cells)))
getcellset(grid, "inlet")

addcellset!(grid, "evaporator", xyz -> xyz[3] >= (1 * (stripped_pipe_length / n_cells)) && xyz[3] <= (19 * (stripped_pipe_length / n_cells)))
getcellset(grid, "evaporator")

addcellset!(grid, "reactor", xyz -> xyz[3] >= (19 * (stripped_pipe_length / n_cells)) && xyz[3] <= (99 * (stripped_pipe_length / n_cells)))
getcellset(grid, "reactor")

addcellset!(grid, "outlet", xyz -> xyz[3] >= (99 * (stripped_pipe_length / n_cells)))
getcellset(grid, "outlet")

struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "..", "..", "physics", "multiphase.jl")) #the .. makes it go up one directory
Revise.includet(joinpath(@__DIR__, "..", "..", "physics", "energy.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "physics", "momentum.jl"))


#this is required because since cp is a cache to allow duals such as u.empty_reactor_thermal_mass to pass through the system
function update_cp!(du, u, cell_id, vol)
    # J/K = (J/kg/K * kg/m^3 * m^3) + (J/(m*K) * m)
    u.cell_combined_stationary_thermal_mass[cell_id] = (u.cp[cell_id] * u.rho[cell_id] * vol) + (u.empty_reactor_thermal_mass[1] * u.per_cell_pipe_length[cell_id])
end

#we use a polynomial due to the fact that the endcaps allow for a lot more heat to escape the ends of the reactor 
#than the cells closer to the center of the reactor 
#we don't need p1 
function update_UA_to_environment!(du, u, cell_id, vol)
    u.UA_to_environment_per_m[cell_id] = u.UA_to_environment_per_m_p2 * (u.per_cell_pipe_length[cell_id] - (u.pipe_length[1] / 2))^2 + u.UA_to_environment_per_m_p0
end

function update_properties!(du, u, cell_id, vol)
    update_cp!(du, u, cell_id, vol)
    mw_avg!(u, cell_id)
    rho_ideal!(u, cell_id)
    #rho_multiphase!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
    update_velocity!(du, u, cell_id, vol)
end

#you could also add any of these functions individually

function sum_and_cap_fluxes!(du, u, cell_id, vol)
    sum_mass_flux_face_to_cell!(du, u, cell_id) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)

    cap_evaporation_rate_to_phase_holdup!(du, u, cell_id, vol)
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
    liquid_holdup = zeros(n_cells)u"m^3/m^3",
    gas_holdup = zeros(n_cells)u"m^3/m^3",
    TC_temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

Revise.includet(joinpath(@__DIR__, "..", "..", "common_reformer_properties", "common_reformer_properties.jl"))
cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)]
reforming_area_properties = return_common_reformer_properties(pipe_length, n_cells, cell_lengths_along_pipe)

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
reaction_names = keys(reforming_area_properties.reactions.reforming_reactions)
species_names = keys(reforming_area_properties.molecular_weights)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        rho = u"kg/m^3",
        cell_combined_stationary_thermal_mass = u"J/K",
        UA_to_environment_per_m = u"W/(m*K)",
        wattage_received_per_m = u"W/m",
        molar_concentrations = u"mol/m^3",
        species_mass_flows = u"kg/s",
        net_rates = u"mol/s",
        mass = u"kg",
        mass_face = u"kg",
        mass_evaporated = u"kg",
        superficial_velocity = u"m/s",
        per_cell_fraction_of_total_heat_received = u"1", 
        #a correlation will be developed to account for the fact that the heating wire inside the reactor is unevenly spaced
        measured_heater_wattage = u"W", 
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
        UA_to_environment_per_m_p0 = 0.0u"W/(m*K)", 
        UA_to_environment_per_m_p2 = 0.0u"W/(m^2*K)", 
        empty_reactor_thermal_mass = 0.0u"J/(m*K)",

        heater_poly_p3 = 0.0u"W/m^3", 
        heater_poly_p2 = 0.0u"W/m^2", 
        heater_poly_p0 = 0.0u"W", 
        
        TCs_UA_to_center_of_reactor = 0.0u"W/K", 
        #TC2_UA_to_center_of_reactor = 0.0u"W/K", 
        #TC3_UA_to_center_of_reactor = 0.0u"W/K", 
        #TC4_UA_to_center_of_reactor = 0.0u"W/K", 
        #TC5_UA_to_center_of_reactor = 0.0u"W/K",
    )
)

function common_physics_functions!(du, u, cell_id, vol)
    #vaporization_model!(du, u, cell_id, vol)
    #ergun_momentum_friction!(du, u, cell_id, vol)
end


inlet_mass_fractions = ComponentVector(
    methanol = 0.0001u"kg/kg",
    water = 1.0u"kg/kg",
    carbon_monoxide = 0.0001u"kg/kg",
    hydrogen = 0.001u"kg/kg",
    carbon_dioxide = 0.0001u"kg/kg",
    air = 0.0001u"kg/kg"
)

total_inlet_mass_fractions = sum(inlet_mass_fractions)
inlet_mass_fractions ./= total_inlet_mass_fractions

empty_mass_fractions = ComponentVector(
    methanol = 1e-20u"kg/kg",
    water = 1e-20u"kg/kg",
    carbon_monoxide = 1e-20u"kg/kg",
    hydrogen = 1e-6u"kg/kg",
    carbon_dioxide = 1e-20u"kg/kg",
    air = 1.0u"kg/kg"
)

total_empty_mass_fractions = sum(empty_mass_fractions)
empty_mass_fractions ./= total_empty_mass_fractions

reforming_area_properties.pipe_mass_flow = 10.0u"g/minute" |> u"kg/s"

Revise.includet(joinpath(@__DIR__, "recorded_data_processing", "values_of_note.jl"))
values_of_note = get_trial_1_values_of_note()

reforming_area_properties.measured_room_temp = values_of_note.room_temperature_at_start
reforming_area_properties.pipe_mass_flow = values_of_note.flow_rate_throughout_trial .* 998u"kg/m^3"

Revise.includet(joinpath(@__DIR__, "recorded_data_processing", "inlet_and_outlet_temperatures.jl"))
inlet_and_outlet_temperatures = get_inlet_and_outlet_temperature_correlations()
inlet_temp_interp = inlet_and_outlet_temperatures.inlet_temp_interp
outlet_temp_interp = inlet_and_outlet_temperatures.outlet_temp_interp

#reforming_area_properties.measured_room_temp = misc_properties.room_temp


#TODO: remember to turn off the heater and to set 
#the inlet mass flow to what was observed in the experiment

add_region!(
    config, "inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = inlet_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 1.0,
        gas_holdup = 0.0,
        TC_temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]

        du.heat[cell_id] += (u.UA_to_environment_per_m[cell_id] * u.per_cell_pipe_length[cell_id]) * (u.measured_room_temp[cell_id] - u.temp[cell_id])

        du.heat[cell_id] += u.per_cell_fraction_of_total_heat_received[cell_id] * u.measured_heater_wattage[cell_id]

        du.heat[cell_id] *= 0.0
        #du.heat[cell_id] += u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id]

        sum_and_cap_fluxes!(du, u, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end
)

add_region!(
    config, "evaporator";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 0.0,
        gas_holdup = 1.0,
        TC_temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function evaporator!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.heat[cell_id] += (u.UA_to_environment_per_m[cell_id] * u.per_cell_pipe_length[cell_id]) * (u.measured_room_temp[cell_id] - u.temp[cell_id])

        du.heat[cell_id] += u.per_cell_fraction_of_total_heat_received[cell_id] * u.measured_heater_wattage[cell_id]

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "reactor";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 0.0,
        gas_holdup = 1.0,
        TC_temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function reactor!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        du.heat[cell_id] += (u.UA_to_environment_per_m[cell_id] * u.per_cell_pipe_length[cell_id]) * (u.measured_room_temp[cell_id] - u.temp[cell_id])

        du.heat[cell_id] += u.per_cell_fraction_of_total_heat_received[cell_id] * u.measured_heater_wattage[cell_id]

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "outlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 0.0,
        gas_holdup = 1.0,
        TC_temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function outlet!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.heat[cell_id] += (u.UA_to_environment_per_m[cell_id] * u.per_cell_pipe_length[cell_id]) * (u.measured_room_temp[cell_id] - u.temp[cell_id])

        du.heat[cell_id] += u.per_cell_fraction_of_total_heat_received[cell_id] * u.measured_heater_wattage[cell_id]

        du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]

        for_fields!(u.mass_fractions, du.species_mass_flows) do species, u_mass_fractions, du_species_mass_flows
            du_species_mass_flows[species[cell_id]] -= u.pipe_mass_flow[cell_id] * u_mass_fractions[species[cell_id]]
        end
        #this is to prevent the concentration of all species from building up at the outlet

        du.heat[cell_id] -= u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id] 
        #this is to prevent the temperature from building up at the outlet

        #du.heat[cell_id] *= 0.0

        sum_and_cap_fluxes!(du, u, cell_id, vol)
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

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

#you can check units by setting check_units = true and du0_vec and u0_vec will be returned as unitful ComponentVectors
#VSCodeServer.@profview 
du0_vec, u0_vec, state_axes, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

#TODO: create an interpolation for heater wattage using data from an arduino measuring voltage and current input

Revise.includet(joinpath(@__DIR__, "..", "..", "common_reformer_properties", "fraction_of_total_heat_received_correl.jl"))

pipe_mass_flow = ustrip(upreferred(values_of_note.flow_rate_throughout_trial .* 998u"kg/m^3"))

inlet_temp_interp = inlet_and_outlet_temperatures.inlet_temp_interp

pump_shutoff_timestamp = ustrip(values_of_note.pump_shut_off_time)
#I wonder if it's fine to put ustrip(values_of_note.pump_shut_off_time in the function itself)

inlet_temp_interp(0.0)
inlet_temp_interp(100.0)

function pump_shut_off(du, u, cell_id, t)
    if t <= pump_shutoff_timestamp #pump on
        #FOR FUTURE REFERENCE JUST SO YOU KNOW WHAT'S HAPPENING:
        #=if eltype(u.temp[1]) <: ForwardDiff.Dual && eltype(t) <: ForwardDiff.Dual
            u.temp[1] = inlet_temp_interp(ForwardDiff.value(t))
        elseif eltype(t) <: ForwardDiff.Dual
            u.temp[1] = inlet_temp_interp(ForwardDiff.value(t))
        elseif eltype(u.temp[1]) <: ForwardDiff.Dual
            u.temp[1] = inlet_temp_interp(t)
        else
            u.temp[1] = inlet_temp_interp(t) 
            #TODO: figure out this nonsense
            #why is this required, this is never required anywhere else
        end=#
        
        u.temp[1] = inlet_temp_interp(ForwardDiff.value(t)) 
        #for anything that uses t for an interpolation, make sure to get the value of it to prevent Dual shenanigans
        
        u.pipe_mass_flow[cell_id] = pipe_mass_flow
    else #pump shut off
        #do nothing to the inlet temp
        u.pipe_mass_flow[cell_id] = 0.0
    end
end

function solve_system!(du, u, p, t, geo, system)
    #VERY IMPORTANT: since most software uses 0-based indexing, you need to adjust the cell id by +1
    #for example, if you mouse over cell_id 5161 in paraview, you need to use 5162 in the code because julia uses 1-based indexing 

    for cell_id in eachindex(geo.cell_volumes)
        update_properties!(du, u, cell_id, geo.cell_volumes[cell_id])
        
        u.measured_heater_wattage[cell_id] = measured_heater_wattage(t)

        pump_shut_off(du, u, cell_id, t)
    end

    for cell_id in 1:length(geo.cell_volumes)-1
        #if u.rho[cell_id] >= 950.0 #if a cell has filled up, allow it to transfer mass into the next cell (which should be empty)
        #wow, this if statement really fucks with the KrylovJL_GMRES solver
        #I guess it likely has something to do with the fact that ever time a cell unfreezes, the matrix changes dramatically and the stiffness becomes highly variable
        #I'm now realizing that this is why the solve times have been so damn long for such a simple system
        #there must be if statements somewhere else (probably in upwinding or somehwere in the physics) that are causing this slowdown
        #this makes my lower solve times in the past make sense dspite the system being more complex 
            du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]
            du.mass_face[cell_id + 1, 1] += u.pipe_mass_flow[cell_id]
        #end
    end

    #this is a mess, ugh, imagine if we had 20 or 50 TCs!
    # K/s = (W/K * K) / (J/K)
    TC1_closest_cell_id = Int(u.TC1_closest_cell_id[1])
    TC1_heat_flux = u.UA_to_environment_per_m[TC1_closest_cell_id] * (u.temp[TC1_closest_cell_id] - u.TC_temp[1])
    du.TC_temp[1] += TC1_heat_flux / (u.empty_reactor_thermal_mass[1] * u.per_cell_pipe_length[TC1_closest_cell_id])
    du.temp[TC1_closest_cell_id] -= TC1_heat_flux / (u.cp[TC1_closest_cell_id] * u.per_cell_pipe_length[TC1_closest_cell_id])

    TC2_closest_cell_id = Int(u.TC2_closest_cell_id[1])
    TC2_heat_flux = u.UA_to_environment_per_m[TC2_closest_cell_id] * (u.temp[TC2_closest_cell_id] - u.TC_temp[2])
    du.TC_temp[2] += TC2_heat_flux / (u.empty_reactor_thermal_mass[1] * u.per_cell_pipe_length[TC2_closest_cell_id])
    du.temp[TC2_closest_cell_id] -= TC2_heat_flux / (u.cp[TC2_closest_cell_id] * u.per_cell_pipe_length[TC2_closest_cell_id])

    TC3_closest_cell_id = Int(u.TC3_closest_cell_id[1])
    TC3_heat_flux = u.UA_to_environment_per_m[TC3_closest_cell_id] * (u.temp[TC3_closest_cell_id] - u.TC_temp[3])
    du.TC_temp[3] += TC3_heat_flux / (u.empty_reactor_thermal_mass[1] * u.per_cell_pipe_length[TC3_closest_cell_id])
    du.temp[TC3_closest_cell_id] -= TC3_heat_flux / (u.cp[TC3_closest_cell_id] * u.per_cell_pipe_length[TC3_closest_cell_id])

    TC4_closest_cell_id = Int(u.TC4_closest_cell_id[1])
    TC4_heat_flux = u.UA_to_environment_per_m[TC4_closest_cell_id] * (u.temp[TC4_closest_cell_id] - u.TC_temp[4])
    du.TC_temp[4] += TC4_heat_flux / (u.empty_reactor_thermal_mass[1] * u.per_cell_pipe_length[TC4_closest_cell_id])
    du.temp[TC4_closest_cell_id] -= TC4_heat_flux / (u.cp[TC4_closest_cell_id] * u.per_cell_pipe_length[TC4_closest_cell_id])

    TC5_closest_cell_id = Int(u.TC5_closest_cell_id[1])
    TC5_heat_flux = u.UA_to_environment_per_m[TC5_closest_cell_id] * (u.temp[TC5_closest_cell_id] - u.TC_temp[5])
    du.TC_temp[5] += TC5_heat_flux / (u.empty_reactor_thermal_mass[1] * u.per_cell_pipe_length[TC5_closest_cell_id])
    du.temp[TC5_closest_cell_id] -= TC5_heat_flux / (u.cp[TC5_closest_cell_id] * u.per_cell_pipe_length[TC5_closest_cell_id])

    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = ustrip.(Vector(ComponentVector(
    UA_to_environment_per_m_p0 = 0.5u"W/(m*K)", 
    UA_to_environment_per_m_p2 = 1.0u"W/(m^2*K)", 
    empty_reactor_thermal_mass = 1.0u"J/(m*K)",

    #is a heater polynomial ever going to be accurate considering how 
    #uneven the spacing between each turn of heating wire is?
    #heater_poly_p3 = 0.0u"W/m^3", 
    #heater_poly_p2 = 0.0u"W/m^2", 

    TCs_UA_to_center_of_reactor = 1.0u"W/K", 
)))

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

Revise.includet(joinpath(@__DIR__, "thermocouple_data_processing", "thermocouple_data.jl"))
path_to_thermocouple_data = joinpath(@__DIR__, "thermocouple_data_processing", "reactor_data_2026_05_02__15_56_13.csv")
thermocouple_data = get_thermocouple_data(path_to_thermocouple_data)

t0 = 0.0
tMax = ustrip(thermocouple_data.timestamps[end])
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization()), callback = approximate_time_to_finish_cb)
#woah, wtf, why is FBDF so much faster all of a sudden?
#even KLUFactorization is way faster than KrylovJL_GMRES
#KrylovJL_BICGSTAB is another option, so is using algebraicmultigrid

@time sol = solve(implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)
#1.366 s (700880 allocations: 78.07 MiB)

#@time sol = solve(implicit_prob, Rodas5P(linsolve = KLUFactorization()), callback = approximate_time_to_finish_cb)

#prob = ODEProblem(ode_func, u0_vec, (0.0, 1e-5), p_guess)
#@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)]

sim_file = @__FILE__

#sol_to_vtk(sol, u_named, grid, sim_file)

timestamps = ustrip.(thermocouple_data.timestamps)

function loss(θ)
    #prob = ODEProblem(ode_func, u0_vec, (0.0, ustrip(thermocouple_data.timestamps[end])), θ)

    #loss_prob = remake(implicit_prob, p = vcat(θ, [100.0, 100.0, 100.0, 100.0, 100.0]))

    loss_prob = remake(implicit_prob, p = θ)

    #FBDF or Rodas5P works well here
    
    sol = solve(
        loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        #InterpolatingAdjoint(autodiff = AutoMooncake()),
        #callback = approximate_time_to_finish_cb
    )

    u_named = [ComponentVector(sol.u[i], state_axes) for i in eachindex(sol.u)]

    mean_squared_error = 0.0

    for i in eachindex(sol.t)
        mean_squared_error += abs2(ustrip(u_named[i].TC_temp[1]) - ustrip(thermocouple_data.TC1_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temp[2]) - ustrip(thermocouple_data.TC2_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temp[3]) - ustrip(thermocouple_data.TC3_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temp[4]) - ustrip(thermocouple_data.TC4_temps_interp(sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temp[5]) - ustrip(thermocouple_data.TC5_temps_interp(sol.t[i])))
        #mean_squared_error += abs2(ustrip(u_named[i].temp[end]) - ustrip(outlet_temp_interp(sol.t[i])))
        #I don't think using the measured outlet temperature is actually useful
    end

    return mean_squared_error
end

p_guess_init = ComponentVector(
    UA_to_environment_per_m_p0 = 0.5u"W/(m*K)", 
    UA_to_environment_per_m_p2 = 0.01u"W/(m^2*K)", 
    empty_reactor_thermal_mass = 1000.0u"J/(m*K)", 

    #heater_poly_p2 = 0.001u"W/m^3", 
    #heater_poly_p1 = 0.0u"W/m^2", 
    #heater_poly_p0 = 1312.0u"W/m", 
    
    TCs_UA_to_center_of_reactor = 0.0005u"W/K", #I think the units here should 1/s
    #TC2_UA_to_center_of_reactor = 0.0005u"W/K",
    #TC3_UA_to_center_of_reactor = 0.0005u"W/K",
    #TC4_UA_to_center_of_reactor = 0.0005u"W/K",
    #TC5_UA_to_center_of_reactor = 0.0005u"W/K"
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
    UA_to_environment_per_m_p0 = 0.001u"W/(m*K)", 
    UA_to_environment_per_m_p2 = 0.0001u"W/(m^2*K)", 
    empty_reactor_thermal_mass = 0.001u"J/(m*K)", 

    #heater_poly_p2 = 0.001u"W/m^3", 
    #heater_poly_p1 = 0.0u"W/m^2", 
    #heater_poly_p0 = 1312.0u"W/m", 
    
    TCs_UA_to_center_of_reactor = 0.0001u"W/K", 
    #TC2_UA_to_center_of_reactor = 0.0001u"W/K", 
    #TC3_UA_to_center_of_reactor = 0.0001u"W/K", 
    #TC4_UA_to_center_of_reactor = 0.0001u"W/K", 
    #TC5_UA_to_center_of_reactor = 0.0001u"W/K",
))))

p_upper_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    UA_to_environment_per_m_p0 = 100.0u"W/(m*K)", 
    UA_to_environment_per_m_p2 = 0.1u"W/(m^2*K)", 
    empty_reactor_thermal_mass = 10000.0u"J/(m*K)", 

    #heater_poly_p2 = 0.001u"W/m^3", 
    #heater_poly_p1 = 0.0u"W/m^2", 
    #heater_poly_p0 = 1312.0u"W/m", 
    
    TCs_UA_to_center_of_reactor = 10.0u"W/K", #maybe more like 0.1
    #TC2_UA_to_center_of_reactor = 10.0u"W/K", 
    #TC3_UA_to_center_of_reactor = 10.0u"W/K", 
    #TC4_UA_to_center_of_reactor = 10.0u"W/K", 
    #TC5_UA_to_center_of_reactor = 10.0u"W/K",
))))

optprob = Optimization.OptimizationProblem(optf, p_guess, lb=p_lower_bounds, ub=p_upper_bounds)

function randomize(lower, upper)
    return lower + (upper - lower) * rand()
end

p_ensemble = [[randomize(p_lower_bounds[i], p_upper_bounds[i]) for i in eachindex(p_lower_bounds)] for _ in 1:Sys.CPU_THREADS]

function prob_func(prob, i, repeat)
    return remake(prob, p = p_ensemble[i])  
end

ensembleprob = EnsembleProblem(optprob; prob_func)

LOSS = Float64[]
PARS = []

cb = function (state, l)
    display(l)
    display(state.u)
    append!(LOSS, l)
    append!(PARS, [state.u])
    false
end

@time res = Optimization.solve(
    optprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    f_abstol=1e-8,
    g_abstol=1e-8,
)

using OptimizationBBO #for BlackBoxOptim
@time res = solve(
    ensembleprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    EnsembleThreads(),
    trajectories = Sys.CPU_THREADS,
    callback = cb,
)

res.u0_vec

loss(ustrip.(upreferred.(Vector(p_vec_fitted))))

p_fitted = ComponentVector(res.u, p_axes)

#0.41741
#0.045
#7050
#8.71

