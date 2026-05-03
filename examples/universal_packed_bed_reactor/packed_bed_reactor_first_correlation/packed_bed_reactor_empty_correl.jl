using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using NonlinearSolve

using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ForwardDiff
using DataFrames
using CSV
using DataInterpolations

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
    TC_temps = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

van_t_hoff_A = ComponentVector(CH3O = 1.7e-6u"s^-1", HCOO = 4.74e-13u"s^-1", OH = 3.32e-14u"s^-1")
van_t_hoff_dH = ComponentVector(CH3O = -46800.0u"J/mol", HCOO = -115000.0u"J/mol", OH = -110000.0u"J/mol")

MSR_rxn = ComponentVector(
    heat_of_reaction = 49500.0u"J/mol", 
    ref_delta_G = -3800.0u"J/mol", 
    ref_temp = 298.15u"K", 
    kf_A = 1.59e10u"s^-1", #sources online point to values around 1.25e7 mol / (kg * s * bar)
    kf_Ea = 103000.0u"J/mol",
    reactant_stoich_coeffs = (methanol = 1, water = 1),
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3), 
    stoich_coeffs = (methanol = -1, water = -1, carbon_dioxide = 1, hydrogen = 3), # CH3OH + H2O ⇋ CO2 + 3(H2)
    van_t_hoff_A = van_t_hoff_A, 
    van_t_hoff_dH = van_t_hoff_dH
)

MD_rxn = ComponentVector(
    heat_of_reaction = 90200.0u"J/mol", 
    ref_delta_G = 24800.0u"J/mol", 
    ref_temp = 298.15u"K", 
    kf_A = 1.46e13u"s^-1", #sources online point to values around 1.15e11 mol / (kg * s * bar)
    kf_Ea = 170000.0u"J/mol",
    reactant_stoich_coeffs = (methanol = 1,), 
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), 
    stoich_coeffs = (methanol = -1, carbon_monoxide = 1, hydrogen = 2), # CH3OH ⇋ CO + 2(H2)
    van_t_hoff_A = van_t_hoff_A,
    van_t_hoff_dH = van_t_hoff_dH
)

WGS_rxn = ComponentVector(
    heat_of_reaction = -41100.0u"J/mol", 
    ref_delta_G = -28600.0u"J/mol", 
    ref_temp = 298.15u"K", 
    kf_A = 4.63e10u"s^-1", #sources online point to values around 3.65e7 mol / (kg * s * bar)
    kf_Ea = 87500.0u"J/mol",
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), 
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), 
    stoich_coeffs = (carbon_monoxide = -1, water = -1, carbon_dioxide = 1, hydrogen = 1), # CO + H2O ⇋ CO2 + H2
    van_t_hoff_A = van_t_hoff_A, 
    van_t_hoff_dH = van_t_hoff_dH
)

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

pipe_area = pi * (pipe_inside_diameter / 2)^2

cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)] 
#this only works because the base of the pipe is at z = 0.0

TC1_position_along_reactor = 1.5u"inch"
TC1_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC1_position_along_reactor))

TC2_position_along_reactor = 3.0u"inch"
TC2_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC2_position_along_reactor))

TC3_position_along_reactor = 5.0u"inch"
TC3_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC3_position_along_reactor))

TC4_position_along_reactor = 7.5u"inch"
TC4_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC4_position_along_reactor))

TC5_position_along_reactor = 10.0u"inch"
TC5_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC5_position_along_reactor))

reforming_area_properties = ComponentVector(
    k = 0.026u"W/(m*K)", 
    cp = 1006u"J/(kg*K)",
    dynamic_viscosity = 1.81e-5u"Pa*s",
    rho = 1.225u"kg/m^3",
    R_gas = 8.314u"J/(mol*K)",

    pipe_mass_flow = 0.0u"g/minute",

    pipe_inside_diameter = pipe_inside_diameter,
    pipe_area = pipe_area,
    per_cell_pipe_length = pipe_length / n_cells,
    cell_lengths_along_pipe = cell_lengths_along_pipe,

    bed_void_fraction = 0.4,
    packing_surface_area = 100.0u"m^2/m^3",
    particle_diameter = 5.0u"mm",

    overall_heat_transfer_coefficient = 1000.0u"W/(m^2*K)",
    measured_room_temp = 18.0u"°C", #change for each experiment

    TC1_bias = 0.0u"°C",
    TC2_bias = 0.0u"°C",
    TC3_bias = 0.0u"°C",
    TC4_bias = 0.0u"°C",
    TC5_bias = 0.0u"°C",

    TC1_closest_cell_id = TC1_closest_cell_id,
    TC2_closest_cell_id = TC2_closest_cell_id,
    TC3_closest_cell_id = TC3_closest_cell_id,
    TC4_closest_cell_id = TC4_closest_cell_id,
    TC5_closest_cell_id = TC5_closest_cell_id,

    external_temp = 300.0u"°C",
    saturation_temp = 72.4u"°C",
    liquid_rho = 791.0u"kg/m^3",
    gas_rho = 1.225u"kg/m^3",
    mass_transfer_coeff_vap = 0.001u"kg/(m^2*s*K)",
    heat_of_vaporization = 1.5u"kJ/g",
    liquid_feed_mass_fractions = empty_mass_fractions,

    diffusion_coefficients = (
        methanol = 1e-5u"m^2/s",
        water = 1e-5u"m^2/s",
        carbon_monoxide = 1e-5u"m^2/s",
        hydrogen = 1e-5u"m^2/s",
        carbon_dioxide = 1e-5u"m^2/s",
        air = 1e-5u"m^2/s"
    ), 
    molecular_weights = (
        methanol = 32.04u"g/mol",
        water = 18.02u"g/mol",
        carbon_monoxide = 28.01u"g/mol",
        hydrogen = 2.02u"g/mol",
        carbon_dioxide = 44.01u"g/mol",
        air = 28.97u"g/mol"
    ), 
    reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
    reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0u"kg/m^3", MD_rxn = 1250.0u"kg/m^3", WGS_rxn = 1250.0u"kg/m^3"),), 
)

permeability = (reforming_area_properties.particle_diameter^2 * reforming_area_properties.bed_void_fraction^3) / (150.0 * (1.0 - reforming_area_properties.bed_void_fraction)^2)

reforming_area_properties = merge_properties(reforming_area_properties, ComponentVector(permeability = permeability))

struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "..", "physics", "multiphase.jl")) #the .. makes it go up one directory
Revise.includet(joinpath(@__DIR__, "..", "physics", "energy.jl"))
Revise.includet(joinpath(@__DIR__, "..", "physics", "momentum.jl"))

function update_properties!(du, u, cell_id, vol)
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

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
reaction_names = keys(reforming_area_properties.reactions.reforming_reactions)
species_names = keys(reforming_area_properties.molecular_weights)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        rho = u"kg/m^3",
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
        #this will probably be measured by a kill-a-watt meter 
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
        UA_to_environment = 1.0u"W/K", 
        empty_reactor_thermal_mass = 1.0u"J/K",

        heater_poly_p2 = 0.0u"W/m^3", 
        heater_poly_p1 = 0.0u"W/m^2", 
        heater_poly_p0 = 0.0u"W/m", 
        
        TC1_UA_to_center_of_reactor = 0.0u"W/K", 
        TC2_UA_to_center_of_reactor = 0.0u"W/K", 
        TC3_UA_to_center_of_reactor = 0.0u"W/K", 
        TC4_UA_to_center_of_reactor = 0.0u"W/K", 
        TC5_UA_to_center_of_reactor = 0.0u"W/K",
    )
)

function common_physics_functions!(du, u, cell_id, vol)
    #vaporization_model!(du, u, cell_id, vol)
    #ergun_momentum_friction!(du, u, cell_id, vol)
end

add_region!(
    config, "inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 0.0,
        gas_holdup = 1.0,
        TC_temps = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.heat[cell_id] = u.UA_to_environment[1] * (u.measured_room_temp[cell_id] - u.temp[cell_id])

        du.heat[cell_id] += u.per_cell_fraction_of_total_heat_received[cell_id] * u.measured_heater_wattage[cell_id]

        sum_and_cap_fluxes!(du, u, cell_id, vol)
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
        TC_temps = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function evaporator!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        #UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        #surface_area = pi * u.pipe_inside_diameter[cell_id] * u.per_cell_pipe_length[cell_id]
        #du.heat[cell_id] += u.overall_heat_transfer_coefficient[cell_id] * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area

        du.heat[cell_id] = u.UA_to_environment[1] * (u.measured_room_temp[cell_id] - u.temp[cell_id])

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
        TC_temps = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function reactor!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        du.heat[cell_id] = u.UA_to_environment[1] * (u.measured_room_temp[cell_id] - u.temp[cell_id])

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
        TC_temps = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    region_function =
    function outlet!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.heat[cell_id] = u.UA_to_environment[1] * (u.measured_room_temp[cell_id] - u.temp[cell_id])

        du.heat[cell_id] += u.per_cell_fraction_of_total_heat_received[cell_id] * u.measured_heater_wattage[cell_id]

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    new_pressure_driven_mass_flux!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
    
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
du0_vec, u0_vec, state_axes, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

Revise.includet(joinpath(@__DIR__, "fraction_of_total_heat_received_correl.jl"))

#experimental_data = CSV.read(joinpath(@__DIR__, "experimental_data.csv"), DataFrame)
#experimental_data_time_to_idx = Dict([experimental_data[1, i] => i for i in eachindex(experimental_data[1, :])])

n_samples = 100

time_cols = [-0.001u"s" + i * 0.001u"s" for i in 1:n_samples]
temp_cols = hcat([fill(273.15u"K", n_samples), fill(300.15u"K", n_samples), fill(350.15u"K", n_samples), fill(300.15u"K", n_samples), fill(273.15u"K", n_samples)]...)
experimental_data = hcat(time_cols, temp_cols)

#TODO: create an interpolation for heater wattage using data from an arduino measuring voltage and current input

function solve_system!(du, u, p, t, geo, system)
    #VERY IMPORTANT: since most software uses 0-based indexing, you need to adjust the cell id by +1
    #for example, if you mouse over cell_id 5161 in paraview, you need to use 5162 in the code because julia uses 1-based indexing 

    for cell_id in eachindex(geo.cell_volumes)
        update_properties!(du, u, cell_id, geo.cell_volumes[cell_id])
        
        u.measured_heater_wattage[cell_id] = measured_heater_wattage(t)
    end

    #this is a mess, ugh, imagine if we had 20 or 50 TCs!
    TC1_closest_cell_id = Int(u.TC1_closest_cell_id[1])
    du.TC_temps[1] += u.TC1_UA_to_center_of_reactor[1] * (u.TC_temps[1] - u.temp[TC1_closest_cell_id])
    #hmm, I think we're going to have to have TC1_UA_to_center_of_reactor be a lumped parameter with thermal mass 

    TC2_closest_cell_id = Int(u.TC2_closest_cell_id[1])
    du.TC_temps[2] += u.TC2_UA_to_center_of_reactor[1] * (u.TC_temps[2] - u.temp[TC2_closest_cell_id])

    TC3_closest_cell_id = Int(u.TC3_closest_cell_id[1])
    du.TC_temps[3] += u.TC3_UA_to_center_of_reactor[1] * (u.TC_temps[3] - u.temp[TC3_closest_cell_id])

    TC4_closest_cell_id = Int(u.TC4_closest_cell_id[1])
    du.TC_temps[4] += u.TC4_UA_to_center_of_reactor[1] * (u.TC_temps[4] - u.temp[TC4_closest_cell_id])

    TC5_closest_cell_id = Int(u.TC5_closest_cell_id[1])
    du.TC_temps[5] += u.TC5_UA_to_center_of_reactor[1] * (u.TC_temps[5] - u.temp[TC5_closest_cell_id])

    calculate_fraction_of_total_heat_received!(du, u, geo)

    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = ustrip.(Vector(ComponentVector(
    UA_to_environment = 1.0u"W/K", 
    empty_reactor_thermal_mass = 1.0u"J/K",

    heater_poly_p2 = 0.001u"W/m^3", 
    heater_poly_p1 = 0.0u"W/m^2", 
    heater_poly_p0 = 1312.0u"W/m", 
    
    TC1_UA_to_center_of_reactor = 1.0u"W/K", 
    TC2_UA_to_center_of_reactor = 1.0u"W/K", 
    TC3_UA_to_center_of_reactor = 1.0u"W/K", 
    TC4_UA_to_center_of_reactor = 1.0u"W/K", 
    TC5_UA_to_center_of_reactor = 1.0u"W/K",
)))

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = ustrip(experimental_data[end, 1])
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), concrete_jac = true), callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, AutoTsit5(FBDF()), callback = approximate_time_to_finish_cb)
#KrylovJL_BICGSTAB is another option, so is using algebraicmultigrid

#prob = ODEProblem(ode_func, u0_vec, (0.0, 1e-5), p_guess)
#@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)]

sim_file = @__FILE__

#sol_to_vtk(sol, u_named, grid, sim_file)

function loss(θ)
    prob = ODEProblem(ode_func, u0_vec, (0.0, ustrip(experimental_data[end, 1])), θ)
    
    sol = solve(
        prob, 
        FBDF(linsolve = KrylovJL_GMRES()), 
        sensealg = InterpolatingAdjoint(autodiff = AutoForwardDiff()),
        callback = approximate_time_to_finish_cb
    )

    u_named = [ComponentVector(sol.u[i], state_axes) for i in eachindex(sol.u)]

    #=mean_squared_error = 0.0

    for i in eachindex(experimental_data[:, 1])
        # Assuming sol.u matches experimental_data indices, but it's better to use interpolation
        # For now, let's just fix the indexing assuming sol.u matches
        mean_squared_error += abs2(ustrip(u_named[i].TC_temps[1]) - ustrip(experimental_data[i, 2]))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temps[2]) - ustrip(experimental_data[i, 3]))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temps[3]) - ustrip(experimental_data[i, 4]))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temps[4]) - ustrip(experimental_data[i, 5]))
        mean_squared_error += abs2(ustrip(u_named[i].TC_temps[5]) - ustrip(experimental_data[i, 6]))
    end

    return mean_squared_error / length(experimental_data[:, 1])=#

    #return [u_named[i].TC_temps[1] for i in eachindex(u_named)]
    return u_named[end].TC_temps[1]
    #divide by number of data points to not favour experiments that ran for longer
end

p_guess = ComponentVector(
    UA_to_environment = 1.0u"W/K", 
    empty_reactor_thermal_mass = 1.0u"J/K",

    heater_poly_p2 = 0.001u"W/m^3", 
    heater_poly_p1 = 0.0u"W/m^2", 
    heater_poly_p0 = 1312.0u"W/m", 
    
    TC1_UA_to_center_of_reactor = 1.0u"W/K", 
    TC2_UA_to_center_of_reactor = 1.0u"W/K", 
    TC3_UA_to_center_of_reactor = 1.0u"W/K", 
    TC4_UA_to_center_of_reactor = 1.0u"W/K", 
    TC5_UA_to_center_of_reactor = 100.0u"W/K",
)

grad = ForwardDiff.gradient(loss, ustrip.(Vector(p_guess)))

jac = ForwardDiff.jacobian(loss, ustrip.(Vector(p_guess)))

using Zygote
grad_zygote = Zygote.gradient(loss, ustrip.(Vector(p_guess)))

jac_zygote = Zygote.jacobian(loss, ustrip.(Vector(p_guess)))

jac[1, :]
jac[round(Int, (length(jac[:, 1])/2)), :]
jac[end, :]

plot(1:length(jac[:, 1]), jac[:, 1], label = "UA_to_environ")
plot!(1:length(jac[:, 2]), jac[:, 2], label = "empty_reactor_thermal_mass")
plot!(1:length(jac[:, 3]), jac[:, 3], label = "heater_poly_p2")
plot!(1:length(jac[:, 4]), jac[:, 4], label = "heater_poly_p1")
plot!(1:length(jac[:, 5]), jac[:, 5], label = "heater_poly_p0")
plot!(1:length(jac[:, 6]), jac[:, 6], label = "TC1_UA_to_center_of_reactor")
plot!(1:length(jac[:, 7]), jac[:, 7], label = "TC2_UA_to_center_of_reactor")
plot!(1:length(jac[:, 8]), jac[:, 8], label = "TC3_UA_to_center_of_reactor")
plot!(1:length(jac[:, 9]), jac[:, 9], label = "TC4_UA_to_center_of_reactor")
plot!(1:length(jac[:, 10]), jac[:, 10], label = "TC5_UA_to_center_of_reactor")

plot(jac[:, 1])





sense_prob = ODEForwardSensitivityProblem(ode_func, u0_vec, tspan, ustrip.(Vector(p_guess)))

@time sol = solve(
    sense_prob, 
    Tsit5(), 
    sensealg = InterpolatingAdjoint(),
    callback = approximate_time_to_finish_cb,
    tspan = (0.0, 1e-5)
)

sensitivities = extract_local_sensitivities(sol)

sensitivities[1]