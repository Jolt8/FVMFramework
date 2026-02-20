using Revise
using Logging

ENV["JULIA_PKG_PRESERVE_TIERED_INSTALLED"] == true
ENV["JULIA_PKG_PRESERVE_TIERED_INSTALLED"] = true

using FVMFramework
#I think this take a long time because Enzyme (despite not being installed) throws a warning from DiffEqBaseEnzyme.jl

using Ferrite
using FerriteGmsh
using OrdinaryDiffEq
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays
using PreallocationTools
using ForwardDiff
using Polyester
using StatsBase

using Unitful

grid = togrid("C://Users//wille//Desktop//FreeCad Projects//Methanol Reformer//output.msh")

grid.cellsets

n_cells = length(grid.cells)

struct OptimizedVariable
    guess::Float64
end

function optimize!(var)
    return OptimizedVariable(var)
end

#optimize!() will be inserted into any variable in the setup functions below
#for example, if we wanted to fit kf_A and kf_Ea to experimental data we would just do this:
#then, our tracer function would check if typeof(var) == OptimizedVariable, and if it is, we would make a pointer to a p value 
#=

MSR_rxn = ComponentVector(
    delta_H = 49500.0, # [J/mol]
    delta_G_ref = -3800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = optimize!(1.25e7), # [s^-1] #if we wanted to fit kf_A to experimental data, we could create a p-pointer here like kf_a = optimized!(p_proto)
    kf_Ea = optimize!(103000.0), # [J/mol]
    reactant_ids = ComponentVector(methanol=1, water=2), # reactant_ids: Methanol, Water
    reactant_stoich_coeffs = ComponentVector(methanol=1, water=1), # reactant_stoich_coeffs
    product_ids = ComponentVector(carbon_dioxide=5, hydrogen=4),     # product_ids: CO2, Hydrogen
    product_stoich_coeffs = ComponentVector(carbon_dioxide=1, hydrogen=3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = ComponentVector(methanol=-1, water=-1, carbon_monoxide=0, hydrogen=3, carbon_dioxide=1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec = van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec = van_t_hoff_dH_vec = van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)
=#


#oh god, I just realized how batshit insane it's going to be combine nested values in state variables
#like if we wanted to use species = ComponentArray(methanol = ComponentVector(mass_fraction = 0.3, mole_fraction = 0.1, molar_mass = 33.03)),
#we would have to inject caches like mole_fraction  or fixed_values like molar_mass into u_merged in the function itself 
u_proto = ComponentArray(
    vel_x=zeros(n_cells), vel_y=zeros(n_cells), vel_z=zeros(n_cells),
    pressure=zeros(n_cells),
    mass_fractions=ComponentVector(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp=zeros(n_cells)
)

config = create_fvm_config(grid, u_proto)

#in the future, we may want to add a method to make all cells that are not under a cell set 
#become part of a default set with no internal physics or variables 

#controller logic
#We're just going to make this monitor one u field for now 
temp_controller = ComponentVector(
    proportional_gain=1.0,
    integral_time=1.0,
    derivative_time=1.0,
    desired_value=ustrip(270.0u"°C" |> u"K"),
    #desired_value_comp_vector=ComponentVector(temp=ustrip(270.0u"°C" |> u"K")), #let's try to implement this back in later, with the tracer it makes it difficult
    initial_volumetric_input=1e7,
    min_volumetric_input=0.0,
    max_volumetric_input=1e8
)

#this will probably have to be converted into an integrator internally to get previous temperature and to prevent excessive cache usage

#TODO: Another thing I'd like to implement in my eventual optimizaiton pipeline is a method to extract a very basic
#correlation between the average_temperature measured in the reforming_area cellset to another temp_sensor cellset 
#that could be fed into an arduino to extrapolate sensor data to the reactor's actual internal reforming temp
add_controller!(config;
    controller_properties=temp_controller,
    monitored_cellset="reforming_area",
    affected_cellset="heating_areas",
    controller_function=
    function pid_temp_controller(du, u, controller_id, monitored_cells, affected_cells, cell_volumes)
        measured_vec = u.temp
        measured_du_vec = du.temp

        measured_avg = 0.0
        measured_du_avg = 0.0

        @batch for monitored_cell_id in monitored_cells
            measured_avg += measured_vec[monitored_cell_id]
            measured_du_avg += measured_du_vec[monitored_cell_id]
        end

        measured_avg /= length(monitored_cells)
        measured_du_avg /= length(monitored_cells)

        error = measured_avg - u.controllers.desired_value[controller_id]

        du.integral_error[controller_id] = error

        corrected_volumetric_addition = (
            u.controllers.initial_volumetric_input[controller_id] +
            (u.controllers.proportional_gain[controller_id] * error) +
            (u.controllers.integral_time[controller_id] * u.integral_error[controller_id]) +
            (u.controllers.derivative_time[controller_id] * measured_du_avg)
        )

        corrected_volumetric_addition = clamp(corrected_volumetric_addition, u.controllers.min_volumetric_input[controller_id], u.controllers.max_volumetric_input[controller_id])

        @batch for affected_cell_id in affected_cells
            du.temp[affected_cell_id] += corrected_volumetric_addition * cell_volumes[affected_cell_id]
        end
    end
)

#for CH3OH, HCOO, OH
van_t_hoff_A_vec = ComponentVector(CH3O=1.7e-6, HCOO=4.74e-13, OH=3.32e-14)
van_t_hoff_dH_vec = ComponentVector(CH3O=-46800.0, HCOO=-115000.0, OH=-110000.0)

MSR_rxn = ComponentVector(
    delta_H=49500.0, # [J/mol]
    delta_G_ref=-3800.0, # [J/mol]
    ref_temp=298.15, # [K]
    kf_A=1.25e7, # [s^-1] #if we wanted to fit kf_A to experimental data, we could create a p-pointer here like kf_a = optimized!(p_proto)
    kf_Ea=103000.0, # [J/mol]
    reactant_ids=ComponentVector(methanol=1, water=2), # reactant_ids: Methanol, Water
    reactant_stoich_coeffs=ComponentVector(methanol=1, water=1), # reactant_stoich_coeffs
    product_ids=ComponentVector(carbon_dioxide=5, hydrogen=4),     # product_ids: CO2, Hydrogen
    product_stoich_coeffs=ComponentVector(carbon_dioxide=1, hydrogen=3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs=ComponentVector(methanol=-1, water=-1, carbon_monoxide=0, hydrogen=3, carbon_dioxide=1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec=van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec=van_t_hoff_dH_vec = van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

MD_rxn = ComponentVector(
    delta_H=90200.0, # [J/mol]
    delta_G_ref=24800.0, # [J/mol]
    ref_temp=298.15, # [K]
    kf_A=1.15e11, # [s^-1]
    kf_Ea=170000.0, # [J/mol]
    reactant_ids=ComponentVector(methanol=1), # reactant_ids: Methanol
    reactant_stoich_coeffs=ComponentVector(methanol=1), # reactant_stoich_coeffs
    product_ids=ComponentVector(carbon_monoxide=3, hydrogen=4), # product_ids: CO, Hydrogen
    product_stoich_coeffs=ComponentVector(carbon_monoxide=1, hydrogen=2), # product_stoich_coeffs: 1 CO + 2 H2
    stoich_coeffs=ComponentVector(methanol=-1, water=0, carbon_monoxide=1, hydrogen=2, carbon_dioxide=0), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec=van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec=van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

WGS_rxn = ComponentVector(
    delta_H=-41100.0, # [J/mol]
    delta_G_ref=-28600.0, # [J/mol]
    ref_temp=298.15, # [K]
    kf_A=3.65e7, # [s^-1]
    kf_Ea=87500.0, # [J/mol]
    reactant_ids=ComponentVector(carbon_monoxide=3, water=2), # reactant_ids: CO, Water
    reactant_stoich_coeffs=ComponentVector(carbon_monoxide=1, water=1), # reactant_stoich_coeffs
    product_ids=ComponentVector(carbon_dioxide=5, hydrogen=4), # product_ids: CO2, Hydrogen
    product_stoich_coeffs=ComponentVector(carbon_dioxide=1, hydrogen=1), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs=ComponentVector(methanol=0, water=-1, carbon_monoxide=-1, hydrogen=1, carbon_dioxide=1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec=van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec=van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

species_molecular_weights = ComponentVector(methanol=0.03204, water=0.01802, carbon_monoxide=0.02801, hydrogen=0.00202, carbon_dioxide=0.04401)

# methanol, water, carbon_monoxide, hydrogen, carbon_dioxide
initial_mass_fractions = ComponentVector(methanol=1.0, water=1.3, carbon_monoxide=0.0001, hydrogen=0.02, carbon_dioxide=0.0001)

initial_mass_fractions = initial_mass_fractions ./ sum(initial_mass_fractions)

reforming_area_properties = ComponentVector(
    k=237.0, # k (W/(m*K))
    cp=3000.0, # cp (J/(kg*K))
    mu=1e-5, # mu (Pa*s)
    permeability=6.0e-12, # permeability (m^2)
    diffusion_coefficients=ComponentVector(methanol=1e-5, water=1e-5, carbon_monoxide=1e-5, hydrogen=1e-5, carbon_dioxide=1e-5), #diffusion coefficients (m^2/s) 
    #FIXME: this is going to be tricky, diffusion coefficients should use the same [:, cell_id] that mass fractions use
    species_molecular_weights=species_molecular_weights,
    reactions=ComponentVector(
        reforming_reactions=ComponentVector(
            MSR_rxn=MSR_rxn,
            MD_rxn=MD_rxn,
            WGS_rxn=WGS_rxn
        )
    ),
    reactions_kg_cat=ComponentVector(MSR_rxn=1250.0, MD_rxn=1250.0, WGS_rxn=1250.0), # cell_kg_cat_per_m3_for_each_reaction
)

#this is another special case, I don't think it would be wise to store the chemical reaction struct within the u_vec
#thus, we could probably use another component array to store the properties of the reaction
#ex. reaction.delta_H, reaction.ref_temp

add_region!(
    config, "reforming_area";
    type=Fluid(),
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=0.0, vel_z=0.0,
        pressure=100000.0,
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    properties=reforming_area_properties,
    region_function=
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        #mw_avg!(u, cell_id)
        #rho_ideal!(u, cell_id)

        #internal physics

        #we have to figure out if we're going to pass in just a singular cell_volume or all cell_volumes and use cell_id
        #power_law_react_cell!(du, u, cell_id, u.reactions.example_reaction, vol)
        #example of how to do a power law reaction

        PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    end
)

config.u_proto

inlet_total_volume = get_cell_set_total_volume(grid, "inlet", config.geo)
#this should be 1e-8 [m^3], which it is

desired_m_dot = ustrip(0.2u"g/s" |> u"kg/s")
corrected_m_dot_per_volume = desired_m_dot / inlet_total_volume

#we should probably create methods to automatically copy the parameters from other already defined regions

add_region!(
    config, "inlet";
    type=Fluid(),
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=-1.0, vel_z=0.0,
        pressure=110000.0,
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    properties=reforming_area_properties,
    region_function=
    function inlet_area!(du, u, cell_id, vol)
        #property retrieval
        #mw_avg!(u, cell_id)
        #rho_ideal!(u, cell_id)

        #internal physics
        #=
        power_law_react_cell!(du, u, cell_id, u.reactions.example_reaction, vol)
        PAM_reforming_react_cell!(du, u, cell_id, vol)
        =#

        #sources
        du.pressure[cell_id] += corrected_m_dot_per_volume * vol

        #boundary conditions

        du.mass_fractions[:, cell_id] .= 0.0

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "outlet";
    type=Fluid(),
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=1.0, vel_z=0.0,
        pressure=90000.0,
        mass_fractions=ComponentVector(methanol=0.0, water=0.0, carbon_monoxide=0.0, hydrogen=0.0, carbon_dioxide=0.0),
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    properties=reforming_area_properties,
    region_function=
    function outlet_area!(du, u, phys, cell_id, vol)
        #property retrieval
        #mw_avg!(u, cell_id, phys.species_molecular_weights)
        #rho_ideal!(u, cell_id)

        #internal physics
        #=
        react_cell!(
            du, u, cell_id,
            vol, 
            phys
        )
            =#

        #sources
        du.pressure[cell_id] *= 0.0

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, phys)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol, phys)
    end
)

add_region!(
    config, "wall";
    type=Solid(),
    initial_conditions=ComponentVector(
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    properties=ComponentVector(
        k = 237.0, # k (W/(m*K))
        rho=2700.0, # rho (kg/m^3)
        cp=921.0, # cp (J/(kg*K))
    ),
    region_function=
    function wall_area!(du, u, phys, cell_id, vol)
        #property retrieval
        #oh, this is an issue, how do we pass in rho for solids?
        #rho = phys.rho

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, phys)
    end
)

#we might want to add something akin to distribute_over_set_volume!(du += input_wattage)
total_heating_volume = get_cell_set_total_volume(grid, "heating_areas", config.geo)

input_wattage = 100.0 # W
corrected_volumetric_heating = input_wattage / total_heating_volume

add_region!(
    config, "heating_areas";
    type=Solid(),
    initial_conditions=ComponentVector(
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    properties=ComponentVector(
        k=237.0, # k (W/(m*K))
        rho=2700.0, # rho (kg/m^3)
        cp=921.0, # cp (J/(kg*K))
    ),
    region_function=
    function heating_area!(du, u, phys, cell_id, vol)
        #property retrieval
        #rho = u.rho[cell_id]

        #internal physics

        #sources
        #corrected_volumetric_heating = 1.2351610076363146e7
        du.temp[cell_id] += corrected_volumetric_heating * vol

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, phys)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    #mw_avg!(u, idx_a)
    #rho_ideal!(u, idx_a)

    #mw_avg!(u, idx_b)
    #rho_ideal!(u, idx_b)

    #mutating-ish, it mutates du.pressure for a
    continuity_and_momentum_darcy!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    all_species_advection!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    enthalpy_advection!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    diffusion_mass_fraction_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function solid_solid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function fluid_solid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    #mw_avg!(u, idx_a, phys_a.species_molecular_weights, mw_avg_cache)
    #rho_ideal!(u, idx_a, rho_cache, mw_avg_cache)
    #rho_b = phys_b.rho

    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function solid_fluid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    #rho_a = phys_a.rho
    #mw_avg!(u, idx_b, phys_b.species_molecular_weights, mw_avg_cache)
    #rho_ideal!(u, idx_b, rho_cache, mw_avg_cache)

    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

#this is the smallest I could make this function
#we use types here because we can create subtypes of types to prevent having to do if a == "reforming_area" && b == "inlet" && return fluid_fluid_flux!
function conneciton_map_function(type_a, type_b)
    type_a <: Fluid && type_b <: Fluid && return fluid_fluid_flux!
    type_a <: Solid && type_b <: Solid && return solid_solid_flux!
    type_a <: Fluid && type_b <: Solid && return fluid_solid_flux!
    type_a <: Solid && type_b <: Fluid && return solid_fluid_flux!
end

du0, u0, geo, system = finish_fvm_config(config, conneciton_map_function)

N::Int = ForwardDiff.pickchunksize(length(u0))

rho_cache = zeros(length(grid.cells))
mw_avg_cache = zeros(length(grid.cells))

n_species = length(u_proto.mass_fractions[:, 1])
change_in_molar_concentrations_cache = zeros(n_species)
molar_concentrations_cache = zeros(n_species) #just using mass fractions for cell 1, this may cause some issues later!

n_reactions = length(reforming_area_physics.chemical_reactions)
net_rates_cache = zeros(n_reactions)

caches = DiffCache(
    ComponentVector(
        rho_cache=rho_cache,
        mw_avg_cache=mw_avg_cache,
        change_in_molar_concentrations_cache=change_in_molar_concentrations_cache,
        molar_concentrations_cache=molar_concentrations_cache,
        net_rates_cache=net_rates_cache
    ), N
)

du_caches = copy(caches)

#I was suprised that this works, if it doesn't in the future, we can still use caches_axes = getaxes(caches_proto)[1]; caches = ComponentVector(caches, cache_axes)

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t,
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,
    system.connection_groups, system.controller_groups, system.region_groups,
    du_caches, caches,
    system.u_merged_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

t0 = 0.0
tMax = 10.0
tspan = (t0, tMax)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

test_prob = ODEProblem(f_closure_implicit, u0, (0.0, 0.00000001), p_guess)
VSCodeServer.@profview sol = solve(test_prob, Tsit5(), callback=approximate_time_to_finish_cb)

#somehow adding reactions makes this less stiff

#holy moly, this is soooo stiff, even with just 0.000005 s of sim time, it has to take 1700 steps! It also took 356 seconds 

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

u_named = rebuild_u_named_vel(sol.u, u_proto)

grid.cellsets["reforming_area"]

mass_fractions_beginning = u_named[1].mass_fractions[:, 6389]

mass_fractions_end = u_named[end].mass_fractions[:, 6389]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end