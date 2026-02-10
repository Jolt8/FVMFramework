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

using Unitful

grid = togrid("C://Users//wille//Desktop//FreeCad Projects//Methanol Reformer//output.msh")

grid.cellsets

#for CH3OH, HCOO, OH
van_t_hoff_A_vec = [1.7e-6, 4.74e-13, 3.32e-14]
van_t_hoff_dH_vec = [-46800.0, -115000.0, -110000.0]

MSR_rxn = MSRReaction(
    49500.0,    # Delta H (Endothermic) [J/mol]
    -3800.0,    # Delta Gibbs free at 298.15K [J/mol]
    298.15,     # ref temp [K]
    1.25e7,     # kf_A (Pre-exponential factor)
    103000.0,   # kf_Ea (Activation Energy) [J/mol]
    [1, 2],     # reactant_ids: Methanol, Water
    [1, 1],     # reactant_stoich_coeffs
    [5, 4],     # product_ids: CO2, Hydrogen
    [1, 3],     # product_stoich_coeffs: 1 CO2 + 3 H2
    [-1, -1, 0, 3, 1], # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

MD_rxn = MDReaction(
    90200.0,    # Delta H (Endothermic) [J/mol]
    24800.0,    # Delta Gibbs free at 298.15K [J/mol]
    298.15,     # ref temp [K]
    1.15e11,    # kf_A
    170000.0,   # kf_Ea [J/mol]
    [1],        # reactant_ids: Methanol
    [1],        # reactant_stoich_coeffs
    [3, 4],     # product_ids: CO, Hydrogen
    [1, 2],     # product_stoich_coeffs: 1 CO + 2 H2
    [-1, 0, 1, 2, 0], # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

WGS_rxn = WGSReaction(
    -41100.0,   # Delta H (Exothermic) [J/mol]
    -28600.0,   # Delta Gibbs free at 298.15K [J/mol]
    298.15,     # ref temp [K]
    3.65e7,     # kf_A (Note: often adjusted depending on specific catalyst)
    87500.0,    # kf_Ea [J/mol]
    [3, 2],     # reactant_ids: CO, Water
    [1, 1],     # reactant_stoich_coeffs
    [5, 4],     # product_ids: CO2, Hydrogen
    [1, 1],     # product_stoich_coeffs
    [0, -1, -1, 1, 1], # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

species_molecular_weights = [0.03204, 0.01802, 0.02801, 0.00202, 0.04401]

# methanol, water, carbon_monoxide, hydrogen, carbon_dioxide
initial_mass_fractions = [1.0, 1.3, 0.0001, 0.02, 0.0001]

initial_mass_fractions = initial_mass_fractions ./ sum(initial_mass_fractions)

n_cells = length(grid.cells)
u_proto = ComponentArray(
    vel_x=zeros(n_cells), vel_y=zeros(n_cells), vel_z=zeros(n_cells),
    pressure=zeros(n_cells),
    mass_fractions=zeros(length(initial_mass_fractions), n_cells),
    temp=zeros(n_cells)
)

config = create_fvm_config(grid, u_proto)

#each region_function will have internal_physics (non-connecting), then sources, then boundary conditions, then capacities

reforming_area_physics = MethanolReformerPhysics(
    237.0, # k (W/(m*K))
    4.184, # cp (J/(kg*K))
    1e-5, # mu (Pa*s)
    0.6e-11, # permeability (m^2)
    [1e-5, 1e-5, 1e-5, 1e-5, 1e-5], #diffusion coefficients (m^2/s)
    species_molecular_weights,
    [MSR_rxn, MD_rxn, WGS_rxn], # chemical_reactions
    [1250.0, 1250.0, 1250.0], # cell_kg_cat_per_m3_for_each_reaction
)

#in the future, we may want to add a method to make all cells that are not under a cell set 
#become part of a default set with no internal physics or variables 


#controller logic
#We're just going to make this monitor one u field for now 
struct PIDController <: AbstractController
    proportional_gain::Float64
    integral_time::Float64
    derivative_time::Float64
    desired_value_comp_vector::ComponentArray
    initial_volumetric_input::Float64
    min_volumetric_input::Float64
    max_volumetric_input::Float64
end

temp_controller = PIDController(
    1, 1, 1, #P I D parameters 
    ComponentVector(temp=ustrip(270.0u"°C" |> u"K")), #field and desired value
    1e7, #initial watts / m^3
    0.0, 1e8 # min/max watts / m^3
)
integral_error = 0.0

#this will probably have to be converted into an integrator internally to get previous temperature and to prevent excessive cache usage

#TODO: Another thing I'd like to implement in my eventual optimizaiton pipeline is a method to extract a very basic
#correlation between the average_temperature measured in the reforming_area cellset to another temp_sensor cellset 
#that could be fed into an arduino to extrapolate sensor data to the reactor's actual internal reforming temp
add_controller!(config,
    controller=temp_controller,
    monitored_cellset="reforming_area",
    affected_cellset="heating_area",
    controller_function=function pid_temp_controller(
        du, u, monitored_cells, affected_cells, controller_id,
        controller, 
        cell_volumes
    )
        field = propertynames(controller.desired_value_comp_vector)[1]
        desired = getproperty(controller.desired_value_comp_vector, field) #desired will be a single scalar value
        measured_vec = getproperty(u, field)
        measured_du_vec = getproperty(du, field)

        measured_avg = 0.0
        measured_du_avg = 0.0

        @batch for monitored_cell_id in monitored_cells
            measured_avg += measured_vec[monitored_cell_id]
            measured_du_avg += measured_du_vec[monitored_cell_id]
        end

        measured_avg /= length(monitored_cells)
        measured_du_avg /= length(monitored_cells)

        error = measured_avg - desired

        du.integral_error[controller_id] = error

        corrected_volumetric_addition = (
            controller.initial_volumetric_input +
            (controller.proportional_gain * error) +
            (controller.integral_time * u.integral_error[controller_id]) +
            (controller.derivative_time * measured_du_avg)
        )

        corrected_volumetric_addition = clamp(corrected_volumetric_addition, controller.min_volumetric_input, controller.max_volumetric_input)

        du_field_vec = getproperty(du, field)

        @batch for affected_cell_id in affected_cells
            du_field_vec[affected_cell_id] += corrected_volumetric_addition * cell_volumes[affected_cell_id]
        end
    end
)

add_region!(
    config, "reforming_area";
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=0.0, vel_z=0.0,
        pressure=100000.0,
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics=reforming_area_physics,
    region_function=function reforming_area!(
        du, u, phys, cell_id, 
        vol,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )
        #property updating/retrieval
        mw_avg!(u, cell_id, phys.species_molecular_weights, mw_avg_cache)
        rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)

        rho = rho_cache[cell_id]
        #rho = @view rho_cache[cell_id][1] #not sure if the above allocates

        #internal physics

        react_cell!(
            du, u, cell_id,
            rho,
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol,
            phys
        )

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol, rho, phys)
    end
)

geo = build_fvm_geo_into_struct(grid)

geo.cell_volumes

geo.cell_face_areas

getcellset(grid, "inlet")
inlet_total_volume = get_cell_set_total_volume(grid, "inlet", geo)
#this should be 1e-8, which it is 

desired_m_dot = ustrip(0.2u"g/s" |> u"kg/s")
corrected_m_dot_per_volume = desired_m_dot / inlet_total_volume

#we should probably create methods to automatically copy the parameters from other already defined regions

add_region!(
    config, "inlet",
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=-1.0, vel_z=0.0,
        pressure=110000.0,
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics=reforming_area_physics,
    region_function=function inlet_area!(
        du, u, phys, cell_id, 
        vol,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )
        #property retrieval
        mw_avg!(u, cell_id, phys.species_molecular_weights, mw_avg_cache)
        rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)

        rho = rho_cache[cell_id]
        #rho = @view rho_cache[cell_id][1] #not sure if the above allocates

        #internal physics
        #=
        react_cell!(
            du, u, cell_id,
            rho_cache[cell_id],
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, 
            phys
        )
            =#

        #sources
        du.pressure[cell_id] += corrected_m_dot_per_volume * vol

        #boundary conditions

        du.mass_fractions[:, cell_id] .= 0.0

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol, rho, phys)
    end
)

add_region!(
    config, "outlet",
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=1.0, vel_z=0.0,
        pressure=90000.0,
        mass_fractions=[0.0, 0.0, 0.0, 0.0, 0.0],
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics=reforming_area_physics,
    region_function=function outlet_area!(
        du, u, phys, cell_id, 
        vol,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )
        #property retrieval
        mw_avg!(u, cell_id, phys.species_molecular_weights, mw_avg_cache)
        rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)

        rho = rho_cache[cell_id]
        #rho = @view rho_cache[cell_id][1] #not sure if the above allocates

        #internal physics
        #=
        react_cell!(
            du, u, cell_id,
            rho_cache[cell_id],
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, 
            phys
        )
            =#

        #sources
        du.pressure[cell_id] *= 0.0

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol, rho, phys)
    end
)

add_region!(
    config, "wall";
    initial_conditions=ComponentVector(
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics=WallPhysics(
        237.0, # k (W/(m*K))
        2700.0, # rho (kg/m^3)
        921.0, # cp (J/(kg*K))
    ),
    region_function=function wall_area!(
        du, u, phys, cell_id, 
        vol,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )
        #property retrieval
        rho = phys.rho

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
    end
)

#we might want to add something akin to distribute_over_set_volume!(du += input_wattage)
total_heating_volume = get_cell_set_total_volume(grid, "heating_areas", geo)

input_wattage = 100.0 # W
corrected_volumetric_heating = input_wattage / total_heating_volume

add_region!(
    config, "heating_areas";
    initial_conditions=ComponentVector(
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics=WallPhysics(
        237.0, # k (W/(m*K))
        2700.0, # rho (kg/m^3)
        921.0, # cp (J/(kg*K))
    ),
    region_function=function heating_area!(
        du, u, phys, cell_id, 
        vol,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )
        #property retrieval
        rho = phys.rho

        #internal physics

        #sources
        corrected_volumetric_heating = 1.2351610076363146e7
        du.temp[cell_id] += corrected_volumetric_heating * vol

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
    end
)

#connection functions
function fluid_fluid_flux!(
        du, u, phys_a, phys_b,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )
    mw_avg!(u, idx_a, phys_a.species_molecular_weights, mw_avg_cache) 
    rho_ideal!(u, idx_a, rho_cache, mw_avg_cache)

    mw_avg!(u, idx_b, phys_b.species_molecular_weights, mw_avg_cache) 
    rho_ideal!(u, idx_b, rho_cache, mw_avg_cache) 

    #mutating-ish, it mutates du.pressure for a
    face_m_dot = continuity_and_momentum_darcy(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        #rho_a, rho_b,
        rho_cache[idx_a], rho_cache[idx_b],
        phys_a, phys_b
    )

    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        phys_a, phys_b
    )

    all_species_advection!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        phys_a, phys_b,
        face_m_dot
    )

    enthalpy_advection!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        phys_a, phys_b,
        face_m_dot
    )

    diffusion_mass_fraction_exchange!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        #rho_a, rho_b,
        rho_cache[idx_a], rho_cache[idx_b],
        phys_a, phys_b
    )
end

function solid_solid_flux!(
        du, u, phys_a, phys_b,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        phys_a, phys_b
    )
end

function fluid_solid_flux!(
        du, u, phys_a, phys_b,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )

    #mw_avg!(u, idx_a, phys_a.species_molecular_weights, mw_avg_cache)
    #rho_ideal!(u, idx_a, rho_cache, mw_avg_cache)
    #rho_b = phys_b.rho

    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        phys_a, phys_b
    )
end

function solid_fluid_flux!(
        du, u, phys_a, phys_b,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )

    #rho_a = phys_a.rho
    #mw_avg!(u, idx_b, phys_b.species_molecular_weights, mw_avg_cache)
    #rho_ideal!(u, idx_b, rho_cache, mw_avg_cache)

    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx],
        phys_a, phys_b
    )
end

#this is the smallest I could make this function
#the reason we use typeof(phys_a) <: AbstractFluidPhysics is just because I think the syntax highlighting makes it look better than phys_a isa AbstractFluidPhysics
#furthermore, since this is a user defined function, you can use isa if you want
function conneciton_map_function(phys_a, phys_b)
    typeof(phys_a) <: AbstractFluidPhysics && typeof(phys_b) <: AbstractFluidPhysics && return fluid_fluid_flux!
    typeof(phys_a) <: AbstractSolidPhysics && typeof(phys_b) <: AbstractSolidPhysics && return solid_solid_flux!
    typeof(phys_a) <: AbstractFluidPhysics && typeof(phys_b) <: AbstractSolidPhysics && return fluid_solid_flux!
    typeof(phys_a) <: AbstractSolidPhysics && typeof(phys_b) <: AbstractFluidPhysics && return solid_fluid_flux!
end

du0, u0, geo, system = finish_fvm_config(config, flux_functions, conneciton_map_function)

N::Int = ForwardDiff.pickchunksize(length(u0))

rho_cache = DiffCache(zeros(length(grid.cells)), N)
mw_avg_cache = DiffCache(zeros(length(grid.cells)), N)

change_in_molar_concentrations_cache = DiffCache(zeros(length(u_proto.mass_fractions[:, 1])), N)
molar_concentrations_cache = DiffCache(zeros(length(u_proto.mass_fractions[:, 1])), N) #just using mass fractions for cell 1, this may cause some issues later!
net_rates_cache = DiffCache(zeros(length(reforming_area_physics.chemical_reactions)), N)

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t,
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,
    system.connection_groups, system.phys, system.cell_phys_id_map,
    system.regions_phys_func_cells,
    system.u_merged_axes,
    rho_cache, mw_avg_cache,
    change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
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
#11721 (10%)

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