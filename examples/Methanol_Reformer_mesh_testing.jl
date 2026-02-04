using Revise

using FVMFramework
#I think this take a long time because Enzyme (despite not being installed) throws a warning from DiffEqBaseEnzyme.jl

using Ferrite
using FerriteGmsh
using DifferentialEquations
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays

using Unitful

grid = togrid("C://Users//wille//Desktop//FreeCad Projects//Methanol Reformer//MeOH for gmsh.msh")

grid.cells[1]

grid.cellsets
grid.facetsets

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
    0.6e-3, # permeability (m^2)
    species_molecular_weights,
    [MSR_rxn, MD_rxn, WGS_rxn], # chemical_reactions
    [1250.0, 1250.0, 1250.0], # cell_kg_cat_per_m3_for_each_reaction
)

#=
add_default_region!(
    config, "default";
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=0.0, vel_z=0.0,
        pressure=ustrip(1.0u"atm" |> u"Pa"),
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics = reforming_area_physics,
    region_function=function reforming_area!(
            du, u, cell_id, 
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, rho, phys
        )
        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
    end
)
=#

add_region!(
    config, "reforming_area";
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=0.0, vel_z=0.0,
        pressure=ustrip(1.0u"atm" |> u"Pa"),
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics = reforming_area_physics,
    region_function=function reforming_area!(
            du, u, cell_id, 
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, rho, phys
        )
        #internal physics
        react_cell!(
            cell_id, du, u,
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, rho,
            phys,
        )

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol, rho, phys)
    end
)


desired_m_dot = 0.2u"g/s" |> u"kg/s"

#we should probably create methods to automatically copy the parameters from other already defined regions
#=
add_facet_region!(
    config, "inlet",
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=1.0, vel_z=0.0,
        pressure=ustrip(1.1u"atm" |> u"Pa"),
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    ),
    region_physics=reforming_area_physics,
    region_function=function inlet_area!(
            du, u, cell_id, 
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, rho, phys
        )
        #internal physics
        react_cell!(
            cell_id, du, u,
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, rho,
            phys,
        )

        #sources
        du.pressure[cell_id] += desired_m_dot

        #boundary conditions
        du.vel_x[cell_id] = 0.0
        u.vel_x[cell_id] = 0.0

        du.vel_y[cell_id] = 0.0
        u.vel_y[cell_id] = 1.0

        du.vel_z[cell_id] = 0.0
        u.vel_z[cell_id] = 0.0

        du.mass_fractions .= 0.0
        u.mass_fractions = initial_mass_fractions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol, rho, phys)
    end
)
=#


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
    region_function=function wall_area!(du, u, cell_id, vol, rho, phys)
    #internal physics
    #no internal physics in this region
    
    #sources
    #no sources in this region

    #boundary conditions
    #no boundary conditions in this region

    #capacities
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
end
)

#we might want to automate this in the future
geo = build_fvm_geo_into_struct(grid)

total_heating_volume = 0.0

for cell_id in eachindex(geo.cell_volumes)
    if cell_id in getcellset(grid, "heating_areas")
        total_heating_volume += geo.cell_volumes[cell_id]
    end
end

total_heating_volume

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
    region_function=function wall_area!(
            du, u, cell_id, 
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            vol, rho, 
            phys
        )
        #internal physics
        #no internal physics in this region
        
        #sources
        du.temp[cell_id] += vol * corrected_volumetric_heating

        #boundary conditions
        #no boundary conditions in this region

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol, rho, phys)
    end
)

connection_groups = methanol_reformer_init_conn_groups()

du0, u0, geo, system = finish_fvm_config(config, connection_catagorizer!, connection_groups)

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t,
    geo.cell_neighbor_map,
    geo.cell_volumes, geo.cell_centroids,
    geo.connection_areas, geo.connection_normals, geo.connection_distances,
    geo.unconnected_areas,

    system.connection_groups, system.phys, system.cell_phys_id_map,
    system.regions_phys_func_cells,
    system.u_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

t0, tMax = 0.0, 1000000.0
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

#test_prob = ODEProblem(f_closure_implicit, u0, tspan, p_guess)
#@time sol = solve(test_prob, Tsit5())

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

#test_prob = ODEProblem(f_closure_implicit, u0, (0.0, 0.0000001), p_guess)

#solve(test_prob, Tsit5())

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

#@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true))
VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true))

rebuild_u_named(sol.u, u_proto)[60].mass_fractions

record_sol = false

sim_file = @__FILE__

u_named = rebuild_u_named_vel(sol.u, u_proto)

sum(u_named[95].mass_fractions[:, 1])

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end