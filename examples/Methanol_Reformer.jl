using Revise

using FVMFramework
#I think this take a long time because Enzyme (despite not being installed) throws a warning from DiffEqBaseEnzyme.jl

using Ferrite
using DifferentialEquations
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays

grid_dimensions = (2, 2, 2)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

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

reaction_physics = MethanolReformerPhysics(
    237.0, # k (W/mK)
    1000.0, # rho (kg/m^3)
    1e-5, # mu (Pa*s)
    4.184, # cp (J/kgK)
    0.6e-3, # permeability (m^2)
    species_molecular_weights,
    [MSR_rxn, MD_rxn, WGS_rxn], # chemical_reactions
    [1250, 1250, 1250], # cell_kg_cat_per_m3_for_each_reaction
    [0.0, 0.0, 0.0, 0.0, 0.0], # chemical_vol_source_term
    0.0 # heat_vol_source_term
)

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

addcellset!(grid, "internal_cells", x -> x != -1e9) # all

using Unitful
add_region!(
    config, "internal_cells";
    physics=reaction_physics,
    initial_conditions=ComponentVector(
        vel_x=0.0, vel_y=0.0, vel_z=0.0,
        pressure=ustrip(1.0u"atm" |> u"Pa"),
        mass_fractions=initial_mass_fractions,
        temp=ustrip(270.0u"°C" |> u"K")
    )
)

add_boundary!(
    config, "internal_cells";
    fixed_conditions=ComponentVector(
    #nothing here means everything is free
    )
)

du0, u0, geo, system = finish_fvm_config(config)

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f!(
    du, u, p, t,
    geo.cell_neighbor_map,
    geo.cell_volumes, geo.cell_centroids,
    geo.connection_areas, geo.connection_normals, geo.connection_distances,
    geo.unconnected_areas,
    system.phys, system.cell_phys_id_map, system.fixed_idxs_and_vals_map,
    system.u_axes
)
#just re-add t to the FVM_iter_f! function above to make it compatible with implicit solving

t0, tMax = 0.0, 100000.0
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

#test_prob = ODEProblem(f_closure_implicit, u0, tspan, p_guess)
#@time sol = solve(test_prob, Tsit5())

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true))
VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true))

rebuild_u_named(sol.u, u_proto)[1].mass_fractions

record_sol = false

sim_file = @__FILE__

u_named = rebuild_u_named_vel(sol.u, u_proto)

typeof(u_named[1].velocity)
#the velocities are constricted into a 3 long vec for each cell

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end