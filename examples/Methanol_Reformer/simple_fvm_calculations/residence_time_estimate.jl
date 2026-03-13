using Revise

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

using Unitful

grid_dimensions = (2, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Tetrahedron, grid_dimensions, left, right)

addcellset!(grid, "reforming_area", x -> true == true)

grid.cellsets

n_cells = length(grid.cells)
u_proto = (
    mass_fractions = (methanol = zeros(n_cells), water = zeros(n_cells), carbon_monoxide = zeros(n_cells), hydrogen = zeros(n_cells), carbon_dioxide = zeros(n_cells)),
)

config = create_fvm_config(grid, u_proto)

#for CH3OH, HCOO, OH
van_t_hoff_A = (CH3O = 1.7e-6, HCOO = 4.74e-13, OH = 3.32e-14)
van_t_hoff_dH = (CH3O = -46800.0, HCOO = -115000.0, OH = -110000.0)

MSR_rxn = (
    heat_of_reaction = 49500.0, # [J/mol]
    ref_delta_G = -3800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 1.25e7, # [s^-1]
    kf_Ea = 103000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = (methanol = -1, water = -1, carbon_monoxide = 0, hydrogen = 3, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

MD_rxn = (
    heat_of_reaction = 90200.0, # [J/mol]
    ref_delta_G = 24800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 1.15e11, # [s^-1]
    kf_Ea = 170000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1,), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), # product_stoich_coeffs: 1 CO + 2 H2
    stoich_coeffs = (methanol = -1, water = 0, carbon_monoxide = 1, hydrogen = 2, carbon_dioxide = 0), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

WGS_rxn = (
    heat_of_reaction = -41100.0, # [J/mol]
    ref_delta_G = -28600.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 3.65e7, # [s^-1]
    kf_Ea = 87500.0, # [J/mol]
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs = (methanol = 0, water = -1, carbon_monoxide = -1, hydrogen = 1, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

#since each mass fraction is modified, they have to be vectors
initial_mass_fractions = (
    methanol = [1.0],
    water = [1.3],
    carbon_monoxide = [0.0001],
    hydrogen = [0.02],
    carbon_dioxide = [0.0001]
)

total_mass_fractions = 0.0

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    total_mass_fractions += initial_mass_fractions[species_name][1]
end

total_mass_fractions

#each mass fraction must be a vector to be mutable
initial_mass_fractions = (
    methanol = [1.0],
    water = [1.3],
    carbon_monoxide = [0.0001],
    hydrogen = [0.02],
    carbon_dioxide = [0.0001]
)

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    initial_mass_fractions[species_name][1] /= total_mass_fractions
end

initial_mass_fractions = NamedTuple{keys(initial_mass_fractions)}(first.(values(initial_mass_fractions)))
#this is not ideal, but we can't use scalars inside the initial mass_fractions because the values of NamedTuples can't be modified 

reforming_area_properties = (
    k = 237.0, # k (W/(m*K))
    cp = 4.184, # cp (J/(kg*K))
    mu = 1e-5, # mu (Pa*s)
    rho = 791.0, # rho (kg/m^3)
    temp = ustrip(270.0u"°C" |> u"K"),
    pressure = ustrip(1.0u"atm" |> u"Pa"),
    species_ids = (methanol = 1, water = 2, carbon_monoxide = 3, hydrogen = 4, carbon_dioxide = 5), #we could use mass_fractions for species loops, but this is just more consistent
    molecular_weights = (
        methanol = 0.03204,
        water = 0.01802,
        carbon_monoxide = 0.02801,
        hydrogen = 0.00202,
        carbon_dioxide = 0.04401
    ), #species_molecular_weights [kg/mol]
    reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
    reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0, MD_rxn = 1250.0, WGS_rxn = 1250.0),), # cell_kg_cat_per_m3_for_each_reaction
)
#remember, for any named tuple with a single field inside it, remember to add a comma to the end 
#ex. (reforming_reactions = (MSR_rxn = MSR_rxn, ),)
#not: (reforming_reaction = (MSR_rxn = MSR_rxn))

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = initial_mass_fractions,
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    state_syms = [:mass_fractions],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :net_rates],
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics

        PAM_reforming_react_cell!(du, u, cell_id, vol)
    end
)
#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
n_reactions = length(config.regions[1].properties.reactions.reforming_reactions)
reaction_names = keys(config.regions[1].properties.reactions.reforming_reactions)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    net_rates = (reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),), 
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
)

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, special_caches)

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t, 

    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,

    system.connection_groups, system.controller_groups, system.region_groups, system.patch_groups,
    system.merged_properties,

    system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)

p_guess = 0.0

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 400000.0), p_guess)
@time sol = solve(prob, Tsit5())

u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)
u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

methanol_initial_mass_fraction = u_named_0.mass_fractions.methanol[1]
methanol_final_mass_fraction = u_named_end.mass_fractions.methanol[1]

conversion = (methanol_initial_mass_fraction - methanol_final_mass_fraction) / methanol_initial_mass_fraction

#400000 seconds for 90% conversion, this needs to be fixed 

hi = 1

t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)=
VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb, saveat = tMax/100)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

u_proto_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

#u_named = rebuild_u_named(sol.u, u_proto_named)

cell = grid.cellsets["inlet"][1]

inlet_pressure = []

for step in eachindex(sol.u)
    push!(inlet_pressure, u_proto_named[step].pressure[cell])
end

inlet_pressure

sol.u[1] == sol.u[end]

mass_fractions_beginning = u_proto_named[1].mass_fractions.methanol[6389]

mass_fractions_end = u_proto_named[end].mass_fractions.methanol[6389]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

length(u_proto_named[1])

if record_sol == true
    sol_to_vtk(sol, u_proto_named, grid, sim_file)
end