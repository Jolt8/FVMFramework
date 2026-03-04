using FVMFramework
using Unitful
using OrdinaryDiffEq
using Ferrite
using SparseConnectivityTracer
import ADTypes

total_pipe_length = 50.0u"cm" |> u"m"
stripped_pipe_length = ustrip(total_pipe_length)
pipe_width = ustrip(1.0u"cm" |> u"m")

grid_dimensions = (100, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((stripped_pipe_length, pipe_width, pipe_width))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

total_pipe_segments = 100

addcellset!(grid, "vaporization_area", xyz -> xyz[1] <= (20 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "vaporization_area")
addcellset!(grid, "reforming_area", xyz -> xyz[1] >= (20 * (stripped_pipe_length / total_pipe_segments)) && xyz[1] <= (100 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "reforming_area")

n_cells = total_pipe_segments
u_proto = (
    mass_fractions = (
        methanol = zeros(n_cells),
        water = zeros(n_cells),
        carbon_monoxide = zeros(n_cells),
        hydrogen = zeros(n_cells),
        carbon_dioxide = zeros(n_cells)
    ),
    pressure = zeros(n_cells),
)

config = create_fvm_config(grid, u_proto)

van_t_hoff_A = (CH3O = 1.7e-6, HCOO = 4.74e-13, OH = 3.32e-14)
van_t_hoff_dH = (CH3O = -46800.0, HCOO = -115000.0, OH = -110000.0)

function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    K_ref = exp(-ΔG_rxn_ref / (R_gas * T_ref))

    ln_K_ratio = (-ΔH_rxn_ref / R_gas) * (1/T_actual - 1/T_ref)
    
    K_T = K_ref * exp(ln_K_ratio)
    
    return ustrip(K_T) # Return a unitless equilibrium constant
end

MSR_rxn = (
    heat_of_reaction = 49500.0, # [J/mol]
    ref_delta_G = -3800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 1.25e12,#1.25e7, # [s^-1] #sources online point to values around 1.25e7 mol / (kg * s * bar)
    kf_Ea = 103000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = (methanol = -1, water = -1, carbon_monoxide = 0, hydrogen = 3, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

K_gibbs_free(MSR_rxn.ref_temp, ustrip(270u"°C" |> u"K"), MSR_rxn.ref_delta_G, MSR_rxn.heat_of_reaction)

MD_rxn = (
    heat_of_reaction = 90200.0, # [J/mol]
    ref_delta_G = 24800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 1.15e16,#1.15e11, # [s^-1] #sources online point to values around 1.15e11 mol / (kg * s * bar)
    kf_Ea = 170000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1,), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), # product_stoich_coeffs: 1 CO + 2 H2
    stoich_coeffs = (methanol = -1, water = 0, carbon_monoxide = 1, hydrogen = 2, carbon_dioxide = 0), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

K_gibbs_free(MD_rxn.ref_temp, ustrip(270u"°C" |> u"K"), MD_rxn.ref_delta_G, MD_rxn.heat_of_reaction)

WGS_rxn = (
    heat_of_reaction = -41100.0, # [J/mol]
    ref_delta_G = -28600.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 3.65e12, #3.65e7, # [s^-1] #sources online point to values around 3.65e7 mol / (kg * s * bar)
    kf_Ea = 87500.0, # [J/mol]
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs = (methanol = 0, water = -1, carbon_monoxide = -1, hydrogen = 1, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

K_gibbs_free(WGS_rxn.ref_temp, ustrip(270u"°C" |> u"K"), WGS_rxn.ref_delta_G, WGS_rxn.heat_of_reaction)

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

total_pipe_length = 50.1u"cm" |> u"m"
n_segments = length(grid.cells)

pipe_mass_flow = 0.0278u"g/s" |> u"kg/s"
pipe_area = 2.0u"cm^2" |> u"m^2"
superficial_mass_velocity = pipe_mass_flow / pipe_area

reforming_area_properties = (
    k = 237.0, # k (W/(m*K))
    cp = 4.184, # cp (J/(kg*K))
    mu = 1e-5, # mu (Pa*s)
    rho = 791.0, # rho (kg/m^3)
    temp = ustrip(270.0u"°C" |> u"K"),
    viscosity = ustrip(1e-5u"Pa*s" |> u"Pa*s"),
    
    pipe_mass_flow = ustrip(pipe_mass_flow |> u"kg/s"),
    pipe_area = ustrip(pipe_area |> u"m^2"),
    pipe_length = ustrip((total_pipe_length / n_segments) |> u"m"),
    superficial_mass_velocity = ustrip(superficial_mass_velocity |> u"kg/(m^2*s)"),
    bed_void_fraction = 0.80,
    catalyst_particle_diameter = ustrip(1.0u"mm" |> u"m"),
    
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

include(joinpath(@__DIR__, "1D physics", "ergun_pressure_drop.jl"))

add_region!(
    config, "vaporization_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = initial_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :net_rates, :mass, :mass_face],
    region_function =
    function vaporization_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        #ergun_pressure_drop!(du, u, cell_id, vol)

        #sum_mass_flux_face_to_cell!(du, u, cell_id)
    end
)

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = initial_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :net_rates, :mass, :mass_face],
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        PAM_reforming_react_cell!(du, u, cell_id, vol)

        #ergun_pressure_drop!(du, u, cell_id, vol)

        #sum_mass_flux_face_to_cell!(du, u, cell_id)
    end
)


#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    #=all_species_advection!(
        du, u, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    #=enthalpy_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
config.regions[1].properties.reactions
n_reactions = length(config.regions[1].properties.reactions.reforming_reactions)
reaction_names = keys(config.regions[1].properties.reactions.reforming_reactions)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))), ), 
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
)

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, special_caches)

system.connection_groups[1].cell_neighbors[end]

f_closure_implicit = (du, u, p, t) -> pipe_f!(
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

function conversion_from_residence_time(tMax)
    prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, tMax), p_guess)
    sol = solve(prob, Tsit5())

    sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)
    sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

    methanol_conversion = ((sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]) / sol_u_named_0.mass_fractions.methanol[80])

    return methanol_conversion
end

using Roots
desired_conversion = 0.95
time = 1.0

find_zero(time -> conversion_from_residence_time(time) - desired_conversion, time)

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 5.0), p_guess)
@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)
sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

sol_u_named_0.pressure
sol_u_named_end.pressure

sol_u_named_0.mass_fractions.methanol
sol_u_named_end.mass_fractions.methanol

sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]

methanol_conversion = ((sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]) / sol_u_named_0.mass_fractions.methanol[80]
)

sol_u_named_0.mass_fractions.hydrogen
sol_u_named_end.mass_fractions.hydrogen

t0 = 0.0
tMax = 500000.0
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

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
VSCodeServer.@profview sol = solve(implicit_prob, FBDF())#, callback = approximate_time_to_finish_cb, saveat = tMax/100)

record_sol = true

sim_file = @__FILE__

u_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

u_named[1].mass_fractions.methanol[80]
u_named[end].mass_fractions.methanol[80]

conversion = ((u_named[1].mass_fractions.methanol[80] - u_named[end].mass_fractions.methanol[80]) / u_named[1].mass_fractions.methanol[80]
)

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end

using Plots

plot(sol.t, [u.mass_fractions.methanol[80] for u in u_named], label = "Methanol")
plot!(sol.t, [u.mass_fractions.water[80] for u in u_named], label = "Water")
plot!(sol.t, [u.mass_fractions.carbon_monoxide[80] for u in u_named], label = "Carbon Monoxide")
plot!(sol.t, [u.mass_fractions.hydrogen[80] for u in u_named], label = "Hydrogen")
plot!(sol.t, [u.mass_fractions.carbon_dioxide[80] for u in u_named], label = "Carbon Dioxide")

methanol_molar_concentration = [(u.mass_fractions.methanol[80] * 1000 / reforming_area_properties.molecular_weights.methanol) for u in u_named]
water_molar_concentration = [(u.mass_fractions.water[80] * 1000 / reforming_area_properties.molecular_weights.water) for u in u_named]
carbon_monoxide_molar_concentration = [(u.mass_fractions.carbon_monoxide[80] * 1000 / reforming_area_properties.molecular_weights.carbon_monoxide) for u in u_named]
hydrogen_molar_concentration = [(u.mass_fractions.hydrogen[80] * 1000 / reforming_area_properties.molecular_weights.hydrogen) for u in u_named]
carbon_dioxide_molar_concentration = [(u.mass_fractions.carbon_dioxide[80] * 1000 / reforming_area_properties.molecular_weights.carbon_dioxide) for u in u_named]

plot(sol.t, methanol_molar_concentration, label = "Methanol")
plot!(sol.t, water_molar_concentration, label = "Water")
plot!(sol.t, carbon_monoxide_molar_concentration, label = "Carbon Monoxide")
plot!(sol.t, hydrogen_molar_concentration, label = "Hydrogen")
plot!(sol.t, carbon_dioxide_molar_concentration, label = "Carbon Dioxide")

total_molar_concentration = methanol_molar_concentration + water_molar_concentration + carbon_monoxide_molar_concentration + hydrogen_molar_concentration + carbon_dioxide_molar_concentration

methanol_molar_fractions = methanol_molar_concentration ./ total_molar_concentration
water_molar_fractions = water_molar_concentration ./ total_molar_concentration
carbon_monoxide_molar_fractions = carbon_monoxide_molar_concentration ./ total_molar_concentration
hydrogen_molar_fractions = hydrogen_molar_concentration ./ total_molar_concentration
carbon_dioxide_molar_fractions = carbon_dioxide_molar_concentration ./ total_molar_concentration

plot(sol.t, methanol_molar_fractions, label = "Methanol")
plot!(sol.t, water_molar_fractions, label = "Water")
plot!(sol.t, carbon_monoxide_molar_fractions, label = "Carbon Monoxide")
plot!(sol.t, hydrogen_molar_fractions, label = "Hydrogen")
plot!(sol.t, carbon_dioxide_molar_fractions, label = "Carbon Dioxide")


1.0e6u"mol/(kg*s*bar)" |> u"mol/(kg*s*Pa)"

1.0e6u"bar/m" |> u"Pa/m"
1.0e6u"m/bar" |> u"m/Pa"