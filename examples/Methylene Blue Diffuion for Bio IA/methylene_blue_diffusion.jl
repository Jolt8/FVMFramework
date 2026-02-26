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

mesh_path = joinpath(@__DIR__, "clamped_end/dialysis_tubing_clamped_output.msh")

grid = togrid(mesh_path)

grid.cellsets
grid.facetsets

n_cells = length(grid.cells)
u_proto = (
    mass_fractions = (methylene_blue = zeros(n_cells), water = zeros(n_cells)),
)

config = create_fvm_config(grid, u_proto)

#=
for (cell_id, facet_idx) in grid.facetsets["dialysis_tubing_surface"]
    println(cell_id, " ", facet_idx)
end
=#

check_cellset_connectivity(config.grid, "dialysis_tubing_interior")



#since each mass fraction is modified, they have to be vectors


function normalize_mass_fractions(mass_fractions)
    total_mass_fractions = 0.0

    for (species_name, mass_fraction) in pairs(mass_fractions)
        total_mass_fractions += mass_fractions[species_name][1]
    end

    for species_name in keys(mass_fractions)
        mass_fractions[species_name][1] /= total_mass_fractions
    end

    return NamedTuple{keys(mass_fractions)}(first.(values(mass_fractions)))
end

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

dialysis_tubing_initial_mass_fractions = (
    methylene_blue = [0.0003846],
    water = [0.9996154]
)

dialysis_tubing_initial_mass_fractions = normalize_mass_fractions(dialysis_tubing_initial_mass_fractions)

add_region!(
    config, "dialysis_tubing_interior";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = dialysis_tubing_initial_mass_fractions,
    ),
    properties = (
        temp = ustrip(21u"°C" |> u"K"),
        k = 0.6, # k (W/(m*K))
        cp = 4184, # cp (J/(kg*K))
        rho = 1000, # rho (kg/m^3)
        mu = 1e-3, # mu (Pa*s)
        species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
        diffusion_coefficients = (
            methylene_blue = 1e-5,
            water = 1e-5
        ), #diffusion coefficients (m^2/s)
        molecular_weights = (
            methylene_blue = 0.31985, 
            water = 0.01802
        ), #species_molecular_weights [kg/mol]
    ), 
    state_syms = [:mass_fractions],
    cache_syms = [:heat, :molar_concentrations], 
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        molar_concentrations!(u, cell_id)

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

surrounding_fluid_initial_mass_fractions = (
    methylene_blue = [0.0],
    water = [1.0]
)

surrounding_fluid_initial_mass_fractions = normalize_mass_fractions(surrounding_fluid_initial_mass_fractions)

add_region!(
    config, "surrounding_fluid";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = surrounding_fluid_initial_mass_fractions,
    ),
    properties = (
        temp = ustrip(21u"°C" |> u"K"),
        k = 0.6, # k (W/(m*K))
        cp = 4184, # cp (J/(kg*K))
        rho = 1000, # rho (kg/m^3)
        mu = 1e-3, # mu (Pa*s)
        species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
        diffusion_coefficients = (
            methylene_blue = 1e-5,
            water = 1e-5
        ), #diffusion coefficients (m^2/s)
        molecular_weights = (
            methylene_blue = 0.31985, 
            water = 0.01802
        ), #species_molecular_weights [kg/mol]
    ), 
    state_syms = [:mass_fractions],
    cache_syms = [:heat, :molar_concentrations], 
    region_function =
    function surrounding_fluid!(du, u, cell_id, vol)
        #property updating/retrieval
        molar_concentrations!(u, cell_id)

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    molar_concentrations!(u, idx_a)

    molar_concentrations!(u, idx_b)

    #=
    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    diffusion_mass_fraction_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

#this is the smallest I could make this function
#I like using <: here because it makes it look nice with syntax highlighting
function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
)

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, special_caches)

u_test = ComponentVector(u0_vec, system.u_proto_axes)

u_nt = NamedTuple(u_test)

create_axes(u_nt, n_cells)

#=

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t, geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals, system.connection_groups, system.controller_groups, system.region_groups,
    system.merged_properties, system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

p_guess = 0.0

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 100000.0), p_guess)
@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

sol_u_named_0 = ComponentVector(sol.u[1], system.u_proto_axes)

sol_u_named_end = ComponentVector(sol.u[end], system.u_proto_axes)

sol_u_named_0.mass_fractions.methylene_blue == sol_u_named_end.mass_fractions.methylene_blue

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

diff_cache = get_tmp(system.u_diff_cache_vec, 1.0)

test_1 = ComponentVector(diff_cache, system.u_cache_axes)

test_1.molar_concentrations[:methylene_blue][1] = 1.0
test_1.molar_concentrations[:methylene_blue][1]

test_2 = (; ComponentVector(diff_cache, system.u_cache_axes)...)

test_2.molar_concentrations[:methylene_blue][1] = 1.0
test_2.molar_concentrations[:methylene_blue][1]

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

t0 = 0.0
tMax = 1000000.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

u_proto_named = [(; ComponentVector(sol.u[i], system.u_proto_axes)...) for i in eachindex(sol.u)]

#u_named = rebuild_u_named(sol.u, u_proto_named)

cell = grid.cellsets["dialysis_tubing_interior"][1]

mass_fractions_dialysis_tubing_interior = []

for step in eachindex(sol.u)
    push!(mass_fractions_dialysis_tubing_interior, u_proto_named[step].mass_fractions.methylene_blue[cell])
end

mass_fractions_dialysis_tubing_interior

sol.u[1] == sol.u[end]

mass_fractions_beginning = u_proto_named[1].mass_fractions.methylene_blue[cell]

mass_fractions_end = u_proto_named[end].mass_fractions.methylene_blue[cell]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

if record_sol == true
    sol_to_vtk(sol, u_proto_named, grid, sim_file)
end