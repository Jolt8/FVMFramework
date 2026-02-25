using Revise
using Logging

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

grid_dimensions = (100, 10, 10)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Tetrahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

n_cells = length(grid.cells)
u_proto = (
    temp = zeros(n_cells),
)

config = create_fvm_config(grid, u_proto)

n_faces = length(config.geo.cell_neighbor_areas[1])

struct Solid <: AbstractPhysics end

add_region!(
    config, "copper";
    type = Solid(),
    initial_conditions = (
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = (
        k = 237.0, # k (W/(m*K))
        rho = 2700.0, # rho (kg/m^3)
        cp = 921.0, # cp (J/(kg*K))
    ),
    state_syms = [:temp],
    cache_syms = [:heat],
    region_function =
    function heat_transfer!(du, u, cell_id, vol)
        #property updating/retrieval

        #variable summation

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "steel";
    type = Solid(),
    initial_conditions = (
        temp = ustrip(350.0u"°C" |> u"K"), #I really dislike that NamedTuples default to vectors if you don't put , at the end for a single field NamedTuples
    ),
    properties = (
        k = 123.0, # k (W/(m*K))
        rho = 7800.0, # rho (kg/m^3)
        cp = 450.0, # cp (J/(kg*K))
    ),
    state_syms = [:temp],
    cache_syms = [:heat],
    region_function=
    function heat_transfer!(du, u, cell_id, vol)
        #property updating/retrieval

        #variable summation

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

function solid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Solid && typeof(type_b) <: Solid && return solid_solid_flux!
end

species_caches = ()

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, species_caches)

f_closure_implicit = (du, u, p, t) -> heat_transfer_f_test!(
    du, u, p, t, 

    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,

    system.connection_groups, system.controller_groups, system.region_groups, 
    system.merged_properties,

    system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

p_guess = 0.0

test_prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1000.0), p_guess)
@time sol = solve(test_prob, Tsit5(), tspan=(0.0, 10.0), callback=approximate_time_to_finish_cb)

sol.u[1] == sol.u[end]

#using BenchmarkTools
#@btime sol = solve(test_prob, Tsit5(), tspan = (0.0, 10.0))

VSCodeServer.@profview sol = solve(test_prob, Tsit5(), tspan=(0.0, 10.0), callback=approximate_time_to_finish_cb)

sol_u_named_0 = ComponentVector(sol.u[1], system.u_proto_axes)

sol_u_named_end = ComponentVector(sol.u[end], system.u_proto_axes)

t0 = 0.0
tMax = 10000.0
tspan = (t0, tMax)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback = approximate_time_to_finish_cb)
#VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback = approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

u_named = rebuild_u_named(sol.u, ComponentVector(u_proto))

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end