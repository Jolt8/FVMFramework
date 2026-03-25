using Revise
using Logging

using FVMFramework

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
u_proto = ComponentVector(
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto)

n_faces = length(config.geo.cell_neighbor_areas[1])

struct Solid <: AbstractPhysics end

add_region!(
    config, "copper";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 270.0u"°C",
    ),
    properties = ComponentVector(
        k = 237.0u"W/(m*K)", 
        rho = 2700.0u"kg/m^3",
        cp = 921.0u"J/(kg*K)",
    ),
    optimized_syms = (),
    cache_syms_and_units = (heat = u"J",),
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
    initial_conditions = ComponentVector(
        temp = 350.0u"°C",
    ),
    properties = ComponentVector(
        k = 123.0u"W/(m*K)",
        rho = 7800.0u"kg/m^3",
        cp = 450.0u"J/(kg*K)",
    ),
    cache_syms_and_units = (heat = u"J",),
    optimized_syms = (),
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
    temp_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Solid && typeof(type_b) <: Solid && return solid_solid_flux!
end

species_caches = ()

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, species_caches, check_units = false);

f_closure_implicit = (du, u, p, t) -> heat_transfer_f_test!(
    du, u, p, t, 

    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,

    system.connection_groups, system.controller_groups, system.region_groups, system.patch_groups,

    system.du_virtual_axes, system.u_virtual_axes,
    system.du_diff_cache, system.u_diff_cache,
    system.merged_properties
)

p_guess = 0.0

test_prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1000.0), p_guess)
@btime sol = solve(test_prob, Tsit5(), tspan = (0.0, 10.0))

# 800.793 ms (1552 allocations: 55.95 MiB) (multithreaded)
# 769.695 ms (1552 allocations: 55.95 MiB) (non-multithreaded)

sol.u[1] == sol.u[end]

#using BenchmarkTools
#@btime sol = solve(test_prob, Tsit5(), tspan = (0.0, 10.0))

#VSCodeServer.@profview sol = solve(test_prob, Tsit5(), tspan=(0.0, 100.0), callback=approximate_time_to_finish_cb)

sol_u_named_0 = ComponentVector(sol.u[1], system.u_proto_axes)

sol_u_named_end = ComponentVector(sol.u[end], system.u_proto_axes)

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

@btime sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
#728.605 ms (341178 allocations: 1.08 GiB) (non-multithreaded)
#787.708 ms (202070 allocations: 1.07 GiB) (only connections multithreading)
#799.862 ms (219386 allocations: 1.07 GiB) (everything multithreaded)
#864.517 ms (211023 allocations: 1.07 GiB) (only regions multithreading)


#VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

u_named = rebuild_u_named(sol.u, ComponentVector(u_proto))

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end