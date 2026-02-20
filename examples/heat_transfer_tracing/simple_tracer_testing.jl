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

grid_dimensions = (5, 5, 5)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Tetrahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

n_cells = length(grid.cells)

config = create_fvm_config(grid)

n_faces = length(config.geo.cell_neighbor_areas[1])

add_region!(
    config, "copper";
    type=Solid(),
    region_function=
    function heat_transfer!(du, u, cell_id, cell_volumes)
        #property updating/retrieval

        #variable summation

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
    end
)

add_region!(
    config, "steel";
    type=Solid(),
    region_function=
    function heat_transfer!(du, u, cell_id, cell_volumes)
        #property updating/retrieval

        #variable summation

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
    end
)

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

function conneciton_map_function(type_a, type_b)
    type_a <: Solid && type_b <: Solid && return solid_solid_flux!
end

geo, system = finish_fvm_config(config, conneciton_map_function)

n_faces = length(geo.cell_neighbor_areas[1])

ctx_loop_context = Dict()

du_tracer, u_tracer, ctx = create_tracer_context(ctx_loop_context)

for conn in system.connection_groups
    idx_a = 1
    neighbor_list = 1
    idx_b = 1
    face_idx = 1

    i = 1

    while neighbor_list == Tuple{Int64,Int64}[]
        idx_a = conn.cell_neighbors[i][1]
        neighbor_list = conn.cell_neighbors[i][2]
        idx_b = neighbor_list[1][1]
        face_idx = neighbor_list[1][2]
        i += 1
    end

    fake_cell_neighbor_areas = Dict(
        conn.name_a => Dict(
            "face_idx_1" => geo.cell_neighbor_areas[idx_a][face_idx],
            "face_idx_2" => geo.cell_neighbor_areas[idx_b][face_idx],
            "face_idx_3" => geo.cell_neighbor_areas[idx_a][face_idx],
            "face_idx_4" => geo.cell_neighbor_areas[idx_b][face_idx],
            "face_idx_5" => geo.cell_neighbor_areas[idx_a][face_idx],
            "face_idx_6" => geo.cell_neighbor_areas[idx_b][face_idx]
        )
    )

    fake_cell_neighbor_normals = Dict(
        conn.name_a => Dict(
            "face_idx_1" => geo.cell_neighbor_normals[idx_a][face_idx],
            "face_idx_2" => geo.cell_neighbor_normals[idx_b][face_idx],
            "face_idx_3" => geo.cell_neighbor_normals[idx_a][face_idx],
            "face_idx_4" => geo.cell_neighbor_normals[idx_b][face_idx],
            "face_idx_5" => geo.cell_neighbor_normals[idx_a][face_idx],
            "face_idx_6" => geo.cell_neighbor_normals[idx_b][face_idx]
        )
    )

    fake_cell_neighbor_distances = Dict(
        conn.name_a => Dict(
            "face_idx_1" => geo.cell_neighbor_distances[idx_a][face_idx],
            "face_idx_2" => geo.cell_neighbor_distances[idx_b][face_idx],
            "face_idx_3" => geo.cell_neighbor_distances[idx_a][face_idx],
            "face_idx_4" => geo.cell_neighbor_distances[idx_b][face_idx],
            "face_idx_5" => geo.cell_neighbor_distances[idx_a][face_idx],
            "face_idx_6" => geo.cell_neighbor_distances[idx_b][face_idx]
        )
    )

    conn.flux_function!(
        du_tracer, u_tracer, conn.name_a, conn.name_b, "face_idx_$face_idx", #we swap out idx_a and idx_b for conn.name_a and conn.name_b for the tracer
        fake_cell_neighbor_areas, fake_cell_neighbor_normals, fake_cell_neighbor_distances
    )
end

#Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
#oh wait, now we don't even need the other fields for the different regions, we only need the region function
for reg in system.region_groups
    fake_cell_volumes = Dict(reg.name => geo.cell_volumes[reg.region_cells[1]])
    reg.region_function!(
        du_tracer, u_tracer, reg.name, #we swap out cell_id for reg.name for the tracer
        fake_cell_volumes
    )
end

var_access_logs, encountered_paths = merge_trace_results(ctx.access_logs)

state_vars, cache_vars, fixed_vars = classify_variables(var_access_logs)

n_cells = length(grid.cells)

region_symbols = Set([Symbol(reg.name) for reg in system.region_groups])
controller_symbols = Set([Symbol(cont.name) for cont in system.controller_groups])

n_faces = length(geo.cell_neighbor_areas[1])

state_setup = region_setup(state_vars, region_symbols, n_cells, controller_symbols, n_faces)
fixed_setup = region_setup(fixed_vars, region_symbols, n_cells, controller_symbols, n_faces)

complete_state = build_component_array_merge_regions(state_vars, region_symbols, n_cells, n_faces)
complete_fixed = build_component_array_merge_regions(fixed_vars, region_symbols, n_cells, n_faces)

setup_script_path = joinpath(dirname(@__FILE__), "heat_transfer_setup.jl")
#generate_setup_script(setup_script_path, state_setup, fixed_setup, system)

#rm(setup_script_path, force=true)

include(setup_script_path)

populate_merged_vector!(complete_state, initial_state, system)

populate_merged_vector!(complete_fixed, fixed_properties, system)

N::Int = ForwardDiff.pickchunksize(length(complete_state))
cache_cv = build_component_array_merge_regions(cache_vars, region_symbols, n_cells, n_faces)

complete_fixed = ComponentVector(complete_fixed; NamedTuple(cache_cv)...)

du_fixed = copy(complete_fixed)
du_fixed .= 0.0
u_fixed = copy(complete_fixed)

fixed_axes = getaxes(complete_fixed)[1]

du_fixed_cache = DiffCache(Vector(du_fixed), N)
u_fixed_cache = DiffCache(Vector(u_fixed), N)

u0_flat = Vector(complete_state)
du0_flat = copy(u0_flat)
du0_flat .= 0.0

state_axes = getaxes(complete_state)[1]

all_vars = ComponentVector(complete_state; NamedTuple(complete_fixed)...)
u_all_vars_vec = Vector(all_vars)
du_all_vars_vec = copy(u_all_vars_vec)
all_axes = getaxes(all_vars)[1]

f_closure = (du, u, p, t) -> heat_transfer_f!(
    du, u, p, t,
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,
    system.connection_groups, system.controller_groups, system.region_groups,
    du_fixed_cache, u_fixed_cache,
    state_axes, fixed_axes
)

p_guess = 0.0

test_prob = ODEProblem(f_closure, u0_flat, (0.0, 100.0), p_guess)
@time sol = solve(test_prob, Tsit5(), callback = approximate_time_to_finish_cb)
#=
detector = SparseConnectivityTracer.TracerSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess, 0.0), u0_flat, du0_flat, detector)

ode_func = ODEFunction(f_closure, jac_prototype=float.(jac_sparsity))

t0 = 0.0
tMax = 10.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_flat, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback = approximate_time_to_finish_cb)

record_sol = true

sim_file = @__FILE__

function rebuild_u_named(u_flat, state_axes)
    return [ComponentVector(u_flat[i], state_axes) for i in eachindex(sol.u)]
end

u_named = rebuild_u_named(sol.u, state_axes)

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end
