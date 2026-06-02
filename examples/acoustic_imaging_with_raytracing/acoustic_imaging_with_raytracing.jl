using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak
using Dates

using SparseArrays
using JLD2

using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ForwardDiff
using DataInterpolations

using FVMFramework

grid_dimensions = (20, 20, 20)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "fluid", xyz -> true)
getcellset(grid, "fluid")

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

#properties
Revise.includet(joinpath(@__DIR__, "properties", "water_properties.jl"))
water_properties = get_water_properties()

function update_fluid_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)
    
    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.fluid_rho[cell_id]

    rho_effective = (u.gas_holdup[cell_id] * u.gas_density[cell_id]) + (u.liquid_holdup[cell_id] * u.fluid_density[cell_id])
    #I wonder if we would use fluid_rho instead of rho_effective 
    compressibility_effective = u.gas_holdup[cell_id] * properties.gas_compressibility[cell_id] + u.liquid_holdup[cell_id] * properties.liquid_compressibility[cell_id]
    
    u.speed_of_sound[cell_id] = 1 / sqrt(rho_effective * compressibility_effective)

    u.n_velocity_updates[cell_id] += 1.0 #because julia uses 1-based indexing, we have to make n_velocity_updates start at 1
end

function fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

u_proto = ComponentVector(
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])

transducer_node_ids = [1, 10, 100, 200, 400, 500, 600, 3000, 5000, 8000]
transducer_ids = collect(1:length(transducer_node_ids))

#Ray mapping
Revise.includet(joinpath(@__DIR__, "transducer_profiles", "ray_mapping", "ray_tracers", "default_ray_tracer.jl"))
Revise.includet(joinpath(@__DIR__, "transducer_profiles", "ray_mapping", "generate_ray_map.jl"))
#generate_ray_map(transducer_ids, transducer_node_ids) #this takes a while

ray_map_intersected_cells, ray_map_distances_through_cells, ray_map_ray_lengths = load_object(
    joinpath(
        @__DIR__, 
        "ray_mapping",
        "saved_ray_maps", 
        "ray_map_$(length(transducer_ids))_transducers_$(length(grid.cells))_cells.jls"
    )
)

#Projections 
transducer_opposing_pairs_ids = [
    (1, 6),
    (2, 7),
    (3, 8),
    (4, 9),
    (5, 10)
]

Revise.includet(joinpath(@__DIR__, "transducer_profiles", "transducer_projecting", "default_transducer_projector.jl"))
Revise.includet(joinpath(@__DIR__, "transducer_profiles", "transducer_projecting", "generate_transducer_projections.jl"))
generate_transducer_projections(transducer_opposing_pairs_ids, transducer_node_ids, grid, geo, most_consistent_timestamp_for_speed_of_sound, experimental_travel_time_interp_matrix, transducer_frequency, transducer_diameter)

transducer_projection_intersected_cells, transducer_projection_volume_in_cells, transducer_projection_cell_distances, transducer_projection_distances, transducer_projection_slant_distances, cell_projection_counts, projection_unit_vectors_to_cell_centers = load_object(
    joinpath(
        @__DIR__, 
        "transducer_profiles",
        "transducer_projecting",
        "saved_transducer_projections", 
        "transducer_projection_$(length(transducer_opposing_pairs_ids))_pairs_$(length(grid.cells))_cells.jls"
    )
)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)",
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        compressibility = u"Pa^-1",
        speed_of_sound = u"m/s",
        mass = u"kg",
        mass_face = u"kg",
        n_velocity_updates = u"1", #this is for keeping track of how many times a cell has been updated by the experimental velocity data
        beam_x_components = u"1",
        beam_y_components = u"1",
        beam_z_components = u"1",
        beam_measured_velocities = u"m/s",
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",
        beam_x_components = [zeros(Float64, cell_projection_counts[cell_id]) for cell_id in grid.cells]u"1",
        beam_y_components = [zeros(Float64, cell_projection_counts[cell_id]) for cell_id in grid.cells]u"1",
        beam_z_components = [zeros(Float64, cell_projection_counts[cell_id]) for cell_id in grid.cells]u"1",
        beam_measured_velocities = [zeros(Float64, cell_projection_counts[cell_id]) for cell_id in grid.cells]u"m/s",
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector()
)

add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = ComponentVector(
        pressure = 1.0u"atm",
        temp = 25.0u"°C",
    ),
    properties = water_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)
#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    #this is the only equation that's still needed
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function fluid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function solid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return fluid_solid_flux!

    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!
end

function solve_system!(du, u, p_vec, t, geo, system)
    for cell_id in grid.cellsets["fluid"]
        update_fluid_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
    end

    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system)
end
    

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = 0.0

#test_prob = ODEProblem(f_closure, u0_vec, (0.0, 0.01), p_guess)
#@time sol = solve(test_prob, Tsit5(), callback = approximate_time_to_finish_cb)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 2300.0

tspan = (t0, tMax)

saveat = 10.0

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

@time sol = solve(implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

#u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)];

#root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

#sol_to_vtk(sol, u_named, grid, @__FILE__, root_dir)


function build_simulated_travel_time_matrix(tMax, saveat)
    n_timesteps = Int(tMax / saveat)
    simulated_travel_time_matrix = [zeros(Float64, n_timesteps) for i in 1:length(transducer_ids), j in 1:length(transducer_ids)]

    return simulated_travel_time_matrix
end

function calculate_simulated_travel_time(u, origin_node_id, destination_node_id, ray_map_intersected_cells, ray_map_distances_through_cells)
    intersected_cells = ray_map_intersected_cells[origin_node_id, destination_node_id]
    intersected_cell_distances = ray_map_distances_through_cells[origin_node_id, destination_node_id]

    travel_time = 0.0

    for (index, cell_id) in enumerate(intersected_cells)
        travel_time += intersected_cell_distances[index] / u.speed_of_sound[cell_id]
    end

    return travel_time
end

function calculate_simulated_received_returned_volume(u, t, origin_node_id, destination_node_id, ray_map_intersected_cells, ray_map_distances_through_cells)
    intersected_cells = ray_map_intersected_cells[origin_node_id, destination_node_id]
    intersected_cell_distances = ray_map_distances_through_cells[origin_node_id, destination_node_id]

    initial_ping_volume = experimental_ping_volumes_interp_matrix[origin_node_id, destination_node_id](t)

    overall_attenuation_coefficient = 0.0

    #just using a simple attenuiation_coefficient model for now 
    for (index, cell_id) in enumerate(intersected_cells)
        overall_attenuation_coefficient += (u.gas_holdup[cell_id] * u.bubble_radii[cell_id]) * intersected_cell_distances[index]
    end

    overall_attenuation_coefficient /= ray_map_ray_lengths[origin_node_id, destination_node_id]

    returned_volume = initial_ping_volume * exp(-overall_attenuation_coefficient * ray_map_ray_lengths[origin_node_id, destination_node_id])

    return returned_volume
end

Revise.includet(joinpath(@__DIR__, "precompute_experimental_velocities.jl"))
experimental_vel_x_interps, experimental_vel_y_interps, experimental_vel_z_interps, = precompute_experimental_velocities_over_time(
    grid, geo,
    #cell_speeds_of_sounds, 
    transducer_ids,
    transducer_opposing_pairs_ids,
    ray_map_intersected_cells,
    ray_map_distances_through_cells,
    transducer_projection_intersected_cells,
    transducer_projection_cell_distances,
    projection_unit_vectors_to_cell_centers,
    cell_projection_counts,

    experimental_double_pulse_echo_history, #TODO: make sure that the two pulses are done at every desired saveat time

    pulse_repetition_interval, #make sure to update this
    samples_per_second, 
    window_size_indices;
    tolerance = 2.0,
    total_simulation_time, 
    simulation_saveat
)

function built_trial_implicit_prob(f_closure, du0_vec, u0_vec, tMax, p_guess)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end

test_interp = LinearInterpolation(
    [0.0, 1.0, 2.0, 3.0],
    [0.0, 1500.0, 4000.0, 8000.0]
)

experimental_travel_time_interp_matrix = [deepcopy(test_interp) for i in 1:length(transducer_ids), j in 1:length(transducer_ids)]

function loss(θ)
    prob, system_copy, geo_copy = get!(task_local_storage(), :implicit_prob_setup) do
        system_copy = deepcopy(system)
        geo_copy = deepcopy(geo)
        
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo_copy, system_copy)
        
        prob = built_trial_implicit_prob(f_closure, du0_vec, u0_vec, experimental_tMax, θ)
        return (prob, system_copy, geo_copy)
    end

    loss_prob = remake(prob, p = θ)

    sol = solve(
        loss_prob,
        FBDF(linsolve = SparspakFactorization()), 
        sensealg = ForwardSensitivity(),
        saveat = saveat,
        callback = approximate_time_to_finish_cb
    )

    n_saves = Int(experimental_tMax / saveat)

    if length(sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(sol.u[i], state_axes) for i in eachindex(sol.u)]

    du_temporary = similar(sol.u[1])

    mean_squared_error = 0.0

    for i in eachindex(sol.t)
        du, u = unpack_fvm_state(du_temporary, sol.u[i], θ, sol.t[i], system_copy)
        
        for cell_id in grid.cellsets["fluid"]
            update_fluid_properties!(du, u, cell_id, geo_copy.cell_volumes[cell_id], system_copy)
        end

        for origin_transducer_idx in transducer_ids
            for destination_transducer_idx in transducer_ids
                if origin_transducer_idx != destination_transducer_idx
                    experimental_travel_time = experimental_travel_time_interp_matrix[origin_transducer_idx, destination_transducer_idx](sol.t[i])
                    sim_travel_time = calculate_simulated_travel_time(
                        u, origin_transducer_idx, destination_transducer_idx, 
                        ray_map_intersected_cells, ray_map_distances_through_cells
                    )
                    mean_squared_error += abs2(sim_travel_time - experimental_travel_time)
                end
            end
        end
    end

    return sol.t, u_named, mean_squared_error
end

function regenerate_fvm_state(sol, system)
    u_list = ComponentVector[]

    u_named = [ComponentVector(sol.u[i], system.state_axes) for i in eachindex(sol.u)]

    du_temporary = similar(merge_properties(u_named[1], deepcopy(ComponentVector(system.properties_vec, system.properties_axes))))

    for i in eachindex(sol.t)
        u = merge_properties(u_named[i], deepcopy(ComponentVector(system.properties_vec, system.properties_axes)))
        u = merge_properties(u, deepcopy(ComponentVector(system.cache_vec, system.cache_axes)))
        
        for cell_id in grid.cellsets["fluid"]
            update_fluid_properties!(du_temporary, u, cell_id, geo.cell_volumes[cell_id], system)
        end

        push!(u_list, u)
    end

    return u_list
end

u_complete_list = regenerate_fvm_state(sol, system);

root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

sol_to_vtk(sol, u_complete_list, grid, @__FILE__, root_dir)

loss(1.0)