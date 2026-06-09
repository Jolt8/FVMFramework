
function build_simulated_travel_time_matrix(tMax, saveat)
    n_timesteps = Int(tMax / saveat)
    simulated_travel_time_matrix = [zeros(Float64, n_timesteps) for i in 1:length(transducer_ids), j in 1:length(transducer_ids)]

    return simulated_travel_time_matrix
end

function calculate_simulated_travel_time(u, origin_transducer_id, destination_transducer_id, ray_map_intersected_cells, ray_map_distances_through_cells)
    intersected_cells = ray_map_intersected_cells[origin_transducer_id, destination_transducer_id]
    intersected_cell_distances = ray_map_distances_through_cells[origin_transducer_id, destination_transducer_id]

    travel_time = 0.0

    for (index, cell_id) in enumerate(intersected_cells)
        travel_time += intersected_cell_distances[index] / u.speed_of_sound[cell_id]
    end

    return travel_time
end

function calculate_simulated_received_returned_volume(u, t, origin_transducer_id, destination_transducer_id, ray_map_intersected_cells, ray_map_distances_through_cells, initial_ping_volume)
    #the way this function is going to work is that in our loss function, we're going to be comparing the simulated return volume from one transducer to another 
    #to the experimentally measured return volume. To satisfy the loss function, variables affecting the void fraction and bubble_radii will be adjusted to make the simulated and experimental 
    #return volumes match.
    
    intersected_cells = ray_map_intersected_cells[origin_transducer_id, destination_transducer_id]
    intersected_cell_distances = ray_map_distances_through_cells[origin_transducer_id, destination_transducer_id]

    overall_attenuation_coefficient = 0.0

    #just using a simple attenuation_coefficient model for now 
    for (index, cell_id) in enumerate(intersected_cells)
        overall_attenuation_coefficient += (u.gas_holdup[cell_id] * u.bubble_radii[cell_id]) * intersected_cell_distances[index]
        #TODO: this may have to be changed to a different model in the future, this model is quite simple
    end

    overall_attenuation_coefficient /= ray_map_ray_lengths[origin_transducer_id, destination_transducer_id]

    returned_volume = initial_ping_volume * exp(-overall_attenuation_coefficient * ray_map_ray_lengths[origin_transducer_id, destination_transducer_id])

    return returned_volume
end

function calculate_bubble_velocity!(du, u, cell_id, vol)
    #= 
    the purpose of this equation is to constrain both the drag_coefficient and the bubble_surface_area
    shouldn't bubble_volume also be constrained by the gas_holdup and cell_volume? We'll deal with that later
    #TODO: find a way to constrain gas_holdup to bubble_radii 
    although maybe this could be used to constrain the shape of the bubble since bubbles are usually not perfectly circular 
    this would make both the bubble volume and bubble surface area fall apart

    also, I'm not sure how different geometries of bubbles observed experimentally affect the actual bubble_radii
    if you had a bubble that was 2 units high, 1 unit wide, and 1 unit deep, what bubble_radii would the transducers pick up
    would they get the bubble volume right with the regular circle volume calculation, would it get the bubble surface_area_right
    how much would they diverge if the bubble shape was more and more irregular
    questions, questions, questions...
        - notes:

        - this really depends on the bubble's orientation relative to the transducer
        if we only saw the bubble with one transducer as a 1 unit x 1 unit object, the bubble could be any length because the transducer 
        cannot see into the bubble
        this means that the only real trustworthy bubble volumes come from when the bubble is observed by 2-3 transducers
        we might have to do a similar method as the bubble velocity where we use linear algebra to account for the different 
        angles a cell is observed at to figure out the bubble's actual dimensions

        - this may not matter because these properties would likely be absorbed into drag_coefficient as a lumped parameter
        however, it should be noted that this is not ideal for bubble-slurry reactors where the bubble surface area 
        is a useful parameter for determining reaction rate based on the surface area of the bubbles
        in addition, in the drift flux model, we're probably going to have to develop an empirical correlation that 
        relates evaporation/condensation rate to bubble surface area
        for example:
            kinetic_constant = u.phase_change_mass_transfer_coefficient[cell_id] * u.bubble_area_per_m3[cell_id]
        however however, a model for determining the phase_change_mass_transfer_coefficient could also lump in the bubble surface area
        although, this doesn't feel great and I'm sure there's a way to roughly determine bubble surface area

        - that being said, the buoyant force is hugely effected by the estimated bubble radius 
        which means that the drag_coefficient is going to become even more distorted to account for the wrongly measured bubble radius

        - I've confimed that the way to handle this issue is to reconstruct the 3D structure of each bubble using linear algebra like 
        in the finalize_velocities! in precompute_experimental_velocities_over_time. 
        - This would make each bubble an elipsoid

    =#

    
    bubble_volume = (4 / 3) * pi * u.bubble_radii[cell_id]^3
    buoyant_force = (u.fluid_density[cell_id] - u.gas_density[cell_id]) * bubble_volume * 9.81

    bubble_surface_area = 4 * pi * u.bubble_radii[cell_id]^2
    drag_force = (1/2) * u.drag_coefficient[cell_id] * u.fluid_density[cell_id] * (u.gas_velocity[cell_id] - u.liquid_velocity[cell_id])^2 * bubble_surface_area

    bubble_mass = u.gas_density[cell_id] * bubble_volume
    
    u.bubble_acceleration[cell_id] = (buoyant_force - drag_force) / bubble_mass
    #Hmm, I guess we're going to have to find the derivative of the simulated_bubble_velocity with respect to time to get the simulated_bubble_acceleration and then 
    #compare them in the losss function to constrain the drag_coefficient and bubble_surface area
    #this is going to be quite strange because bubble_velocity is already constrained by the z component of precompute_experimental_velocities_over_time
    #I think the best way to do this would be to get the derivative of the experimental_vel_z_interps and compare that to bubble_acceleration in the loss function

    #ok yes, I've confirmed that the correct way to use this in a loss function is to find the derivative of experimental_vel_z_interps and compare that to the 
    #bubble acceleration found here
    #we just have to do DataInterpolations.derivative(experimental_vel_z_interps[cell_id], t) and compare that value to du.bubble_velocity[cell_id]
end

function reconstruct_simulated_speeds_of_sound(
    u, geo,
    transducer_opposing_pairs_ids,
    center_ray_intersected_cells, center_ray_distances_through_cells,
    ray_map_intersected_cells, ray_map_distances_through_cells;
    regularization_lambda = 1e-4
)
    # Build sparse A matrix (ray distances through cells) and d vector (travel times)
    # Each row is a ray, each column is a cell, A[i,j] = distance ray i travels through cell j
    # The system A * s = d solves for per-cell slowness s = 1/speed_of_sound
    n_cells = length(geo.cell_volumes)
    row_indices = Int[]
    col_indices = Int[]
    a_values = Float64[]
    travel_times = Float64[]
    reconstructed_speeds_of_sound = zeros(Float64, n_cells)
    
    ray_idx = 0
    for (transducer_a_id, transducer_b_id) in transducer_opposing_pairs_ids
        # Ray from A to B
        ray_idx += 1
        ray_cells = center_ray_intersected_cells[transducer_a_id]
        ray_dists = center_ray_distances_through_cells[transducer_a_id]
        
        delta_t = calculate_simulated_travel_time(u, transducer_a_id, transducer_b_id, ray_map_intersected_cells, ray_map_distances_through_cells)
        push!(travel_times, delta_t)
        
        for (j, cell_id) in enumerate(ray_cells)
            push!(row_indices, ray_idx)
            push!(col_indices, cell_id)
            push!(a_values, ray_dists[j])
        end
        
        # Ray from B to A
        ray_idx += 1
        ray_cells = center_ray_intersected_cells[transducer_b_id]
        ray_dists = center_ray_distances_through_cells[transducer_b_id]
        
        delta_t = calculate_simulated_travel_time(u, transducer_b_id, transducer_a_id, ray_map_intersected_cells, ray_map_distances_through_cells)
        push!(travel_times, delta_t)
        
        for (j, cell_id) in enumerate(ray_cells)
            push!(row_indices, ray_idx)
            push!(col_indices, cell_id)
            push!(a_values, ray_dists[j])
        end
    end
    
    n_rays = ray_idx
    A = sparse(row_indices, col_indices, a_values, n_rays, n_cells)
    d = travel_times
    
    # Tikhonov regularized solution: min ||As - d||² + λ||s - s_prior||²
    # Normal equations: (AᵀA + λI)s = Aᵀd + λ * s_prior
    AᵀA = A' * A
    Aᵀd = A' * d
    
    s_prior_vec = 1.0 ./ u.speed_of_sound[1:n_cells]
    regularized_rhs = Aᵀd .+ regularization_lambda .* s_prior_vec
    
    # Solve the regularized system (AᵀA + λI is sparse with diagonal dominance)
    s = (AᵀA + regularization_lambda * I) \ regularized_rhs
    
    # Convert slowness to speed of sound
    for cell_id in 1:n_cells
        if s[cell_id] > 0.0
            reconstructed_speeds_of_sound[cell_id] = 1.0 / s[cell_id]
        end
    end
    
    # Neighbor-smooth any cells that still have zero speed (isolated from all rays)
    #TODO: we may just want to exclude comparisons between cells whose values are derived from neighboring cells in the loss function
    #this would be done to prevent cells derived from neighbors from having a higher weighting in the loss function
    while any(==(0.0), reconstructed_speeds_of_sound)
        new_speeds = copy(reconstructed_speeds_of_sound)
        
        updates_made = false
        
        for cell_id in eachindex(reconstructed_speeds_of_sound)
            if reconstructed_speeds_of_sound[cell_id] == 0.0
                sum_speed_of_sound = 0.0
                count = 0
                
                neighbor_cells = geo.cell_neighbors[cell_id][2]
                for (neighbor_id, face_idx) in neighbor_cells
                    if neighbor_id != 0 && reconstructed_speeds_of_sound[neighbor_id] != 0.0
                        sum_speed_of_sound += reconstructed_speeds_of_sound[neighbor_id]
                        count += 1
                    end
                end
                
                if count > 0
                    new_speeds[cell_id] = sum_speed_of_sound / count
                    updates_made = true
                end
            end
        end
        
        if !updates_made
            break
        end
        
        reconstructed_speeds_of_sound .= new_speeds
    end

    return reconstructed_speeds_of_sound
end