using DSP
using DataInterpolations
using LinearAlgebra
using Ferrite

function precompute_experimental_speeds_of_sound(transducer_opposing_pairs_ids, experimental_double_pulse_echo_data, ray_map_intersected_cells, ray_map_ray_lengths, geo)
    cell_speeds_of_sound = [zeros(Float64, length(geo.cell_volumes)) for _ in eachindex(experimental_double_pulse_echo_data.timestamps)]
    cell_n_times_visited = [zeros(Int, length(geo.cell_volumes)) for _ in eachindex(experimental_double_pulse_echo_data.timestamps)]
    
    for timestamp_idx in eachindex(experimental_double_pulse_echo_data.timestamps)
        for (transducer_a_id, transducer_b_id) in transducer_opposing_pairs_ids
            pulse_a_echo_b_data = experimental_double_pulse_echo_data.transducer_to_transducer_pulse_receive_amplitudes[timestamp_idx][transducer_a_id][transducer_b_id]

            for cell_id in ray_map_intersected_cells[transducer_a_id, transducer_b_id]
                delta_t = pulse_a_echo_b_data.pulse_echo_arrival_time - pulse_a_echo_b_data.pulse_send_time
                distance = ray_map_ray_lengths[transducer_a_id, transducer_b_id]

                speed_of_sound = distance / delta_t

                if cell_speeds_of_sound[timestamp_idx][cell_id] == 0.0
                    cell_speeds_of_sound[timestamp_idx][cell_id] = speed_of_sound
                    cell_n_times_visited[timestamp_idx][cell_id] += 1
                else
                    cell_speeds_of_sound[timestamp_idx][cell_id] = (cell_speeds_of_sound[timestamp_idx][cell_id] * cell_n_times_visited[timestamp_idx][cell_id] + speed_of_sound) / (cell_n_times_visited[timestamp_idx][cell_id] + 1)
                    cell_n_times_visited[timestamp_idx][cell_id] += 1
                end
            end
        end
    end

    for timestamp_idx in eachindex(experimental_double_pulse_echo_data.timestamps)
        while any(==(0.0), cell_speeds_of_sound[timestamp_idx])
            new_speeds = copy(cell_speeds_of_sound[timestamp_idx]) #this was suggested to me to prevent direction-biased smoothing, don't know how that would work

            updates_made = false 

            for cell_id in eachindex(cell_speeds_of_sound[timestamp_idx])
                if cell_speeds_of_sound[timestamp_idx][cell_id] == 0.0        
                    sum_speed_of_sound = 0.0
                    count = 0
                
                    neighbor_cells = geo.cell_neighbors[cell_id][2]
                    for (neighbor_id, face_idx) in neighbor_cells
                        if neighbor_id != 0 && cell_speeds_of_sound[timestamp_idx][neighbor_id] != 0.0
                            sum_speed_of_sound += cell_speeds_of_sound[timestamp_idx][neighbor_id]
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

            cell_speeds_of_sound[timestamp_idx] .= new_speeds
        end
    end

    return cell_speeds_of_sound
end
    

    

function build_time_to_distance_map(speed_of_sound::Vector{Float64}, ray_cells, ray_distances)
    cumulative_time = [0.0]
    cumulative_distance = [0.0]
    
    current_time = 0.0
    current_distance = 0.0
    
    for (i, cell_id) in enumerate(ray_cells)
        # Pulse-Echo means sound goes out AND comes back, so distance counts twice!
        travel_time = (ray_distances[i] * 2) / speed_of_sound[cell_id]
        
        current_time += travel_time
        current_distance += ray_distances[i]
        
        push!(cumulative_time, current_time)
        push!(cumulative_distance, current_distance)
    end
    
    # Creates a function: f(time) = distance
    return LinearInterpolation(cumulative_distance, cumulative_time)
end

function extract_velocities_from_data(pulse_1_data, pulse_2_data, pulse_repetition_interval, samples_per_second, window_size_indices)
    n_samples = length(pulse_1_data)
    n_windows = floor(Int, n_samples / window_size_indices)
    
    window_times = Float64[]
    window_partial_velocities = Float64[]
    
    for w in 1:n_windows
        start_idx = (w - 1) * window_size_indices + 1
        end_idx = w * window_size_indices
        
        # 1. Chop out the window from both pulses
        window_1 = pulse_1_data[start_idx:end_idx]
        window_2 = pulse_2_data[start_idx:end_idx]
        
        # 2. Cross-correlate to find the lag (index shift)
        lags = xcorr(window_1, window_2)
        shift_index = argmax(lags) - length(window_1) # Center the lag around 0
        
        # 3. Calculate Time Shift (tau)
        tau = shift_index / samples_per_second
        
        # 4. Calculate what time this window represents (to map to distance later)
        center_time = ((start_idx + end_idx) / 2.0) / samples_per_second
        
        # We can't calculate velocity yet because we need the speed of sound (c) from the FVM!
        # So we just store the tau and the time.
        push!(window_times, center_time)
        push!(window_partial_velocities, tau / (2 * pulse_repetition_interval)) # Partial velocity equation: v = c * (tau / 2*T_PRI)
    end
    
    return window_times, window_partial_velocities
end

function map_experimental_velocities_to_cells!(
    beam_velocities::Vector{Vector{Float64}},
    cell_update_counts::Vector{Int},
    speed_of_sound::Vector{Float64},
    window_times,
    window_partial_velocities,
    time_to_distance_map,
    projection_cells,
    projection_distances,
    tolerance = 2.0
)
    # 1. Loop through the time windows we got from the cross-correlation
    for (i, t_arrival) in enumerate(window_times)
        
        # 2. Use the speed of sound map to find exactly how far away this window is!
        bubble_distance = time_to_distance_map(t_arrival)
        
        for cell_idx in eachindex(projection_cells)
            cell_dist = projection_distances[cell_idx]
            cell_id = projection_cells[cell_idx]
            
            # 3. If the cell is physically located at this bubble distance...
            if abs(cell_dist - bubble_distance) <= tolerance
                # 4. Finish the velocity calculation (multiply by the local speed of sound!)
                idx = cell_update_counts[cell_id] + 1
                if idx <= length(beam_velocities[cell_id])
                    beam_velocities[cell_id][idx] = speed_of_sound[cell_id] * window_partial_velocities[i]
                    cell_update_counts[cell_id] += 1
                end
            end
            
            if cell_dist > bubble_distance + tolerance
                break
            end
        end
    end
end

function finalize_velocities!(
    experimental_vel_x::Vector{Float64},
    experimental_vel_y::Vector{Float64},
    experimental_vel_z::Vector{Float64},
    beam_directions::Vector{Vector{Vec{3, Float64}}},
    beam_velocities::Vector{Vector{Float64}},
    cell_update_counts::Vector{Int}
)
    n_cells = length(experimental_vel_x)
    for cell_id in 1:n_cells
        n_updates = cell_update_counts[cell_id]
        
        if n_updates >= 3
            # Extract ONLY the valid rows that actually received data
            A = zeros(Float64, n_updates, 3)
            for i in 1:n_updates
                A[i, 1] = beam_directions[cell_id][i][1]
                A[i, 2] = beam_directions[cell_id][i][2]
                A[i, 3] = beam_directions[cell_id][i][3]
            end
            b = beam_velocities[cell_id][1:n_updates]
            
            # Check if the matrix has a rank of 3 (meaning 3 distinct angles)
            if rank(A) >= 3
                cell_velocity_vector = A \ b
                experimental_vel_x[cell_id] = cell_velocity_vector[1]
                experimental_vel_y[cell_id] = cell_velocity_vector[2]
                experimental_vel_z[cell_id] = cell_velocity_vector[3]
                continue # Successfully solved!
            end
        end
        
        # If we reach here, the cell doesn't have enough 3D data.
        # Mark it with NaNs so our patching function knows to fix it later.
        experimental_vel_x[cell_id] = NaN
        experimental_vel_y[cell_id] = NaN
        experimental_vel_z[cell_id] = NaN
    end
end

function patch_unaccounted_experimental_velocities!(experimental_vel_x::Vector{Float64}, experimental_vel_y::Vector{Float64}, experimental_vel_z::Vector{Float64}, geo)
    # Loop until no NaNs are left
    while any(isnan, experimental_vel_x)
        # Create copies so we don't pollute the averages while updating
        new_exp_x = copy(experimental_vel_x)
        new_exp_y = copy(experimental_vel_y)
        new_exp_z = copy(experimental_vel_z)
        
        for cell_id in eachindex(experimental_vel_x)
            if isnan(experimental_vel_x[cell_id])
                neighbor_cells = geo.cell_neighbors[cell_id][2]
                
                sum_vec = Vec{3, Float64}((0.0, 0.0, 0.0))
                count = 0
                
                for (neighbor_id, face_idx) in neighbor_cells
                    if neighbor_id != 0 && !isnan(experimental_vel_x[neighbor_id])
                        sum_vec += Vec{3, Float64}((experimental_vel_x[neighbor_id], experimental_vel_y[neighbor_id], experimental_vel_z[neighbor_id]))
                        count += 1
                    end
                end
                
                if count > 0
                    new_exp_x[cell_id] = sum_vec[1] / count
                    new_exp_y[cell_id] = sum_vec[2] / count
                    new_exp_z[cell_id] = sum_vec[3] / count
                end
            end
        end
        
        # If no updates were made in a pass, break to avoid an infinite loop (e.g. disjoint component)
        if new_exp_x == experimental_vel_x && new_exp_y == experimental_vel_y && new_exp_z == experimental_vel_z
            break
        end
        experimental_vel_x .= new_exp_x
        experimental_vel_y .= new_exp_y
        experimental_vel_z .= new_exp_z
    end
end

struct DoublePulseEchoData
    double_pulse_timestamps::Vector{Float64}
    # transducer_self_pulse_echo_pairs[transducer_A_id][transducer_B_id] = (pulse_1_amplitudes, pulse_2_amplitudes)
    transducer_self_pulse_echo_pairs::Vector{Vector{Tuple{Vector{Float64}, Vector{Float64}}}} 
    #transducer_to_transducer_pulse_receive[timestamp_idx][transducer_A_id][transducer_B_id] = transducer_B received signal when transducer_A pulsed
    transducer_to_transducer_pulse_receive_amplitudes::Vector{Vector{Vector{Float64}}} 
end

function precompute_experimental_velocities_over_time(
    grid, geo,
    #cell_speeds_of_sound::Vector{Float64}, we'll add this in later once we start using two-loop optimization
    transducer_ids::Vector{Int},
    transducer_opposing_pairs_ids::Vector{Tuple{Int, Int}},
    ray_map_intersected_cells,
    ray_map_distances_through_cells,
    ray_map_ray_lengths,
    transducer_projection_intersected_cells,
    transducer_projection_cell_distances,
    projection_unit_vectors_to_cell_centers,
    cell_projection_counts::Vector{Int},
    experimental_double_pulse_echo_data::DoublePulseEchoData,
    pulse_repetition_interval::Float64,
    samples_per_second::Float64,
    window_size_indices::Int;
    tolerance = 2.0
)
    n_cells = length(grid.cells)
    cell_speeds_of_sound = precompute_experimental_speeds_of_sound(
        transducer_opposing_pairs_ids, 
        experimental_double_pulse_echo_data, 
        ray_map_intersected_cells, 
        ray_map_ray_lengths, 
        geo
    )
    
    # 3. Populate beam_directions for each cell from the transducer projections
    beam_directions = [Vector{Vec{3, Float64}}(undef, cell_projection_counts[cell_id]) for cell_id in 1:n_cells]
    projection_index = zeros(Int, n_cells)
    for transducer_id in transducer_ids
        intersected_cells = transducer_projection_intersected_cells[transducer_id]
        unit_vectors = projection_unit_vectors_to_cell_centers[transducer_id]
        for (transducer_cell_idx, cell_id) in enumerate(intersected_cells)
            projection_index[cell_id] += 1
            p = projection_index[cell_id]
            beam_directions[cell_id][p] = unit_vectors[transducer_cell_idx]
        end
    end

    n_timesteps = length(experimental_double_pulse_echo_data.double_pulse_timestamps)

    experimental_discrete_vel_x = [zeros(Float64, n_cells) for _ in 1:n_timesteps]
    experimental_discrete_vel_y = [zeros(Float64, n_cells) for _ in 1:n_timesteps]
    experimental_discrete_vel_z = [zeros(Float64, n_cells) for _ in 1:n_timesteps]
    
    # 4. Loop over each timestep
    for timestamp_idx in 1:n_timesteps
        cell_update_counts = zeros(Int, n_cells)
        beam_velocities = [zeros(Float64, cell_projection_counts[cell_id]) for cell_id in 1:n_cells]
        
        # Process each transducer beam
        for (transducer_A_id, transducer_B_id) in transducer_opposing_pairs_ids
            # Get center ray for this transducer's beam
            ray_cells = ray_map_intersected_cells[transducer_A_id, transducer_B_id]
            ray_distances = ray_map_distances_through_cells[transducer_A_id, transducer_B_id]
            
            # Build mapping
            time_to_distance_map = build_time_to_distance_map(cell_speeds_of_sound, ray_cells, ray_distances)
            
            # Extract velocities
            pulse_1 = experimental_double_pulse_echo_data.transducer_self_pulse_echo_pairs[transducer_A_id][timestamp_idx][1]
            pulse_2 = experimental_double_pulse_echo_data.transducer_self_pulse_echo_pairs[transducer_A_id][timestamp_idx][2]
            
            window_times, window_partial_velocities = extract_velocities_from_data(
                pulse_1, pulse_2, pulse_repetition_interval, samples_per_second, window_size_indices
            )
            
            # Map experimental velocities to cells
            proj_cells = transducer_projection_intersected_cells[transducer_A_id]
            proj_dists = transducer_projection_cell_distances[transducer_A_id]
            
            map_experimental_velocities_to_cells!(
                beam_velocities, cell_update_counts, cell_speeds_of_sound,
                window_times, window_partial_velocities, time_to_distance_map,
                proj_cells, proj_dists, tolerance
            )
        end
        
        # Solve least squares system for this timestep
        finalize_velocities!(
            experimental_discrete_vel_x[timestamp_idx], experimental_discrete_vel_y[timestamp_idx], experimental_discrete_vel_z[timestamp_idx], 
            beam_directions, beam_velocities, cell_update_counts
        )
        
        # Patch unresolved cells
        patch_unaccounted_experimental_velocities!(
            experimental_discrete_vel_x[timestamp_idx], experimental_discrete_vel_y[timestamp_idx], experimental_discrete_vel_z[timestamp_idx], 
            geo
        )
    end

    experimental_vel_x_interp = Vector{LinearInterpolation}(undef, n_cells)
    experimental_vel_y_interp = Vector{LinearInterpolation}(undef, n_cells)
    experimental_vel_z_interp = Vector{LinearInterpolation}(undef, n_cells)

    for cell_id in 1:n_cells
        cell_vel_x_history = [experimental_discrete_vel_x[t][cell_id] for t in 1:n_timesteps]
        cell_vel_y_history = [experimental_discrete_vel_y[t][cell_id] for t in 1:n_timesteps]
        cell_vel_z_history = [experimental_discrete_vel_z[t][cell_id] for t in 1:n_timesteps]

        experimental_vel_x_interp[cell_id] = LinearInterpolation(
            experimental_double_pulse_echo_data.double_pulse_timestamps, 
            cell_vel_x_history, 
            extrapolate = true
        )

        experimental_vel_y_interp[cell_id] = LinearInterpolation(
            experimental_double_pulse_echo_data.double_pulse_timestamps, 
            cell_vel_y_history, 
            extrapolate = true
        )

        experimental_vel_z_interp[cell_id] = LinearInterpolation(
            experimental_double_pulse_echo_data.double_pulse_timestamps, 
            cell_vel_z_history, 
            extrapolate = true
        )
    end
    
    return experimental_vel_x_interp, experimental_vel_y_interp, experimental_vel_z_interp
end