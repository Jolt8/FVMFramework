using DSP
using DataInterpolations
using LinearAlgebra
using SparseArrays
using FFTW
using Ferrite

function precompute_experimental_speeds_of_sound(
    transducer_opposing_pairs_ids, experimental_sonic_data,
    center_ray_intersected_cells, center_ray_distances_through_cells,
    n_cells, geo;
    regularization_lambda = 1e-4, #this controls how much you trust the data vs the prior data. 
    #higher = more trust in prior data
    #lower = less trust in prior data
    prior_speed_of_sound = 1500.0,
    timestamp_tolerance = Inf #set to Inf to always include the closest measurement
)
    # 1. Collect ALL unique timestamps across all pairs
    #    Each transducer pair fires independently, so we build a global time grid
    #    from the union of all send timestamps across all pairs.
    all_timestamps = Float64[]

    for (transducer_a_id, transducer_b_id) in transducer_opposing_pairs_ids
        pair_data_ab = experimental_sonic_data.transducer_to_transducer_send_receive_timestamps[transducer_a_id, transducer_b_id]
        for ts in pair_data_ab
            push!(all_timestamps, ts.send)
        end
        pair_data_ba = experimental_sonic_data.transducer_to_transducer_send_receive_timestamps[transducer_b_id, transducer_a_id]
        for ts in pair_data_ba
            push!(all_timestamps, ts.send)
        end
    end

    unique_timestamps = sort(unique(all_timestamps))
    n_timestamps = length(unique_timestamps)
    cell_speeds_of_sound = [zeros(Float64, n_cells) for _ in 1:n_timestamps]
    cell_speeds_of_sound_uncertainties = [zeros(Float64, n_cells) for _ in 1:n_timestamps]
    s_prior = 1.0 / prior_speed_of_sound

    # 2. For each global timestamp, build a linear system from all pairs that fired near it
    for (t_idx, t_global) in enumerate(unique_timestamps)
        row_indices = Int[]
        col_indices = Int[]
        a_values = Float64[]
        travel_times = Float64[]

        ray_idx = 0
        for (transducer_a_id, transducer_b_id) in transducer_opposing_pairs_ids
            # --- Ray A→B ---
            pair_data_ab = experimental_sonic_data.transducer_to_transducer_send_receive_timestamps[transducer_a_id, transducer_b_id]
            send_times_ab = [ts.send for ts in pair_data_ab]
            best_idx_ab = argmin(abs.(send_times_ab .- t_global))
            best_time_ab = send_times_ab[best_idx_ab]

            if abs(best_time_ab - t_global) <= timestamp_tolerance
                delta_t_ab = pair_data_ab[best_idx_ab].receive - pair_data_ab[best_idx_ab].send

                ray_idx += 1
                ray_cells = center_ray_intersected_cells[transducer_a_id]
                ray_dists = center_ray_distances_through_cells[transducer_a_id]
                push!(travel_times, delta_t_ab)
                for (j, cell_id) in enumerate(ray_cells)
                    push!(row_indices, ray_idx)
                    push!(col_indices, cell_id)
                    push!(a_values, ray_dists[j])
                end

                #uncertainty: trust measurement less when more cells are along the path
                cell_speeds_of_sound_uncertainties[t_idx][ray_cells] .+= length(ray_cells)
            end

            # --- Ray B→A ---
            pair_data_ba = experimental_sonic_data.transducer_to_transducer_send_receive_timestamps[transducer_b_id, transducer_a_id]
            send_times_ba = [ts.send for ts in pair_data_ba]
            best_idx_ba = argmin(abs.(send_times_ba .- t_global))
            best_time_ba = send_times_ba[best_idx_ba]

            if abs(best_time_ba - t_global) <= timestamp_tolerance
                delta_t_ba = pair_data_ba[best_idx_ba].receive - pair_data_ba[best_idx_ba].send

                ray_idx += 1
                ray_cells = center_ray_intersected_cells[transducer_b_id]
                ray_dists = center_ray_distances_through_cells[transducer_b_id]
                push!(travel_times, delta_t_ba)
                for (j, cell_id) in enumerate(ray_cells)
                    push!(row_indices, ray_idx)
                    push!(col_indices, cell_id)
                    push!(a_values, ray_dists[j])
                end

                # uncertainty: trust measurement less when more cells are along the path
                cell_speeds_of_sound_uncertainties[t_idx][ray_cells] .+= length(ray_cells)
            end
        end

        # Solve Tikhonov system for this timestamp
        n_rays = ray_idx
        if n_rays == 0
            # No rays contributed to this timestamp — fill with prior
            cell_speeds_of_sound[t_idx] .= prior_speed_of_sound
            continue
        end

        A = sparse(row_indices, col_indices, a_values, n_rays, n_cells)
        d = travel_times

        # Tikhonov regularized solution: min ||As - d||² + λ||s - s_prior||²
        # Normal equations: (AᵀA + λI)s = Aᵀd + λ * s_prior
        AᵀA = A' * A
        Aᵀd = A' * d

        s_prior_vec = fill(s_prior, n_cells)
        regularized_rhs = Aᵀd .+ regularization_lambda .* s_prior_vec

        # Solve the regularized system
        s = (AᵀA + regularization_lambda * I) \ regularized_rhs

        # Convert slowness to speed of sound
        for cell_id in 1:n_cells
            if s[cell_id] > 0.0
                cell_speeds_of_sound[t_idx][cell_id] = 1.0 / s[cell_id]
            end
        end

        # Neighbor-smooth any cells that still have zero speed (isolated from all rays)
        #TODO: we may just want to exclude comparisons between cells whose values are derived 
        #from neighboring cells in the loss function, to prevent neighbor-derived cells from 
        #having an inflated weighting
        while any(==(0.0), cell_speeds_of_sound[t_idx])
            new_speeds = copy(cell_speeds_of_sound[t_idx])

            updates_made = false

            for cell_id in eachindex(cell_speeds_of_sound[t_idx])
                if cell_speeds_of_sound[t_idx][cell_id] == 0.0
                    sum_speed_of_sound = 0.0
                    count = 0

                    neighbor_cells = geo.cell_neighbors[cell_id][2]
                    for (neighbor_id, face_idx) in neighbor_cells
                        if neighbor_id != 0 && cell_speeds_of_sound[t_idx][neighbor_id] != 0.0
                            sum_speed_of_sound += cell_speeds_of_sound[t_idx][neighbor_id]
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

            cell_speeds_of_sound[t_idx] .= new_speeds
        end
    end

    # 3. Build per-cell interpolations over the global timestamps
    #    Each cell gets an interpolation: f(t) -> speed_of_sound at that cell
    speeds_of_sound_interps = Vector{LinearInterpolation}(undef, n_cells)
    speeds_of_sound_uncertainty_interps = Vector{LinearInterpolation}(undef, n_cells)
    #actually, maybe creating an uncertainty interp is incorrect because the uncertainty of a single cell can vary significanlty depending on which 
    #pairs of transducers were used to determine its velocity
    #this is because I'm guessing that the uncertainty would vary significanlty almost every second
    
    for cell_id in 1:n_cells
        cell_speed_history = [cell_speeds_of_sound[t_idx][cell_id] for t_idx in 1:n_timestamps]
        cell_speeds_of_sound_uncertainty_history = [cell_speeds_of_sound_uncertainties[t_idx][cell_id] for t_idx in 1:n_timestamps]
        speeds_of_sound_interps[cell_id] = LinearInterpolation(cell_speed_history, unique_timestamps)
        speeds_of_sound_uncertainty_interps[cell_id] = LinearInterpolation(cell_speeds_of_sound_uncertainty_history, unique_timestamps)
    end

    #TODO: remove this later, this is just to see how much the uncertainty varies with time and if it needs further processing
    @display plot(unique_timestamps, [speeds_of_sound_interps[1](t) for t in unique_timestamps])

    return cell_speeds_of_sound, speeds_of_sound_interps, speeds_of_sound_uncertainty_interps, unique_timestamps
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
    return LinearInterpolation(cumulative_distance, cumulative_time) #DataInterpolations is (Dependent, Independent)
end

function extract_velocities_from_data(pulse_1_data, pulse_2_data, pulse_repetition_interval, samples_per_second, window_size_indices)
    n_samples = length(pulse_1_data)
    n_windows = floor(Int, n_samples / window_size_indices)
    
    window_times = Float64[]
    window_partial_velocities = Float64[]
    window_backscatter_amplitudes = Float64[]
    window_turbulence_widths = Float64[]
    
    for w in 1:n_windows
        start_idx = (w - 1) * window_size_indices + 1
        end_idx = w * window_size_indices
        
        # 1. Chop out the window from both pulses
        window_1 = pulse_1_data[start_idx:end_idx]
        window_2 = pulse_2_data[start_idx:end_idx]
        
        # 2. Cross-correlate to find the lag (index shift)
        lags = xcorr(window_1, window_2)
        peak_idx = argmax(lags)
        
        # 3. Sub-sample refinement via parabolic interpolation around the peak
        if peak_idx > 1 && peak_idx < length(lags)
            α = lags[peak_idx - 1]
            β = lags[peak_idx]
            γ = lags[peak_idx + 1]
            denominator = α - 2β + γ
            if abs(denominator) > 1e-12
                refined_offset = 0.5 * (α - γ) / denominator
            else
                refined_offset = 0.0
            end
            shift_index = (peak_idx - length(window_1)) + refined_offset
        else
            shift_index = peak_idx - length(window_1)
        end
        
        # 4. Calculate Time Shift (tau)
        tau = shift_index / samples_per_second
        
        # 5. Calculate what time this window represents (to map to distance later)
        center_time = ((start_idx + end_idx) / 2.0) / samples_per_second
        
        # We can't calculate velocity yet because we need the speed of sound (c) from the FVM!
        # So we just store the tau and the time.
        push!(window_times, center_time)
        push!(window_partial_velocities, tau / (2 * pulse_repetition_interval)) # Partial velocity equation: v = c * (tau / 2*T_PRI)

        # ── Feature 1: Backscatter amplitude (RMS envelope of the window) ────────
        # The RMS amplitude at each depth window is proportional to local scatterer
        # reflectivity (gas_holdup × bubble_cross_section). This gives an independent
        # constraint on void fraction from the self-to-self echo path.
        rms_amplitude = sqrt(sum(window_1 .^ 2) / length(window_1))
        push!(window_backscatter_amplitudes, rms_amplitude)

        # ── Feature 2: Spectral broadening / turbulence intensity ────────────────
        # The width of the cross-correlation peak contains velocity variance info.
        # A sharp peak = uniform flow; a broad peak = turbulent mixing.
        # Fit log(|corr|) to a parabola over ±3 points around the peak to get σ.
        half_fit_width = 3
        if peak_idx > half_fit_width && peak_idx <= length(lags) - half_fit_width
            fit_indices = (peak_idx - half_fit_width):(peak_idx + half_fit_width)
            fit_values = lags[fit_indices]

            peak_val = lags[peak_idx]
            if peak_val > 1e-12
                # Normalize by peak and take log (clamp to avoid log(0))
                normalized = max.(fit_values ./ peak_val, 1e-10)
                log_corr = log.(normalized)

                # Fit parabola: log(G(x)) ≈ a × x² + b, where a = -1/(2σ²)
                x = collect(-half_fit_width:half_fit_width)
                x2 = x .^ 2
                n_fit = length(x)
                mean_log = sum(log_corr) / n_fit
                mean_x2 = sum(x2) / n_fit
                a_coeff = sum(x2 .* (log_corr .- mean_log)) / sum(x2 .* (x2 .- mean_x2))

                if a_coeff < -1e-12
                    σ_indices = sqrt(-1.0 / (2.0 * a_coeff))
                    # Store as "partial" turbulence width (multiply by c later)
                    σ_partial = σ_indices / (2.0 * samples_per_second * pulse_repetition_interval)
                    push!(window_turbulence_widths, σ_partial)
                else
                    push!(window_turbulence_widths, 0.0)
                end
            else
                push!(window_turbulence_widths, 0.0)
            end
        else
            push!(window_turbulence_widths, 0.0)
        end
    end
    
    return window_times, window_partial_velocities, window_backscatter_amplitudes, window_turbulence_widths
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
    cell_velocity_uncertainties = [[0.0, 0.0, 0.0] for _ in 1:n_cells]
    for cell_id in 1:n_cells
        n_updates = cell_update_counts[cell_id]
        
        if n_updates >= 1
            # Build the direction matrix with only valid entries
            A = zeros(Float64, n_updates, 3)
            for i in 1:n_updates
                A[i, 1] = beam_directions[cell_id][i][1]
                A[i, 2] = beam_directions[cell_id][i][2]
                A[i, 3] = beam_directions[cell_id][i][3]
            end
            b = beam_velocities[cell_id][1:n_updates]

            # Use condition number of the direction matrix as uncertainty metric
            # High condition number = beams are nearly parallel = poorly constrained
            cell_velocity_uncertainties[cell_id][1:3] .= cond(A[1:n_updates, :])
            
            # A \ b handles all cases automatically:
            #   1 beam:  minimum-norm solution (velocity projected onto beam direction)
            #   2 beams: minimum-norm solution (velocity in the plane of both beams)
            #   3+ beams: least-squares solution (full 3D reconstruction)
            cell_velocity_vector = A \ b
            experimental_vel_x[cell_id] = cell_velocity_vector[1]
            experimental_vel_y[cell_id] = cell_velocity_vector[2]
            experimental_vel_z[cell_id] = cell_velocity_vector[3]
        else
            # No beam data at all — mark for neighbor patching
            experimental_vel_x[cell_id] = NaN
            experimental_vel_y[cell_id] = NaN
            experimental_vel_z[cell_id] = NaN
        end
    end

    return cell_velocity_uncertainties
end

function patch_unaccounted_experimental_velocities!(experimental_vel_x::Vector{Float64}, experimental_vel_y::Vector{Float64}, experimental_vel_z::Vector{Float64}, geo)
    # Loop until no NaNs are left
    while any(isnan, experimental_vel_x)
        # Create copies so we don't pollute the averages while updating
        new_exp_x = copy(experimental_vel_x)
        new_exp_y = copy(experimental_vel_y)
        new_exp_z = copy(experimental_vel_z)

        updates_made = false
        
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
                    updates_made = true
                end
            end
        end
        
        # If no updates were made in a pass, break to avoid an infinite loop (e.g. disjoint component)
        if !updates_made
            break
        end

        experimental_vel_x .= new_exp_x
        experimental_vel_y .= new_exp_y
        experimental_vel_z .= new_exp_z
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 1: Backscatter amplitude — mapping and finalization
# ═══════════════════════════════════════════════════════════════════════════════

function map_experimental_backscatter_to_cells!(
    cell_backscatter_sum::Vector{Float64},
    cell_backscatter_counts::Vector{Int},
    window_times,
    window_backscatter_amplitudes,
    time_to_distance_map,
    projection_cells,
    projection_distances,
    tolerance = 2.0
)
    for (i, t_arrival) in enumerate(window_times)
        bubble_distance = time_to_distance_map(t_arrival)

        for cell_idx in eachindex(projection_cells)
            cell_dist = projection_distances[cell_idx]
            cell_id = projection_cells[cell_idx]

            if abs(cell_dist - bubble_distance) <= tolerance
                cell_backscatter_sum[cell_id] += window_backscatter_amplitudes[i]
                cell_backscatter_counts[cell_id] += 1
            end

            if cell_dist > bubble_distance + tolerance
                break
            end
        end
    end
end

function finalize_backscatter!(
    cell_backscatter_intensity::Vector{Float64},
    cell_backscatter_sum::Vector{Float64},
    cell_backscatter_counts::Vector{Int}
)
    for cell_id in eachindex(cell_backscatter_intensity)
        if cell_backscatter_counts[cell_id] > 0
            # Average backscatter across all beam observations
            cell_backscatter_intensity[cell_id] = cell_backscatter_sum[cell_id] / cell_backscatter_counts[cell_id]
        else
            cell_backscatter_intensity[cell_id] = NaN
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 2: Spectral broadening / turbulence intensity — mapping and finalization
# ═══════════════════════════════════════════════════════════════════════════════

function map_experimental_turbulence_to_cells!(
    cell_turbulence_sum_sq::Vector{Float64},
    cell_turbulence_counts::Vector{Int},
    speed_of_sound::Vector{Float64},
    window_times,
    window_turbulence_widths,
    time_to_distance_map,
    projection_cells,
    projection_distances,
    tolerance = 2.0
)
    for (i, t_arrival) in enumerate(window_times)
        bubble_distance = time_to_distance_map(t_arrival)

        for cell_idx in eachindex(projection_cells)
            cell_dist = projection_distances[cell_idx]
            cell_id = projection_cells[cell_idx]

            if abs(cell_dist - bubble_distance) <= tolerance
                # Convert partial turbulence width to velocity units: σ_v = c × σ_partial
                # Sum of squares for RMS averaging across beams
                σ_v = speed_of_sound[cell_id] * window_turbulence_widths[i]
                cell_turbulence_sum_sq[cell_id] += σ_v^2
                cell_turbulence_counts[cell_id] += 1
            end

            if cell_dist > bubble_distance + tolerance
                break
            end
        end
    end
end

function finalize_turbulence!(
    cell_turbulence_intensity::Vector{Float64},
    cell_turbulence_sum_sq::Vector{Float64},
    cell_turbulence_counts::Vector{Int}
)
    for cell_id in eachindex(cell_turbulence_intensity)
        if cell_turbulence_counts[cell_id] > 0
            # RMS of velocity standard deviations across beams
            cell_turbulence_intensity[cell_id] = sqrt(cell_turbulence_sum_sq[cell_id] / cell_turbulence_counts[cell_id])
        else
            cell_turbulence_intensity[cell_id] = NaN
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Generic scalar field NaN-patching (for backscatter and turbulence)
# ═══════════════════════════════════════════════════════════════════════════════

function patch_unaccounted_scalar_field!(field::Vector{Float64}, geo)
    while any(isnan, field)
        new_field = copy(field)
        updates_made = false

        for cell_id in eachindex(field)
            if isnan(field[cell_id])
                neighbor_cells = geo.cell_neighbors[cell_id][2]
                total = 0.0
                count = 0

                for (neighbor_id, face_idx) in neighbor_cells
                    if neighbor_id != 0 && !isnan(field[neighbor_id])
                        total += field[neighbor_id]
                        count += 1
                    end
                end

                if count > 0
                    new_field[cell_id] = total / count
                    updates_made = true
                end
            end
        end

        if !updates_made
            break
        end

        field .= new_field
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 3: Spectral shape extraction from received signals
# Normalizes by the center frequency to cancel the transmitted spectrum shape.
# ═══════════════════════════════════════════════════════════════════════════════

function extract_normalized_spectral_shape(
    received_amplitude_signal::Vector{Float64},
    samples_per_second::Float64,
    transducer_center_frequency::Float64;
    n_frequency_bins = 10,
    bandwidth_factor = 0.25
)
    n_samples = length(received_amplitude_signal)

    # FFT of received signal (magnitude only)
    spectrum = abs.(fft(received_amplitude_signal))
    freqs = collect(0:n_samples - 1) .* (samples_per_second / n_samples)

    # Positive frequencies only
    n_positive = div(n_samples, 2)
    spectrum = spectrum[1:n_positive]
    freqs = freqs[1:n_positive]

    # Usable frequency range: within ±2Γ of center (Lorentzian half-width Γ = f₀ × bandwidth_factor)
    Γ = transducer_center_frequency * bandwidth_factor
    f_min = transducer_center_frequency - 2.0 * Γ
    f_max = transducer_center_frequency + 2.0 * Γ

    usable_mask = (freqs .>= f_min) .& (freqs .<= f_max)
    usable_freqs = freqs[usable_mask]
    usable_spectrum = spectrum[usable_mask]

    if length(usable_freqs) < 2
        return Float64[], Float64[]
    end

    # Bin the usable range into n_frequency_bins
    bin_edges = range(f_min, f_max, length = n_frequency_bins + 1)

    bin_center_freqs = Float64[]
    bin_log_amplitudes = Float64[]

    for bin in 1:n_frequency_bins
        f_lo = bin_edges[bin]
        f_hi = bin_edges[bin + 1]
        bin_mask = (usable_freqs .>= f_lo) .& (usable_freqs .< f_hi)

        n_in_bin = sum(bin_mask)
        if n_in_bin > 0
            avg_amplitude = sum(usable_spectrum[bin_mask]) / n_in_bin
            if avg_amplitude > 1e-30
                push!(bin_center_freqs, (f_lo + f_hi) / 2.0)
                push!(bin_log_amplitudes, log(avg_amplitude))
            end
        end
    end

    if length(bin_log_amplitudes) < 2
        return Float64[], Float64[]
    end

    # Normalize by center frequency bin — this cancels the transmitted spectrum shape
    center_bin_idx = argmin(abs.(bin_center_freqs .- transducer_center_frequency))
    center_log_amplitude = bin_log_amplitudes[center_bin_idx]
    normalized_shape = bin_log_amplitudes .- center_log_amplitude

    return bin_center_freqs, normalized_shape
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 4: Bubble passage frequency from backscatter time series
# Uses zero-crossing rate of detrended backscatter in a sliding window.
# ═══════════════════════════════════════════════════════════════════════════════

function compute_bubble_passage_frequencies(
    cell_backscatter_history,  # Vector of Vector{Float64}: [timestamp_idx] → per-cell backscatter
    timestamps::Vector{Float64},
    n_cells::Int;
    passage_half_window = 5
)
    n_timestamps = length(timestamps)
    cell_passage_frequencies = [zeros(Float64, n_cells) for _ in 1:n_timestamps]

    for cell_id in 1:n_cells
        # Extract backscatter time series for this cell
        bs_series = [cell_backscatter_history[t][cell_id] for t in 1:n_timestamps]

        # Detrend by subtracting the mean
        mean_val = sum(bs_series) / n_timestamps
        detrended = bs_series .- mean_val

        for t_idx in 1:n_timestamps
            w_start = max(1, t_idx - passage_half_window)
            w_end = min(n_timestamps, t_idx + passage_half_window)

            if w_end - w_start < 2
                cell_passage_frequencies[t_idx][cell_id] = 0.0
                continue
            end

            # Count zero-crossings in window
            n_crossings = 0
            for i in (w_start + 1):w_end
                if detrended[i] * detrended[i - 1] < 0.0
                    n_crossings += 1
                end
            end

            window_duration = timestamps[w_end] - timestamps[w_start]
            if window_duration > 0.0
                # Each full cycle has 2 zero crossings
                cell_passage_frequencies[t_idx][cell_id] = n_crossings / (2.0 * window_duration)
            end
        end
    end

    return cell_passage_frequencies
end

struct SonicData
    #= self_to_self_echo_pairs[transducer_id] is a Vector of timestamp entries,
       each entry is a Tuple of two RF signal vectors (pulse_1, pulse_2)
       Access: self_to_self_echo_pairs[transducer_id][timestamp_idx] -> (pulse_1_vec, pulse_2_vec)
    =#
    self_to_self_echo_pairs::Vector{Vector{Tuple{Vector{Float64}, Vector{Float64}}}}

    #= self_to_self_send_echo_timestamps[transducer_id] is a Vector of NamedTuples over timestamps
       Access: self_to_self_send_echo_timestamps[transducer_id][timestamp_idx].send
               self_to_self_send_echo_timestamps[transducer_id][timestamp_idx].echo
    =#
    self_to_self_send_echo_timestamps::Vector{Vector{NamedTuple{(:send, :echo), Tuple{Float64, Float64}}}}

    #= self_to_self_echo_frequencies[transducer_id] is a Vector of frequencies over timestamps
       Access: self_to_self_echo_frequencies[transducer_id][timestamp_idx]
    =#
    self_to_self_echo_frequencies::Vector{Vector{Float64}}

    #= transducer_to_transducer_receive_amplitudes[transducer_a_id, transducer_b_id] is a Vector of 
       amplitude vectors over timestamps
       Access: transducer_to_transducer_receive_amplitudes[transducer_a_id, transducer_b_id][timestamp_idx]
    =#
    transducer_to_transducer_receive_amplitudes::Matrix{Vector{Vector{Float64}}}

    #= transducer_to_transducer_send_receive_timestamps[transducer_a_id, transducer_b_id] is a Vector of NamedTuples
       Access: transducer_to_transducer_send_receive_timestamps[transducer_a_id, transducer_b_id][timestamp_idx].send
               transducer_to_transducer_send_receive_timestamps[transducer_a_id, transducer_b_id][timestamp_idx].receive
    =#
    transducer_to_transducer_send_receive_timestamps::Matrix{Vector{NamedTuple{(:send, :receive), Tuple{Float64, Float64}}}}

    pulse_repetition_interval::Float64
    samples_per_second::Float64
end

function precompute_experimental_velocities_over_time(
    grid, geo,
    #cell_speeds_of_sound::Vector{Float64}, we'll add this in later once we start using two-loop optimization
    transducer_ids::Vector{Int},
    transducer_opposing_pairs_ids::Vector{Tuple{Int, Int}},
    transducer_projection_intersected_cells,
    transducer_projection_cell_distances,
    transducer_projection_center_ray_distances,
    transducer_projection_center_ray_intersected_cells,
    transducer_projection_center_ray_distances_through_cells,
    projection_unit_vectors_to_cell_centers,
    cell_projection_counts,
    experimental_sonic_data::SonicData,
    window_size_indices::Int;
    reltol = 1.5
)
    n_cells = length(grid.cells)
    
    # Compute grid-relative tolerance if not explicitly specified
    if isnothing(tolerance)
        avg_cell_volume = sum(geo.cell_volumes) / n_cells
        tolerance = reltol * cbrt(avg_cell_volume) #default reltol is 1.5
    end
    
    cell_speeds_of_sound, speeds_of_sound_interps, speed_of_sound_uncertainty_interps, speed_of_sound_timestamps = precompute_experimental_speeds_of_sound(
        transducer_opposing_pairs_ids,
        experimental_sonic_data,
        transducer_projection_center_ray_intersected_cells,
        transducer_projection_center_ray_distances_through_cells,
        n_cells, geo
    )
    
    # Populate beam_directions for each cell from the transducer projections
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

    # Velocity timesteps come from self-to-self echo pairs
    # Use the first transducer's timestamps as the reference time grid for velocities
    # (all transducers should fire at the same self-to-self rate)
    n_velocity_timesteps = length(experimental_sonic_data.self_to_self_send_echo_timestamps[1])

    experimental_discrete_vel_x = [zeros(Float64, n_cells) for _ in 1:n_velocity_timesteps]
    experimental_discrete_vel_y = [zeros(Float64, n_cells) for _ in 1:n_velocity_timesteps]
    experimental_discrete_vel_z = [zeros(Float64, n_cells) for _ in 1:n_velocity_timesteps]
    cell_velocity_uncertainties = Vector{Vector{Vector{Float64}}}(undef, n_velocity_timesteps)

    # Storage for backscatter and turbulence (Features 1 & 2)
    experimental_discrete_backscatter = [zeros(Float64, n_cells) for _ in 1:n_velocity_timesteps]
    experimental_discrete_turbulence = [zeros(Float64, n_cells) for _ in 1:n_velocity_timesteps]
    
    # Loop over each velocity timestep
    for timestamp_idx in 1:n_velocity_timesteps
        cell_update_counts = zeros(Int, n_cells)
        beam_velocities = [zeros(Float64, cell_projection_counts[cell_id]) for cell_id in 1:n_cells]

        # Feature 1 & 2: accumulators for backscatter and turbulence
        cell_backscatter_sum = zeros(Float64, n_cells)
        cell_backscatter_counts = zeros(Int, n_cells)
        cell_turbulence_sum_sq = zeros(Float64, n_cells)
        cell_turbulence_counts = zeros(Int, n_cells)
        
        # Find the speed-of-sound snapshot closest to this velocity timestamp
        # Use the average of send and echo times as the representative time for this measurement
        reference_transducer_id = transducer_ids[1]
        velocity_time = (experimental_sonic_data.self_to_self_send_echo_timestamps[reference_transducer_id][timestamp_idx].send +
                         experimental_sonic_data.self_to_self_send_echo_timestamps[reference_transducer_id][timestamp_idx].echo) / 2
        
        # Find closest speed-of-sound timestamp
        sos_t_idx = argmin(abs.(speed_of_sound_timestamps .- velocity_time))
        
        # Process each transducer beam
        for transducer_id in transducer_ids
            # Get center ray for this transducer's beam
            ray_cells = transducer_projection_center_ray_intersected_cells[transducer_id]
            ray_distances = transducer_projection_center_ray_distances_through_cells[transducer_id]
            
            # Build mapping using speed of sound at the closest timestamp
            time_to_distance_map = build_time_to_distance_map(cell_speeds_of_sound[sos_t_idx], ray_cells, ray_distances)
            
            # Extract velocities from this transducer's echo pair at this timestamp
            pulse_1 = experimental_sonic_data.self_to_self_echo_pairs[transducer_id][timestamp_idx][1]
            pulse_2 = experimental_sonic_data.self_to_self_echo_pairs[transducer_id][timestamp_idx][2]
            
            window_times, window_partial_velocities, window_backscatter_amplitudes, window_turbulence_widths = extract_velocities_from_data(
                pulse_1, pulse_2, experimental_sonic_data.pulse_repetition_interval, experimental_sonic_data.samples_per_second, window_size_indices
            )
            
            # Map experimental velocities to cells
            proj_cells = transducer_projection_intersected_cells[transducer_id]
            proj_dists = transducer_projection_cell_distances[transducer_id]
            
            map_experimental_velocities_to_cells!(
                beam_velocities, cell_update_counts, cell_speeds_of_sound[sos_t_idx],
                window_times, window_partial_velocities, time_to_distance_map,
                proj_cells, proj_dists, tolerance
            )

            # Feature 1: Map backscatter amplitudes to cells
            map_experimental_backscatter_to_cells!(
                cell_backscatter_sum, cell_backscatter_counts,
                window_times, window_backscatter_amplitudes, time_to_distance_map,
                proj_cells, proj_dists, tolerance
            )

            # Feature 2: Map turbulence widths to cells
            map_experimental_turbulence_to_cells!(
                cell_turbulence_sum_sq, cell_turbulence_counts,
                cell_speeds_of_sound[sos_t_idx],
                window_times, window_turbulence_widths, time_to_distance_map,
                proj_cells, proj_dists, tolerance
            )
        end
        
        # Solve least squares system for this timestep
        cell_velocity_uncertainties[timestamp_idx] = finalize_velocities!(
            experimental_discrete_vel_x[timestamp_idx], experimental_discrete_vel_y[timestamp_idx], experimental_discrete_vel_z[timestamp_idx], 
            beam_directions, beam_velocities, cell_update_counts
        )
        
        # Patch unresolved cells
        patch_unaccounted_experimental_velocities!(
            experimental_discrete_vel_x[timestamp_idx], 
            experimental_discrete_vel_y[timestamp_idx], 
            experimental_discrete_vel_z[timestamp_idx], 
            geo
        )

        # Feature 1: Finalize and patch backscatter
        finalize_backscatter!(
            experimental_discrete_backscatter[timestamp_idx],
            cell_backscatter_sum, cell_backscatter_counts
        )
        patch_unaccounted_scalar_field!(experimental_discrete_backscatter[timestamp_idx], geo)

        # Feature 2: Finalize and patch turbulence intensity
        finalize_turbulence!(
            experimental_discrete_turbulence[timestamp_idx],
            cell_turbulence_sum_sq, cell_turbulence_counts
        )
        patch_unaccounted_scalar_field!(experimental_discrete_turbulence[timestamp_idx], geo)
    end

    # Build the time axis from the reference transducer's self-to-self timestamps
    reference_transducer_id = transducer_ids[1]
    velocity_timestamps = [
        (experimental_sonic_data.self_to_self_send_echo_timestamps[reference_transducer_id][t].send +
         experimental_sonic_data.self_to_self_send_echo_timestamps[reference_transducer_id][t].echo) / 2
        for t in 1:n_velocity_timesteps
    ]

    # Feature 4: Compute bubble passage frequency from backscatter time series
    cell_passage_frequencies = compute_bubble_passage_frequencies(
        experimental_discrete_backscatter, velocity_timestamps, n_cells
    )

    # Build per-cell velocity interpolations over the velocity timestamps
    experimental_vel_x_interp = Vector{LinearInterpolation}(undef, n_cells)
    experimental_vel_y_interp = Vector{LinearInterpolation}(undef, n_cells)
    experimental_vel_z_interp = Vector{LinearInterpolation}(undef, n_cells)

    # Feature 1 & 2 & 4: interpolations for backscatter, turbulence, passage frequency
    backscatter_intensity_interps = Vector{LinearInterpolation}(undef, n_cells)
    turbulence_intensity_interps = Vector{LinearInterpolation}(undef, n_cells)
    bubble_passage_frequency_interps = Vector{LinearInterpolation}(undef, n_cells)

    for cell_id in 1:n_cells
        cell_vel_x_history = [experimental_discrete_vel_x[t][cell_id] for t in 1:n_velocity_timesteps]
        cell_vel_y_history = [experimental_discrete_vel_y[t][cell_id] for t in 1:n_velocity_timesteps]
        cell_vel_z_history = [experimental_discrete_vel_z[t][cell_id] for t in 1:n_velocity_timesteps]

        experimental_vel_x_interp[cell_id] = LinearInterpolation(cell_vel_x_history, velocity_timestamps)
        experimental_vel_y_interp[cell_id] = LinearInterpolation(cell_vel_y_history, velocity_timestamps)
        experimental_vel_z_interp[cell_id] = LinearInterpolation(cell_vel_z_history, velocity_timestamps)

        # Feature 1: Backscatter intensity per cell over time
        cell_backscatter_history = [experimental_discrete_backscatter[t][cell_id] for t in 1:n_velocity_timesteps]
        backscatter_intensity_interps[cell_id] = LinearInterpolation(cell_backscatter_history, velocity_timestamps)

        # Feature 2: Turbulence intensity per cell over time
        cell_turbulence_history = [experimental_discrete_turbulence[t][cell_id] for t in 1:n_velocity_timesteps]
        turbulence_intensity_interps[cell_id] = LinearInterpolation(cell_turbulence_history, velocity_timestamps)

        # Feature 4: Bubble passage frequency per cell over time
        cell_passage_history = [cell_passage_frequencies[t][cell_id] for t in 1:n_velocity_timesteps]
        bubble_passage_frequency_interps[cell_id] = LinearInterpolation(cell_passage_history, velocity_timestamps)
    end

    cell_vel_x_uncertainty_interps = Vector{LinearInterpolation}(undef, n_cells)
    cell_vel_y_uncertainty_interps = Vector{LinearInterpolation}(undef, n_cells)
    cell_vel_z_uncertainty_interps = Vector{LinearInterpolation}(undef, n_cells)

    for cell_id in 1:n_cells 
        cell_vel_x_uncertainty_interps[cell_id] = LinearInterpolation(
            [cell_velocity_uncertainties[t][cell_id][1] for t in 1:n_velocity_timesteps],
            velocity_timestamps
        )
        cell_vel_y_uncertainty_interps[cell_id] = LinearInterpolation(
            [cell_velocity_uncertainties[t][cell_id][2] for t in 1:n_velocity_timesteps],
            velocity_timestamps
        )
        cell_vel_z_uncertainty_interps[cell_id] = LinearInterpolation(
            [cell_velocity_uncertainties[t][cell_id][3] for t in 1:n_velocity_timesteps],
            velocity_timestamps
        )
    end

    return experimental_vel_x_interp, experimental_vel_y_interp, experimental_vel_z_interp, 
           cell_vel_x_uncertainty_interps, cell_vel_y_uncertainty_interps, cell_vel_z_uncertainty_interps, 
           speeds_of_sound_interps, speeds_of_sound_uncertainty_interps,
           backscatter_intensity_interps, turbulence_intensity_interps, bubble_passage_frequency_interps
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 3: Spectral attenuation precomputation
# Processes transducer-to-transducer received amplitude signals to extract
# frequency-dependent attenuation (normalized spectral shape).
# This is separate from the velocity pipeline because it uses different data
# (transducer-to-transducer paths, not self-to-self echo pairs).
# ═══════════════════════════════════════════════════════════════════════════════

function precompute_experimental_spectral_attenuation(
    transducer_ids::Vector{Int},
    transducer_opposing_pairs_ids::Vector{Tuple{Int, Int}},
    experimental_sonic_data::SonicData,
    transducer_center_frequency::Float64;
    n_frequency_bins = 10,
    bandwidth_factor = 0.25
)
    n_transducers = length(transducer_ids)

    # Determine the frequency bins from the first available measurement
    first_pair = transducer_opposing_pairs_ids[1]
    first_signal = experimental_sonic_data.transducer_to_transducer_receive_amplitudes[first_pair[1], first_pair[2]][1]
    spectral_frequencies, _ = extract_normalized_spectral_shape(
        first_signal, experimental_sonic_data.samples_per_second, transducer_center_frequency;
        n_frequency_bins, bandwidth_factor
    )
    n_freq_bins = length(spectral_frequencies)

    if n_freq_bins == 0
        return spectral_frequencies, [Vector{LinearInterpolation}() for _ in 1:n_transducers, _ in 1:n_transducers]
    end

    spectral_attenuation_interps = [Vector{LinearInterpolation}(undef, n_freq_bins) for _ in 1:n_transducers, _ in 1:n_transducers]

    for (transducer_a_id, transducer_b_id) in transducer_opposing_pairs_ids
        for (src_id, dst_id) in [(transducer_a_id, transducer_b_id), (transducer_b_id, transducer_a_id)]
            amplitude_signals = experimental_sonic_data.transducer_to_transducer_receive_amplitudes[src_id, dst_id]
            timestamps_data = experimental_sonic_data.transducer_to_transducer_send_receive_timestamps[src_id, dst_id]

            n_measurements = length(amplitude_signals)
            measurement_timestamps = [timestamps_data[t].send for t in 1:n_measurements]

            # Extract spectral shape at each timestamp
            spectral_shapes = [zeros(Float64, n_freq_bins) for _ in 1:n_measurements]

            for t_idx in 1:n_measurements
                _, shape = extract_normalized_spectral_shape(
                    amplitude_signals[t_idx], experimental_sonic_data.samples_per_second, transducer_center_frequency;
                    n_frequency_bins, bandwidth_factor
                )
                if length(shape) == n_freq_bins
                    spectral_shapes[t_idx] = shape
                end
            end

            # Build per-frequency-bin interpolations over time
            for f_idx in 1:n_freq_bins
                freq_history = [spectral_shapes[t][f_idx] for t in 1:n_measurements]
                spectral_attenuation_interps[src_id, dst_id][f_idx] = LinearInterpolation(freq_history, measurement_timestamps)
            end
        end
    end

    return spectral_frequencies, spectral_attenuation_interps
end