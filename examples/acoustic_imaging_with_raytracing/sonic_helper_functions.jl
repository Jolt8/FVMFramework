
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

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 1: Backscatter amplitude simulation
# The pulse-echo signal carries envelope amplitude at each depth, which is
# proportional to gas_holdup × bubble_cross_section, attenuated by the
# round-trip path through the medium.
# ═══════════════════════════════════════════════════════════════════════════════

function calculate_simulated_backscatter_amplitudes(
    u, transducer_id,
    center_ray_intersected_cells, center_ray_distances_through_cells
)
    ray_cells = center_ray_intersected_cells[transducer_id]
    ray_dists = center_ray_distances_through_cells[transducer_id]

    n_ray_cells = length(ray_cells)
    cell_backscatter_amplitudes = zeros(Float64, n_ray_cells)

    cumulative_one_way_attenuation = 0.0

    for (idx, cell_id) in enumerate(ray_cells)
        # Local scatterer reflectivity ∝ gas_holdup × geometric bubble cross-section
        scatterer_cross_section = pi * u.bubble_radii[cell_id]^2
        local_reflectivity = u.gas_holdup[cell_id] * scatterer_cross_section

        # Round-trip attenuation (sound travels out to the cell and back)
        round_trip_attenuation = exp(-2.0 * cumulative_one_way_attenuation)

        cell_backscatter_amplitudes[idx] = local_reflectivity * round_trip_attenuation

        # Accumulate one-way attenuation for the next cell
        local_attenuation_coefficient = u.gas_holdup[cell_id] * u.bubble_radii[cell_id]
        cumulative_one_way_attenuation += local_attenuation_coefficient * ray_dists[idx]
    end

    return ray_cells, cell_backscatter_amplitudes
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 3: Frequency-dependent attenuation → bubble size distribution
# Uses a simplified Commander-Prosperetti model for single-bubble extinction
# cross-section, integrated over a log-normal bubble size distribution.
#
# Acoustic attenuation vs. frequency follows different scaling laws:
#   ka << 1 (thermal/viscous regime): α ∝ f⁴ a³
#   ka ~ 1  (resonance regime):       α peaks near Minnaert frequency
#   ka >> 1 (geometric regime):        α ∝ a²
# ═══════════════════════════════════════════════════════════════════════════════

function minnaert_resonance_frequency(bubble_radius, ambient_pressure, liquid_density; polytropic_exponent = 1.4)
    # Minnaert frequency: f_res = (1 / 2πa) √(3γP₀ / ρ_l)
    return (1.0 / (2.0 * pi * bubble_radius)) * sqrt(3.0 * polytropic_exponent * ambient_pressure / liquid_density)
end

function single_bubble_damping(frequency, bubble_radius, liquid_density, liquid_viscosity, speed_of_sound)
    ω = 2.0 * pi * frequency

    # Viscous damping
    δ_viscous = 4.0 * liquid_viscosity / (liquid_density * ω * bubble_radius^2)

    # Radiation damping
    δ_radiation = ω * bubble_radius / speed_of_sound

    # Thermal damping (simplified constant approximation, ~25% of viscous for air-water)
    δ_thermal = 0.25 * δ_viscous

    return δ_viscous + δ_radiation + δ_thermal
end

function single_bubble_extinction_cross_section(
    frequency, bubble_radius,
    ambient_pressure, liquid_density, liquid_viscosity, speed_of_sound;
    polytropic_exponent = 1.4
)
    f_res = minnaert_resonance_frequency(bubble_radius, ambient_pressure, liquid_density; polytropic_exponent)
    δ_total = single_bubble_damping(frequency, bubble_radius, liquid_density, liquid_viscosity, speed_of_sound)

    # Commander-Prosperetti extinction cross-section:
    # σ_ext = 4πa² / ((f_res²/f² - 1)² + δ_total²)
    f_ratio = f_res / frequency
    σ_ext = 4.0 * pi * bubble_radius^2 / ((f_ratio^2 - 1.0)^2 + δ_total^2)

    return σ_ext
end

function lognormal_pdf(a, a_median, σ_ln)
    # Log-normal PDF where ln(a) ~ N(ln(a_median), σ_ln²)
    if a <= 0.0
        return 0.0
    end
    ln_a = log(a)
    ln_median = log(a_median)
    return exp(-(ln_a - ln_median)^2 / (2.0 * σ_ln^2)) / (a * σ_ln * sqrt(2.0 * pi))
end

function lognormal_mean_extinction_cross_section(
    frequency, a_median, σ_ln,
    ambient_pressure, liquid_density, liquid_viscosity, speed_of_sound;
    polytropic_exponent = 1.4,
    n_quadrature_points = 20
)
    # Integrate σ_ext(f, a) × p(a; a_median, σ_ln) da via trapezoidal rule in log-space
    # Integration bounds: ±4σ around the median in ln(a)-space
    ln_median = log(a_median)
    ln_min = ln_median - 4.0 * σ_ln
    ln_max = ln_median + 4.0 * σ_ln

    d_ln = (ln_max - ln_min) / n_quadrature_points

    integral = 0.0

    for i in 0:n_quadrature_points
        ln_a = ln_min + i * d_ln
        a = exp(ln_a)

        σ_ext = single_bubble_extinction_cross_section(
            frequency, a,
            ambient_pressure, liquid_density, liquid_viscosity, speed_of_sound;
            polytropic_exponent
        )

        pdf_val = lognormal_pdf(a, a_median, σ_ln)

        # Trapezoidal weighting
        weight = (i == 0 || i == n_quadrature_points) ? 0.5 : 1.0

        # Change of variables: da = a × d(ln a)
        integral += weight * σ_ext * pdf_val * a * d_ln
    end

    return integral
end

function calculate_frequency_dependent_attenuation_coefficient(
    frequency, gas_holdup, bubble_radii_median, bubble_radii_sigma,
    ambient_pressure, liquid_density, liquid_viscosity, speed_of_sound;
    polytropic_exponent = 1.4
)
    # Number density from gas holdup:
    # For log-normal: E[a³] = a_median³ × exp(9σ_ln²/2)
    mean_a3 = bubble_radii_median^3 * exp(4.5 * bubble_radii_sigma^2)
    V_avg = (4.0 / 3.0) * pi * mean_a3
    n_density = gas_holdup / max(V_avg, 1e-30)

    if bubble_radii_sigma > 0.01
        # Distribution-averaged cross-section
        mean_σ_ext = lognormal_mean_extinction_cross_section(
            frequency, bubble_radii_median, bubble_radii_sigma,
            ambient_pressure, liquid_density, liquid_viscosity, speed_of_sound;
            polytropic_exponent
        )
    else
        # Narrow distribution → single-bubble approximation
        mean_σ_ext = single_bubble_extinction_cross_section(
            frequency, bubble_radii_median,
            ambient_pressure, liquid_density, liquid_viscosity, speed_of_sound;
            polytropic_exponent
        )
    end

    # Attenuation coefficient α (Np/m)
    return n_density * mean_σ_ext
end

function calculate_simulated_total_attenuation_per_frequency(
    u, origin_transducer_id, destination_transducer_id,
    ray_map_intersected_cells, ray_map_distances_through_cells,
    frequencies;
    ambient_pressure = 101325.0,
    liquid_viscosity = 1e-3,
    polytropic_exponent = 1.4
)
    intersected_cells = ray_map_intersected_cells[origin_transducer_id, destination_transducer_id]
    intersected_cell_distances = ray_map_distances_through_cells[origin_transducer_id, destination_transducer_id]

    n_freqs = length(frequencies)
    total_attenuation = zeros(Float64, n_freqs)

    for f_idx in 1:n_freqs
        for (cell_idx, cell_id) in enumerate(intersected_cells)
            α = calculate_frequency_dependent_attenuation_coefficient(
                frequencies[f_idx],
                u.gas_holdup[cell_id],
                u.bubble_radii[cell_id],
                u.bubble_radii_sigma[cell_id],
                ambient_pressure,
                u.fluid_rho[cell_id],
                liquid_viscosity,
                u.speed_of_sound[cell_id];
                polytropic_exponent
            )

            total_attenuation[f_idx] += α * intersected_cell_distances[cell_idx]
        end
    end

    return total_attenuation
end

# ═══════════════════════════════════════════════════════════════════════════════
# Feature 4: Bubble passage frequency / number density
# The passage rate of bubbles through a cell is:
#   f_passage = n_density × |v_gas| × A_measurement
# where n_density = gas_holdup / bubble_volume.
# Combined with bubble_radii and gas_velocity, this directly closes the
# constraint between gas_holdup and bubble_radii.
# ═══════════════════════════════════════════════════════════════════════════════

function calculate_simulated_bubble_passage_frequency(u, cell_id, measurement_cross_section)
    # Number density from gas holdup and mean bubble volume
    # For log-normal: E[a³] = a_median³ × exp(9σ_ln²/2)
    mean_a3 = u.bubble_radii[cell_id]^3 * exp(4.5 * u.bubble_radii_sigma[cell_id]^2)
    bubble_volume = (4.0 / 3.0) * pi * mean_a3
    n_density = u.gas_holdup[cell_id] / max(bubble_volume, 1e-30)

    # Passage frequency = number density × |gas velocity| × measurement cross-section
    gas_speed = abs(u.local_gas_drift_velocity[cell_id])

    passage_frequency = n_density * gas_speed * measurement_cross_section

    return passage_frequency, n_density
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