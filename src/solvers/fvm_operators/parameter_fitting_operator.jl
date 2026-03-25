function methylene_blue_diffuion_parameter_fitting_f!(
    du_vec, u_vec, p, t, 
    state_data, state_time,

    p_axes,

    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    connection_groups, controller_groups, patch_groups, region_groups, 
    
    properties,

    du_diff_cache, u_diff_cache,
    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes
)
    du_cache_vec = get_tmp(du_diff_cache, first(u_vec) + first(p))
    u_cache_vec = get_tmp(u_diff_cache, first(u_vec) + first(p))

    #println(typeof(du_cache_vec))

    #do you know what is fucking crazy?
    #the only reason this worked in the past was because I zeroed out the u_cache_vec
    #otherwise, it just returns #undef for everything
    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = create_views_inline(du_cache_vec, du_cache_axes)
    u_cache_nt = create_views_inline(u_cache_vec, u_cache_axes)

    du_vec .= 0.0
    du_nt = create_views_inline(du_vec, du_proto_axes)
    u_nt = create_views_inline(u_vec, u_proto_axes)

    du = (; du_nt..., du_cache_nt...)
    u = (; properties..., u_nt..., u_cache_nt...) 
    #remember, the right-most fields overwrite other fields with the same name

    u.rho .= (u_cache_nt.rho .+ properties.rho)
    #the reason we add it to u_cache_nt is because u_cache_nt is a DiffCache

    #get state variables from experimental data from an interp
    temp_at_t = state_data.temp(t)
    u.temp .= temp_at_t
    #we don't need a DiffCache here

    p_named = create_views_inline(p, p_axes)

    #get optimized parameters from p 
    #since single numbers are stored as one element vectors in create_views_inline, we have to get the first 
    u.diffusion_pre_exponential_factor .= p_named.diffusion_pre_exponential_factor[1]
    u.diffusion_activation_energy .= p_named.diffusion_activation_energy[1]

    #println(u.diffusion_pre_exponential_factor[1], ", ", u.diffusion_activation_energy[1])
    
    #Connection Loops
    for conn in connection_groups
        solve_connection_group!(
            du, u, 
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        )
    end

    #patches must be solved before regions because the regions contain the capacities
    for patch in patch_groups
        solve_patch_group!(
            du, u,
            patch.patch_function!, patch.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            cell_volumes
        )
    end

    for reg in region_groups
        solve_region_group!(
            du, u,
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end

    #custom method for a well-mixed assumption
    for reg in region_groups
        if reg.name == "surrounding_fluid"
            total_methylene_blue_dm_dt = 0.0
            total_water_dm_dt = 0.0
            total_reservoir_mass = 0.0
            total_reservoir_volume = 0.0

            for cell_id in reg.region_cells
                total_methylene_blue_dm_dt += du.mass_fractions.methylene_blue[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]
                total_water_dm_dt += du.mass_fractions.water[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]

                total_reservoir_mass += u.rho[cell_id] * cell_volumes[cell_id]
                total_reservoir_volume += cell_volumes[cell_id]
            end

            well_mixed_methylene_blue_dt = total_methylene_blue_dm_dt / total_reservoir_mass
            well_mixed_water_dt = total_water_dm_dt / total_reservoir_mass
            
            for cell_id in reg.region_cells
                du.mass_fractions.methylene_blue[cell_id] = well_mixed_methylene_blue_dt
                du.mass_fractions.water[cell_id] = well_mixed_water_dt
            end
        end
        if reg.name == "dialysis_tubing_interior"
            total_methylene_blue_dm_dt = 0.0
            total_water_dm_dt = 0.0
            total_reservoir_mass = 0.0
            total_reservoir_volume = 0.0

            for cell_id in reg.region_cells
                total_methylene_blue_dm_dt += du.mass_fractions.methylene_blue[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]
                total_water_dm_dt += du.mass_fractions.water[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]

                total_reservoir_mass += u.rho[cell_id] * cell_volumes[cell_id]
                total_reservoir_volume += cell_volumes[cell_id]
            end

            well_mixed_methylene_blue_dt = total_methylene_blue_dm_dt / total_reservoir_mass
            well_mixed_water_dt = total_water_dm_dt / total_reservoir_mass
            
            for cell_id in reg.region_cells
                du.mass_fractions.methylene_blue[cell_id] = well_mixed_methylene_blue_dt
                du.mass_fractions.water[cell_id] = well_mixed_water_dt
            end
        end
    end
end

#instead of updating rho in each flux and internal physics, we could do an initial property retrieval loop
#=
@batch for cell_id in eachindex(cell_volumes)
    mw_avg!(u, cell_id, molecular_weights, mw_avg_cache)
    rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)
end
=#