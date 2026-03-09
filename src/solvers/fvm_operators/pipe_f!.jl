function pipe_f!(
    du_vec, u_vec, p, t, 

    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    connection_groups, controller_groups, region_groups, patch_groups,
    
    properties,

    du_diff_cache_vec, u_diff_cache_vec,
    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes
)
    du_cache_vec = get_tmp(du_diff_cache_vec, u_vec)
    u_cache_vec = get_tmp(u_diff_cache_vec, u_vec)

    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = create_views_inline(du_cache_vec, du_cache_axes)
    u_cache_nt = create_views_inline(u_cache_vec, u_cache_axes)

    du_vec .= 0.0
    du_nt = create_views_inline(du_vec, du_proto_axes)
    u_nt = create_views_inline(u_vec, u_proto_axes)

    du = (; du_nt..., du_cache_nt...)
    #u = (; u_nt..., u_cache_nt..., properties...)
    u = (; properties..., u_nt..., u_cache_nt...)
    #remember, the right-most fields overwrite other fields with the same name

    u.rho .= (u_cache_nt.rho .+ properties.rho)

    #1st cell neighbors: 
        #- (1, [(2, 3)])
        #- (cell_id, [(neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor)])
    #2:n_cells-1 neigbors: 
        #- (2, [(3, 3), (1, 5)])
        #- (cell_id, [(neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor), (neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor)])
    #last cell neighbors: 
        #- (20, [(19, 5)])
        #- (cell_id, [(neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor)])

    #du.mass_face[1][3] -= u.pipe_mass_flow[1] #inlet
    
    for cell_id in 1:length(cell_volumes)-1
        du.mass_face[cell_id][3] -= u.pipe_mass_flow[cell_id]
        du.mass_face[cell_id + 1][5] += u.pipe_mass_flow[cell_id]
    end

    #du.mass_face[end][5] += u.pipe_mass_flow[end] #outlet

    for conn in connection_groups
        solve_connection_group!(
            du, u,
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        )
    end

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
end