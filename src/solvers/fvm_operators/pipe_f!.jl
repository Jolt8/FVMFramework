function pipe_f!(
    du_vec, u_vec, p, t, 

    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    connection_groups, controller_groups, region_groups, patch_groups,

    du_virtual_axes, u_virtual_axes,
    du_diff_cache, u_diff_cache,
    properties_vec,
)   
    du, u = unpack_fvm_state(du_vec, u_vec, p, t, du_virtual_axes, u_virtual_axes, du_diff_cache, u_diff_cache, properties_vec)

    #1st cell neighbors: 
        #- (1, [(2, 3)])
        #- (cell_id, [(neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor)])
    #2:n_cells-1 neigbors: 
        #- (2, [(3, 3), (1, 5)])
        #- (cell_id, [(neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor), (neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor)])
    #last cell neighbors: 
        #- (20, [(19, 5)])
        #- (cell_id, [(neighbor_cell_id, current_cell's_face_index_pointing_to_neighbor)])
    
    for cell_id in 1:length(cell_volumes)-1
        du.mass_face[cell_id, 3] -= u.pipe_mass_flow[cell_id]
        du.mass_face[cell_id + 1, 5] += u.pipe_mass_flow[cell_id]
    end

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