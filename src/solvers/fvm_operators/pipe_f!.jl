function pipe_f!(
    du_vec, u_vec, p, t, 

    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    connection_groups, controller_groups, region_groups, patch_groups,

    du_virtual_axes, u_virtual_axes,
    du_diff_cache_vec, u_diff_cache_vec,
    properties_vec,
)
    
    du_vec .= 0.0
    
    if (first(u_vec) + first(p)) isa SparseConnectivityTracer.Dual{Float64}
        u = VirtualFVMArray((u_vec, (get_tmp(u_diff_cache_vec, first(u_vec) + first(p)) .= 0.0), properties_vec), u_virtual_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache_vec, first(u_vec) + first(p)) .= 0.0)), du_virtual_axes)
    else
        u = VirtualFVMArray((u_vec, get_tmp(u_diff_cache_vec, first(u_vec) + first(p)), properties_vec), u_virtual_axes)
        du = VirtualFVMArray((du_vec, get_tmp(du_diff_cache_vec, first(u_vec) + first(p))), du_virtual_axes)
    end

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

    @show du.temp[:]
    @show du.pressure[:]

    for patch in patch_groups
        solve_patch_group!(
            du, u,
            patch.patch_function!, patch.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            cell_volumes
        )
    end

    @show du.temp[:]
    @show du.pressure[:]
    
    
    for reg in region_groups
        solve_region_group!(
            du, u,
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end
    
    @show du.temp[:]
    @show du.pressure[:]
end