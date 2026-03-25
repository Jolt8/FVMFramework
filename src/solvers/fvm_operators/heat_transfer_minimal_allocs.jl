function heat_transfer_f_test!(
    du_vec, u_vec, p, t, 

    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    connection_groups, controller_groups, region_groups, patch_groups,

    du_virtual_axes, u_virtual_axes,
    du_diff_cache, u_diff_cache,
    properties_vec,
)
    
    du_vec .= 0.0
    
    if (first(u_vec) + first(p)) isa SparseConnectivityTracer.Dual{Float64}
        u = VirtualFVMArray((u_vec, (get_tmp(u_diff_cache, first(u_vec) + first(p)) .= 0.0), properties_vec), u_virtual_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache, first(u_vec) + first(p)) .= 0.0)), du_virtual_axes)
    else
        u = VirtualFVMArray((u_vec, get_tmp(u_diff_cache, first(u_vec) + first(p)), properties_vec), u_virtual_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache, first(u_vec) + first(p)) .= 0.0)), du_virtual_axes)
    end
    
    #Connection Loops
    for conn in connection_groups
        solve_connection_group!(
            du, u, 
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        )
    end

    for cont in controller_groups
        solve_controller_group!(
            du, u, cont.id, cont.controller, cont.controller_function!, cont.monitored_cells, cont.affected_cells,
            cell_volumes
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

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    for reg in region_groups
        solve_region_group!(
            du, u, 
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end

    #for reg in region_groups
    #    debug_region!(du, u, reg)
    #end
end