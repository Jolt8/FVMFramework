function heat_transfer_f!(
    du_flat, u_flat, p, t,
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,
    connection_groups, controller_groups, region_groups,
    du_fixed_vec, u_fixed_vec,
    state_axes, fixed_axes
)
    du_cv = ComponentVector(du_flat, state_axes)
    u_cv = ComponentVector(u_flat, state_axes)

    if eltype(du_flat) == Float64
        du_fixed_vec = get_tmp(du_fixed_vec, Float64)
        u_fixed_vec = get_tmp(u_fixed_vec, Float64)
    else
        du_fixed_vec = get_tmp(du_fixed_vec, ForwardDiff.Dual(1.0, 1.0))
        u_fixed_vec = get_tmp(u_fixed_vec, ForwardDiff.Dual(1.0, 1.0))
    end

    du_fixed_cv = ComponentVector(du_fixed_vec, fixed_axes)
    u_fixed_cv = ComponentVector(u_fixed_vec, fixed_axes)

    #du_merged = ComponentVector(merge(NamedTuple(du_cv), NamedTuple(du_fixed_cv)))
    #u_merged = ComponentVector(merge(NamedTuple(u_cv), NamedTuple(u_fixed_cv)))

    du_merged = ComponentVector(du_cv; NamedTuple(du_fixed_cv)...)
    u_merged = ComponentVector(u_cv; NamedTuple(u_fixed_cv)...)

    du_merged .= 0.0

    # Connection Loops
    for conn in connection_groups
        @batch for (idx_a, neighbor_list) in conn.cell_neighbors
            for (idx_b, face_idx) in neighbor_list
                conn.flux_function!(
                    du_merged, u_merged, idx_a, idx_b, face_idx,
                    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
                )
            end
        end
    end

    # Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    for reg in region_groups
        @batch for cell_id in reg.region_cells
            reg.region_function!(
                du_merged, u_merged, cell_id,
                cell_volumes
            )
        end
    end
    @views du_flat .= du_merged[1:length(du_flat)]
end