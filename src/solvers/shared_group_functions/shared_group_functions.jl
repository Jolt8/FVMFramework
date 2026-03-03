function solve_connection_group!(
    du, u,
    flux!::F, cell_neighbors,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
) where {F}

    #=@batch=# for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx) in neighbor_list
            flux!(
                du, u,
                idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
            )
        end
    end
end

function solve_controller_group!(
    du, u,
    controller::C, controller_id,
    control!::F, monitored_cells, affected_cells,
    cell_volumes
) where {C, F}
    control!(
        du, u, controller, controller_id,
        monitored_cells, affected_cells,
        cell_volumes
    )
end

function solve_region_group!(
    du, u,
    internal_physics!::F, region_cells,
    cell_volumes
) where {F}
    #=@batch=# for cell_id in region_cells
        internal_physics!(
            du, u, cell_id,
            cell_volumes[cell_id]
        )
    end
end

function solve_patch_group!(
    du, u,
    patch_physics!::F, cell_neighbors,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
) where {F}
    #=@batch=# for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx) in neighbor_list
            patch_physics!(
                du, u,
                idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                cell_volumes
            )
        end
    end
end