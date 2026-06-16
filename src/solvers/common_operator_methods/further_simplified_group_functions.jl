
function solve_connection_groups!(du, u, p, t, geo, system)
    for conn in system.connection_groups
        solve_connection_group!(
            du, u, p, t,
            conn.flux_function!, conn.cell_neighbors,
            geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
            geo.cell_volumes
        )
    end
end

#I sometimes wonder if controllers would be niche enough to not warrant their own group.
#I'll keep them for now

function solve_controller_groups!(du, u, p, t, geo, system)
    for cont in system.controller_groups 
        solve_controller_group!(
            du, u, p, t, 
            cont.id, cont.controller, cont.controller_function!, cont.monitored_cells, cont.affected_cells,
            geo.cell_volumes
        )
    end
end

function solve_patch_groups!(du, u, p, t, geo, system)
    for patch in system.patch_groups
        solve_patch_group!(
            du, u, p, t, 
            patch.patch_function!, patch.cell_neighbors,
            geo.cell_face_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
            geo.cell_volumes
        )
    end
end

function update_region_groups!(du, u, p, t, geo, system)
    for reg in system.region_groups
        update_region_group!(
            du, u, p, t,
            reg.property_update_function!, reg.region_cells,
            geo.cell_volumes
        )
    end
end

function solve_region_groups!(du, u, p, t, geo, system)
    for reg in system.region_groups
        solve_region_group!(
            du, u, p, t,
            reg.region_function!, reg.region_cells,
            geo.cell_volumes
        )
    end
end

function default_order_solve_all_groups!(du, u, p, t, geo, system)
    update_region_groups!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, p, t, geo, system)
    solve_controller_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end