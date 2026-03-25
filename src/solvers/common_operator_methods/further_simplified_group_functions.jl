
function solve_connection_groups!(du, u, geo, system)
    for conn in system.connection_groups
        solve_connection_group!(
            du, u,
            conn.flux_function!, conn.cell_neighbors,
            geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
        )
    end
end

#I sometimes wonder if controllers would be niche enough to not warrant their own group.
#I'll keep them for now

function solve_controller_groups!(du, u, geo, system)
    for cont in system.controller_groups 
        solve_controller_group!(
            du, u, cont.id, cont.controller, cont.controller_function!, cont.monitored_cells, cont.affected_cells,
            geo.cell_volumes
        )
    end
end

function solve_patch_groups!(du, u, geo, system)
    for patch in system.patch_groups
        solve_patch_group!(
            du, u,
            patch.patch_function!, patch.cell_neighbors,
            geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
            geo.cell_volumes
        )
    end
end

function solve_region_groups!(du, u, geo, system)
    for reg in system.region_groups
        solve_region_group!(
            du, u,
            reg.region_function!, reg.region_cells,
            geo.cell_volumes
        )
    end
end

function default_order_solve_all_groups!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system)
end