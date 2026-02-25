#=
function solve_connection_group!(
    du, u,
    flux!::F, neighbors,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
) where {F}

    @batch for (idx_a, neighbor_list) in neighbors
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
    @batch for cell_id in region_cells
        internal_physics!(
            du, u, cell_id,
            cell_volumes[cell_id]
        )
    end
end
=#

function heat_transfer_f_test!(
        du_vec, u_vec, p, t,
        
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        unconnected_cell_face_map, cell_face_areas, cell_face_normals,
        connection_groups, controller_groups, region_groups,

        properties,
        du_diff_cache_vec, u_diff_cache_vec,
        du_proto_axes, u_proto_axes,
        du_cache_axes, u_cache_axes
    )
    du_cache_vec = get_tmp(du_diff_cache_vec, u_vec)
    u_cache_vec = get_tmp(u_diff_cache_vec, u_vec)

    #do you know what is fucking crazy?
    #the only reason this worked in the past was because I zeroed out the caches
    #otherwise, it just returns #undef for everything
    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    #V1 
        # 93.499 ms, 573 allocations, 49.92 MiB (10000 cells)
        # 9.227 ms, 573 allocations, 5.01 MiB (1000 cells)
    # requires: du_vec .= Vector(du.temp)
    #=
    du_cache_nt = (; NamedTuple(ComponentArray(du_cache_vec, du_cache_axes))...)
    u_cache_nt = (; NamedTuple(ComponentArray(u_cache_vec, u_cache_axes))...)
    
    du_vec .= 0.0
    du_nt = (; NamedTuple(ComponentArray(du_vec, du_proto_axes))...)
    u_nt = (; NamedTuple(ComponentArray(u_vec, u_proto_axes))...)
    =#

    # V2 
        # 103.895 ms, 1278 allocations, 15.62 MiB (10000 cells)
        # 11.616 ms, 1278 allocations, 1.61 MiB
    # does not require: du_vec .= Vector(du.temp)
    #=
    du_cache_ca = ComponentArray(du_cache_vec, du_cache_axes)
    du_cache_nt = NamedTuple{propertynames(du_cache_ca)}(Tuple(getproperty(du_cache_ca, p) for p in propertynames(du_cache_ca)))

    u_cache_ca = ComponentArray(u_cache_vec, u_cache_axes)
    u_cache_nt = NamedTuple{propertynames(u_cache_ca)}(Tuple(getproperty(u_cache_ca, p) for p in propertynames(u_cache_ca)))
    
    du_vec .= 0.0
    du_ca = ComponentArray(du_vec, du_proto_axes)
    du_nt = NamedTuple{propertynames(du_ca)}(Tuple(getproperty(du_ca, p) for p in propertynames(du_ca)))
    
    u_ca = ComponentArray(u_vec, u_proto_axes)
    u_nt = NamedTuple{propertynames(u_ca)}(Tuple(getproperty(u_ca, p) for p in propertynames(u_ca)))
    =#

    #V3 (Raw) 
        # 84.729 ms, 348 allocations, 15.58 MiB (10000 cells)
        # 7.330 ms, 348 allocations, 1.57 MiB
    # does not require: du_vec .= Vector(du.temp)
    #
    du_cache_nt = (heat = du_cache_vec,)
    u_cache_nt = (heat = u_cache_vec,)
    
    du_vec .= 0.0
    du_nt = (temp = du_vec,)
    u_nt = (temp = u_vec,)
    #

    du = (; du_nt..., du_cache_nt...)
    u = (; u_nt..., u_cache_nt..., properties...)
    
    #Connection Loops
    for conn in connection_groups
        solve_connection_group!(
            du, u, 
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        )
    end

    #=
    #Controller Loops
    for cont in controller_groups
        solve_controller_group!(
            du, u, cont.id, cont.controller, cont.controller_function!, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
    end
    =#

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    for reg in region_groups
        solve_region_group!(
            du, u, 
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end

    u.heat[1] += u.temp[1]

    #du_vec .= Vector(du.temp)

    #for reg in region_groups
    #    debug_region!(du, u, reg)
    #end
end