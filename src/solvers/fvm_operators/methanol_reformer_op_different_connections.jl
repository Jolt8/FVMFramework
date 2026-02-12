function solve_connection_group!(
        du, u, phys_a::TA, phys_b::TB, flux!, neighbors,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
    ) where {TA, TB}
    #look into storing this data into a CSR in the future
    @batch for (idx_a, neighbor_list) in neighbors
        for (idx_b, face_idx) in neighbor_list
            flux!(
                du, u, phys_a, phys_b,
                idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
            )
        end
    end
end

function solve_controller_group!(
        du, u, controller::T, controller_id, control!, monitored_cells, affected_cells, 
        cell_volumes
    ) where T
    #no loop here because control! already loops over monitored_cells and affected_cells
    control!(
        du, u, controller, controller_id, monitored_cells, affected_cells,
        cell_volumes
    )
end

function solve_region_group!(
        du, u, phys::T, internal_physics!, region_cells
    ) where T
    @batch for cell_id in region_cells
        internal_physics!(
            du, u, phys, cell_id,
            cell_volumes[cell_id]
        )
    end
end

function FVM_Tracer_Operator!(
    du_trace, u_trace,
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,
    connection_groups, controller_groups, region_groups
)

    #Connection Loops
    #only neighbors and flux_function! is needed now
    for conn in connection_groups
        @batch for (idx_a, neighbor_list) in conn.cell_neighbors
            for (idx_b, face_idx) in neighbor_list
                conn.flux_function!(
                    du_trace, u_trace, idx_a, idx_b, face_idx,
                    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
                )
            end
        end
    end
    
    #Controller Loops
    #only controller_funciton, monitored_cells, and affected_cells are needed now
    for cont in controller_groups
        cont.controller_function!(
            du_trace, u_trace, cont.controller, cont.controller_id, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
    end

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    #oh wait, now we don't even need the other fields for the different regions, we only need the region function
    for reg in region_groups
        @batch for cell_id in reg.region_cells
            reg.internal_physics!(
                du_trace, u_trace, cell_id,
                cell_volumes[cell_id]
            )
        end
    end
end

#oh, since there's no more physics types, the methods above with where T are no longer needed, nice!

function methanol_reformer_f_test!(
    du, u, p, t,
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,
    connection_groups, controller_groups, region_groups,
    #phys::Vector{AbstractPhysics}, cell_phys_id_map::Vector{Int}, #we probably won't use these again, but they might be useful in the future
    du_caches, caches,
    u_merged_axes
)

    u = ComponentVector(u, u_merged_axes)
    du = ComponentVector(du, u_merged_axes)

    caches = get_tmp(caches, u)
    du_caches = get_tmp(du_caches, u)

    u = ComponentVector(u; NamedTuple(caches)...)
    du = ComponentVector(du; NamedTuple(du_caches)...)
    du .= 0.0

    #u = [u; caches] this is another option, but it might allocate

    #even though react_cell!() already sets change_in_molar_concentrations_cache to .= 0.0, we're just doing .= 0.0 to all to make sure
    #NOTE: I'm pretty sure that .*= 0.0 is ever so slightly less performant than .= 0.0

    #I think one of the hardest things we're going to have to do is to create methods to construct the merged u and merged du from input data

    #Connection Loops
    #only neighbors and flux_function! is needed now
    for conn in connection_groups
        @batch for (idx_a, neighbor_list) in conn.cell_neighbors
            for (idx_b, face_idx) in neighbor_list
                conn.flux_function!(
                    du, u, idx_a, idx_b, face_idx,
                    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
                )
            end
        end
    end
    
    #Controller Loops
    #only controller_funciton, monitored_cells, and affected_cells are needed now
    for cont in controller_groups
        cont.controller_function!(
            du, u, cont.controller_id, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
    end

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    #oh wait, now we don't even need the other fields for the different regions, we only need the region function
    for reg in region_groups
        @batch for cell_id in region_cells
            reg.internal_physics!(
                du, u, cell_id,
                cell_volumes[cell_id]
            )
        end
    end
end

#VERY IMPORTANT!!!!
#= For future reference when getting properties using u[field]:
    u_cv.temp[cell_id] = val
        works!
    u_cv[field][cell_id] = val
        fails :(
    view(u_cv, field)[cell_id] = val
        works!
    getproperty(u_cv, field)[cell_id] = val
        works!
=#

#instead of updating rho in each flux and internal physics, we could do an initial property retrieval loop
#the only disadvantage is that once again, solids can't do rho_ideal!() and just have phys.rho
#=
@batch for cell_id in eachindex(cell_volumes)
    mw_avg!(u, cell_id, molecular_weights, mw_avg_cache)
    rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)
end
=#