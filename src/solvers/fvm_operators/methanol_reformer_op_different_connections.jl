function solve_connection_group!(
    du, u, phys_a::TA, phys_b::TB, flux!, neighbors,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
) where {TA,TB}
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
    du, u, controller::T, id, control!, monitored_cells, affected_cells,
    cell_volumes
) where T
    #no loop here because control! already loops over monitored_cells and affected_cells
    control!(
        du, u, controller, id, monitored_cells, affected_cells,
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
            du_trace, u_trace, cont.id, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
    end

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    #oh wait, now we don't even need the other fields for the different regions, we only need the region function
    for reg in region_groups
        @batch for cell_id in reg.region_cells
            reg.region_function!(
                du_trace, u_trace, cell_id,
                cell_volumes[cell_id]
            )
        end
    end
end

#oh, since there's no more physics types, the methods above with where T are no longer needed, nice!

function methanol_reformer_f_test!(
    du_flat, u_flat, p, t,
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,
    connection_groups, controller_groups, region_groups,
    du_fixed, u_fixed,
    state_axes
)
    # Wrap incoming flat vectors in ComponentVectors (these are views of the solver arrays)
    du_cv = ComponentVector(du_flat, state_axes)
    u_cv = ComponentVector(u_flat, state_axes)

    # Get the caches (which are usually stored in DiffCaches)
    if eltype(u) == Float64
        du_fixed_cv = get_tmp(du_fixed, Float64)
        u_fixed_cv = get_tmp(u_fixed, Float64)
    else
        du_fixed_cv = get_tmp(du_fixed, ForwardDiff.Dual(1.0))
        u_fixed_cv = get_tmp(u_fixed, ForwardDiff.Dual(1.0))
    end

    # Create merged ComponentVectors for the physics functions to use.
    # Note: ComponentVector(cv; NamedTuple(cv_fixed)...) creates a NEW array.
    du_merged = ComponentVector(du_cv; NamedTuple(du_fixed_cv)...)
    u_merged = ComponentVector(u_cv; NamedTuple(u_fixed_cv)...)

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

    # Controller Loops
    for cont in controller_groups
        cont.controller_function!(
            du_merged, u_merged, cont.id, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
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