#=
#this is an example of what else we could do to construct du_merged, but I've found that it doesn't actually perform better


N::Int = ForwardDiff.pickchunksize(length(complete_state))
cache_cv = build_component_array_merge_regions(cache_vars, region_symbols, n_cells, n_faces)

complete_fixed = ComponentVector(complete_fixed; NamedTuple(cache_cv)...)

du_fixed = copy(complete_fixed)
du_fixed .= 0.0
u_fixed = copy(complete_fixed)

fixed_axes = getaxes(complete_fixed)[1]

u_fixed_vec = DiffCache(Vector(complete_fixed), N)
du_fixed_vec = DiffCache(Vector(du_fixed), N)

u0_flat = Vector(complete_state)
du0_flat = copy(u0_flat)
du0_flat .= 0.0

state_axes = getaxes(complete_state)[1]

all_vars = ComponentVector(complete_state; NamedTuple(complete_fixed)...)
u_all_vars_vec = Vector(all_vars)
du_all_vars_vec = copy(u_all_vars_vec)
all_axes = getaxes(all_vars)[1]

f_closure = (du, u, p, t) -> heat_transfer_f!(
    du, u, p, t,
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,
    system.connection_groups, system.controller_groups, system.region_groups,
    Vector(du_fixed), Vector(u_fixed),
    state_axes, fixed_axes,
    #u_all_vars_vec, du_all_vars_vec, all_axes
)
    
function heat_transfer_f!(
    du_flat, u_flat, p, t,
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,
    connection_groups, controller_groups, region_groups,
    du_fixed_vec, u_fixed_vec,
    state_axes, fixed_axes,
    du_all_vars_vec, u_all_vars_vec, all_vars_axes
)
    du_cv = ComponentVector(du_flat, state_axes)
    u_cv = ComponentVector(u_flat, state_axes)

    #println(u_cv)

    if eltype(du_flat) == Float64
        du_fixed_vec = get_tmp(du_fixed_vec, Float64)
        u_fixed_vec = get_tmp(u_fixed_vec, Float64)
    else
        du_fixed_vec = get_tmp(du_fixed_vec, ForwardDiff.Dual(1.0))
        u_fixed_vec = get_tmp(u_fixed_vec, ForwardDiff.Dual(1.0))
    end

    #du_fixed_cv = ComponentVector(du_fixed_vec, fixed_axes)
    #u_fixed_cv = ComponentVector(u_fixed_vec, fixed_axes)

    #du_merged = ComponentVector(du_cv; NamedTuple(du_fixed_cv)...)
    #u_merged = ComponentVector(u_cv; NamedTuple(u_fixed_cv)...)

    du_all_vars_vec[1:length(du_flat)] .= du_flat
    du_all_vars_vec[length(du_flat)+1:length(du_flat)+length(du_fixed_vec)] .= du_fixed_vec
    u_all_vars_vec[1:length(u_flat)] .= u_flat
    u_all_vars_vec[length(u_flat)+1:length(u_flat)+length(u_fixed_vec)] .= u_fixed_vec

    #println(u_all_vars_vec)

    du_merged = ComponentVector(du_all_vars_vec, all_vars_axes)
    u_merged = ComponentVector(u_all_vars_vec, all_vars_axes)

    #println(u_merged)

    du_merged .= 0.0

    #println("u_merged[1], ", u_merged[1])

    # Connection Loops
    for conn in connection_groups
        for (idx_a, neighbor_list) in conn.cell_neighbors
            for (idx_b, face_idx) in neighbor_list
                conn.flux_function!(
                    du_merged, u_merged, idx_a, idx_b, face_idx,
                    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
                )
            end
        end
    end

    #println(du_merged)

    # Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    for reg in region_groups
        for cell_id in reg.region_cells
            reg.region_function!(
                du_merged, u_merged, cell_id,
                cell_volumes
            )
        end
    end
    #println(u_merged)
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
=#