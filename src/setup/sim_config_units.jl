
function run_and_check_units(du0_vec_units, u0_vec_units, geo, system, du_unitful_cache_vec, u_unitful_cache_vec)
    du_unitful_cache_vec = upreferred.(du_unitful_cache_vec)
    du_unitful_cache_vec .*= 0.0
    du_unitful_cache_vec ./= 1.0u"s"
    u_unitful_cache_vec = upreferred.(u_unitful_cache_vec)
    u_unitful_cache_vec .*= 0.0

    du_cache_nt = create_views_inline(upreferred.(du_unitful_cache_vec), system.du_cache_axes)
    u_cache_nt = create_views_inline(upreferred.(u_unitful_cache_vec), system.u_cache_axes)

    du0_vec_units .*= 0.0
    du0_vec_units ./= 1.0u"s"

    du_nt = create_views_inline(upreferred.(du0_vec_units), system.du_proto_axes)
    u_nt = create_views_inline(upreferred.(u0_vec_units), system.u_proto_axes)

    du = (; du_nt..., du_cache_nt...)
    u = (; system.merged_properties..., u_nt..., u_cache_nt...)
    #remember, the right-most fields overwrite other fields with the same name

    u.rho .= (u_cache_nt.rho .+ system.merged_properties.rho)


    #applying units 
    cell_volumes = geo.cell_volumes .* u"m^3"
    cell_centroids = geo.cell_centroids .* u"m"
    
    cell_neighbor_areas = geo.cell_neighbor_areas .* u"m^2"
    cell_neighbor_normals = geo.cell_neighbor_normals .* u"m"
    cell_neighbor_distances = geo.cell_neighbor_distances .* u"m"

    cell_face_areas = geo.cell_face_areas .* u"m^2"
    cell_face_normals = geo.cell_face_normals .* u"m"

    for conn in system.connection_groups
        solve_connection_group!(
            du, u,
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        )
    end

    for patch in system.patch_groups
        solve_patch_group!(
            du, u,
            patch.patch_function!, patch.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            cell_volumes
        )
    end
    
    for reg in system.region_groups
        solve_region_group!(
            du, u,
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end

    return du, u
end

#just an idea, but we could also run a separate run-through of the solver that doesn't contain units 
#and then at the end strip the units from the unitful run-through and compare the two
#while I'm pretty sure that upreferred would take care of most of this, it would make you extra
#confident that you're not missing any unit conversions

