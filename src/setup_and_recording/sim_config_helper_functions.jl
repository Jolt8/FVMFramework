function get_cell_ids_in_facet_set(grid, set_name)
    cell_ids_in_facet_set = Int[]

    for (cell_id, facet_idx) in getfacetset(grid, set_name)
        push!(cell_ids_in_facet_set, cell_id)
    end

    return cell_ids_in_facet_set
end