function get_cell_set_total_volume(grid, set_name, geo)
    total_volume = 0.0
    for cell_id in eachindex(geo.cell_volumes)
        if cell_id in getcellset(grid, set_name)
            total_volume += geo.cell_volumes[cell_id]
        end
    end

    return total_volume 
end

function get_facet_set_total_area(grid, set_name, geo)
    #although there's 6 face areas per cell, all we need is the total area of all the facets of the cell that belongs in the facet_set
    total_area = 0.0

    for (cell_id, facet_idx) in getfacetset(grid, set_name)
        face_area = geo.cell_face_areas[cell_id][facet_idx]
        total_area += face_area
    end

    return total_area
end

#this style of function is not needed for cell volumes because it's very easy to access cell_volumes with just [cell_id]
function get_facet_set_cells_respective_areas(grid, set_name, geo)
    n_cells = length(grid.cells)
    cell_ids_respective_areas = zeros(Float64, n_cells) 

    for (cell_id, facet_idx) in getfacetset(grid, set_name)
        face_area = geo.cell_face_areas[cell_id][facet_idx]
        cell_ids_respective_areas[cell_id] += face_area #we += because a cell might have 2 faces in this facet_set
    end

    return cell_ids_respective_areas
end

function get_cell_ids_in_facet_set(grid, set_name)
    cell_ids_in_facet_set = Int[]

    for (cell_id, facet_idx) in getfacetset(grid, set_name)
        push!(cell_ids_in_facet_set, cell_id)
    end

    return cell_ids_in_facet_set
end