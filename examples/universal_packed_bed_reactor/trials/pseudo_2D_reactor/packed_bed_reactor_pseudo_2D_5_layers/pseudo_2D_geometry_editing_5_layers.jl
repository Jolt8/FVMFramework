#Script for artificially enlarging the area and volume of each region to match the physical reactor
#this was done because a full 3D reactor would be too computationally expensive and because visualizing/working with 
#a 2D model is tedious and annoying
#so we use 1D + volume and area corrections 

function edit_pseudo_2D_geometry_5_layers!(geo, grid, common_properties)


region_cell_areas = Dict(
    ("pipe_inlet", "steel_pipe_wall") => common_properties.reforming_area_to_steel_pipe_wall_cell_areas,
    ("silicon_carbide_preheater", "steel_pipe_wall") => common_properties.reforming_area_to_steel_pipe_wall_cell_areas,
    ("copper_mesh_reformer", "steel_pipe_wall") => common_properties.reforming_area_to_steel_pipe_wall_cell_areas,
    ("pipe_outlet", "steel_pipe_wall") => common_properties.reforming_area_to_steel_pipe_wall_cell_areas,
    ("steel_pipe_wall", "thermocouple_and_jacket") => common_properties.steel_pipe_wall_to_thermocouple_cell_areas,
    ("thermocouple_and_jacket", "heating_wire") => common_properties.thermocouple_region_to_heating_wire_cell_areas,
    ("heating_wire", "insulation") => common_properties.heating_wire_to_insulation_cell_areas,
    #("insulation", "environment") => insulation_to_environment_cell_areas #the parameter that we're using to fit heat transfer to the environment is lumped 
    #in with the area of the insulation that's exposed to the environment
)

region_cell_distances = Dict(
    ("pipe_inlet", "steel_pipe_wall") => ((common_properties.steel_pipe_wall_outside_diameter / 2) + (common_properties.steel_pipe_wall_inside_diameter / 2)) / 2 - 
    ((common_properties.reforming_area_outside_diameter / 2) + 0.0u"m") / 2,
    ("silicon_carbide_preheater", "steel_pipe_wall") => ((common_properties.steel_pipe_wall_outside_diameter / 2) + (common_properties.steel_pipe_wall_inside_diameter / 2)) / 2 - 
    ((common_properties.reforming_area_outside_diameter / 2) + 0.0u"m") / 2,
    ("copper_mesh_reformer", "steel_pipe_wall") => ((common_properties.steel_pipe_wall_outside_diameter / 2) + (common_properties.steel_pipe_wall_inside_diameter / 2)) / 2 - 
    ((common_properties.reforming_area_outside_diameter / 2) + 0.0u"m") / 2,
    ("pipe_outlet", "steel_pipe_wall") => ((common_properties.steel_pipe_wall_outside_diameter / 2) + (common_properties.steel_pipe_wall_inside_diameter / 2)) / 2 - 
    ((common_properties.reforming_area_outside_diameter / 2) + 0.0u"m") / 2,
    ("steel_pipe_wall", "thermocouple_and_jacket") => ((common_properties.thermocouple_region_outside_diameter / 2) + (common_properties.thermocouple_region_inside_diameter / 2)) / 2 - 
    ((common_properties.steel_pipe_wall_outside_diameter / 2) + (common_properties.steel_pipe_wall_inside_diameter / 2)) / 2,
    ("thermocouple_and_jacket", "heating_wire") => ((common_properties.heating_wire_outside_diameter / 2) + (common_properties.heating_wire_inside_diameter / 2)) / 2 - 
    ((common_properties.thermocouple_region_outside_diameter / 2) + (common_properties.thermocouple_region_inside_diameter / 2)) / 2,
    ("heating_wire", "insulation") => ((common_properties.insulation_outside_diameter / 2) + (common_properties.insulation_inside_diameter / 2)) / 2 - 
    ((common_properties.heating_wire_outside_diameter / 2) + (common_properties.heating_wire_inside_diameter / 2)) / 2
    #("insulation", "environment") => insulation_to_environment_cell_distances #the parameter that we're using to fit heat transfer to the environment is lumped 
    #in with the area of the insulation that's exposed to the environment
)

#the way this works is that we find the "centroid"/average radius of a region like the insulation by taking the (outside radius + inside radius) / 2
#then we find the same for the region adjacent to it (like the heating wire) by taking the (outside radius + inside radius) / 2
#we then take the difference between these two values to get teh distance between the average radii

for (region_a, region_b) in keys(region_cell_areas)
    region_a_cells = collect(grid.cellsets[region_a])
    region_b_cells = collect(grid.cellsets[region_b])
    
    for (idx_a, neighbor_list) in geo.cell_neighbors[region_a_cells]
        for (idx_b, face_idx) in neighbor_list
            if idx_b in region_b_cells
                geo.cell_neighbor_areas[idx_a][face_idx] = ustrip(upreferred(region_cell_areas[(region_a, region_b)]))
                geo.cell_face_areas[idx_a][face_idx] = ustrip(upreferred(region_cell_areas[(region_a, region_b)]))

                geo.cell_neighbor_distances[idx_a][face_idx] = ustrip(upreferred(region_cell_distances[(region_a, region_b)]))
            end
        end
    end

    for (idx_a, neighbor_list) in geo.cell_neighbors[region_b_cells]
        for (idx_b, face_idx) in neighbor_list
            if idx_b in region_a_cells
                geo.cell_neighbor_areas[idx_a][face_idx] = ustrip(upreferred(region_cell_areas[(region_a, region_b)]))
                geo.cell_face_areas[idx_a][face_idx] = ustrip(upreferred(region_cell_areas[(region_a, region_b)]))

                geo.cell_neighbor_distances[idx_a][face_idx] = ustrip(upreferred(region_cell_distances[(region_a, region_b)]))
            end
        end
    end
end

for (cell_id, face_idx) in getfacetset(grid, "insulation_to_air")
    geo.cell_face_areas[cell_id][face_idx] = ustrip(upreferred(common_properties.insulation_to_environment_cell_areas))
end

region_cell_volumes = Dict(
    "pipe_inlet" => common_properties.reforming_area_cell_volumes,
    "silicon_carbide_preheater" => common_properties.reforming_area_cell_volumes,
    "copper_mesh_reformer" => common_properties.reforming_area_cell_volumes,
    "pipe_outlet" => common_properties.reforming_area_cell_volumes,
    "steel_pipe_wall" => common_properties.steel_pipe_wall_cell_volumes,
    "thermocouple_and_jacket" => common_properties.thermocouple_region_cell_volumes,
    "heating_wire" => common_properties.heating_wire_cell_volumes,
    "insulation" => common_properties.insulation_cell_volumes
)

for (region_name, corrected_cell_volume) in region_cell_volumes
    region_cells = collect(grid.cellsets[region_name])
    for cell_id in region_cells
        geo.cell_volumes[cell_id] = ustrip(upreferred(corrected_cell_volume))
    end
end

end