function generate_ray_map_between_nodes(transducer_ids, node_ids, grid, geo)
    ray_map_intersected_cells = [Int[] for i in 1:length(transducer_ids), j in 1:length(transducer_ids)]
    ray_map_distances_through_cells = [Float64[] for i in 1:length(transducer_ids), j in 1:length(transducer_ids)]
    ray_map_ray_lengths = zeros(Float64, length(transducer_ids), length(transducer_ids))

    for origin_transducer_id in transducer_ids
        origin_node_id = transducer_node_ids[origin_transducer_id]
        for destination_transducer_id in transducer_ids
            destination_node_id = transducer_node_ids[destination_transducer_id]

            if origin_transducer_id != destination_transducer_id #I'm lazy, so we're not slicing the array
                intersected_cells, corresponding_intersection_lengths, ray_length = find_traced_cells(origin_node_id, destination_node_id, grid)

                ray_map_intersected_cells[origin_transducer_id, destination_transducer_id] = intersected_cells
                ray_map_distances_through_cells[origin_transducer_id, destination_transducer_id] = corresponding_intersection_lengths
                ray_map_ray_lengths[origin_transducer_id, destination_transducer_id] = ray_length
            end
        end
    end
    return ray_map_intersected_cells, ray_map_distances_through_cells, ray_map_ray_lengths
end


function generate_ray_map(transducer_ids, node_ids, grid, geo)
    ray_map_intersected_cells, ray_map_distances_through_cells, ray_map_ray_lengths = generate_ray_map_between_nodes(transducer_ids, node_ids, grid, geo)
    
    save_object(joinpath(
            @__DIR__, 
            "saved_ray_maps", 
            "ray_map_$(length(transducer_ids))_transducers_$(length(grid.cells))_cells.jls"
        ), 
        (ray_map_intersected_cells, ray_map_distances_through_cells, ray_map_ray_lengths)
    )
end