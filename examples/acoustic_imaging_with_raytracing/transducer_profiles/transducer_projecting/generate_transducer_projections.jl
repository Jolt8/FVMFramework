
function count_cell_projections(transducer_projection_intersected_cells, n_cells)
    projection_counts = zeros(Int, n_cells)
    for cell_list in transducer_projection_intersected_cells
        for cell_id in cell_list
            projection_counts[cell_id] += 1
        end
    end
    return projection_counts
end

function generate_transducer_projections_between_transducers(
    transducer_opposing_pairs_ids::Vector{Tuple{Int, Int}}, transducer_node_ids::Vector{Int}, grid, geo,
    most_consistent_timestamp_for_speed_of_sound,
    experimental_travel_time_interp_matrix, #used to find the speed of sound
    transducer_frequency,
    transducer_diameter
)
    n_transducers = length(transducer_node_ids)
    transducer_projection_intersected_cells = [Int[] for _ in 1:n_transducers]
    transducer_projection_volume_in_cells = [Float64[] for _ in 1:n_transducers]
    transducer_projection_cell_distances = [Float64[] for _ in 1:n_transducers]
    transducer_projection_distances = zeros(Float64, n_transducers)
    transducer_projection_slant_distances = zeros(Float64, n_transducers)
    projection_unit_vectors_to_cell_centers = [Float64[] for _ in 1:n_transducers]
    
    for (transducer_A_id, transducer_B_id) in transducer_opposing_pairs_ids #A and B capitalized for clarity
        transducer_A_node_id = transducer_node_ids[transducer_A_id]
        transducer_B_node_id = transducer_node_ids[transducer_B_id]

        transducer_projection_intersected_cells[transducer_A_id], 
        transducer_projection_volume_in_cells[transducer_A_id], 
        transducer_projection_cell_distances[transducer_A_id], 
        transducer_projection_distances[transducer_A_id], 
        transducer_projection_slant_distances[transducer_A_id],
        projection_unit_vectors_to_cell_centers[transducer_A_id] = get_cells_in_transducer_projection(
            transducer_A_node_id, #Transducer A is source
            transducer_B_node_id, #Transducer B opposes
            grid, geo,
            most_consistent_timestamp_for_speed_of_sound,
            experimental_travel_time_interp_matrix,
            transducer_frequency,
            transducer_diameter
        )

        transducer_projection_intersected_cells[transducer_B_id], 
        transducer_projection_volume_in_cells[transducer_B_id], 
        transducer_projection_cell_distances[transducer_B_id],
        transducer_projection_distances[transducer_B_id],  
        transducer_projection_slant_distances[transducer_B_id],
        projection_unit_vectors_to_cell_centers[transducer_B_id] = get_cells_in_transducer_projection(
            transducer_B_node_id, #Transducer B is source
            transducer_A_node_id, #Transducer A opposes
            grid, geo,
            most_consistent_timestamp_for_speed_of_sound,
            experimental_travel_time_interp_matrix,
            transducer_frequency,
            transducer_diameter
        )
    end

    cell_projection_counts = count_cell_projections(transducer_projection_intersected_cells, length(grid.cells))

    return transducer_projection_intersected_cells, transducer_projection_volume_in_cells, transducer_projection_distances, transducer_projection_slant_distances, cell_projection_counts, projection_unit_vectors_to_cell_centers
end

function generate_transducer_projections(
    transducer_opposing_pairs_ids::Vector{Tuple{Int, Int}}, transducer_node_ids::Vector{Int}, grid, geo,
    most_consistent_timestamp_for_speed_of_sound,
    experimental_travel_time_interp_matrix, #used to find the speed of sound
    transducer_frequency,
    transducer_diameter
)
    transducer_projection_intersected_cells, 
    transducer_projection_volume_in_cells, 
    transducer_projection_cell_distances, 
    transducer_projection_distances, 
    transducer_projection_slant_distances,
    cell_projection_counts,
    projection_unit_vectors_to_cell_centers = generate_transducer_projections_between_transducers(
        transducer_opposing_pairs_ids, 
        transducer_node_ids, 
        grid, geo,
        most_consistent_timestamp_for_speed_of_sound,
        experimental_travel_time_interp_matrix, #used to find the speed of sound
        transducer_frequency,
        transducer_diameter
    )
    
    save_object(joinpath(
            @__DIR__, 
            "saved_transducer_projections", 
            "transducer_projection_$(length(transducer_opposing_pairs_ids))_pairs_$(length(grid.cells))_cells.jls"
        ), 
        (
            transducer_projection_intersected_cells, 
            transducer_projection_volume_in_cells, 
            transducer_projection_cell_distances, 
            transducer_projection_distances, 
            transducer_projection_slant_distances, 
            cell_projection_counts,
            projection_unit_vectors_to_cell_centers
        )
    )
end
