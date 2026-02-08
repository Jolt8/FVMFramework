#COMMITING BEFORE SWITCHING TO A CELL CENTRIC CONNECTION MAP

function calculate_tetra_volume(p)
    volume = (1 / 6) * abs(dot(p[1] - p[4], cross_product((p[2] - p[4]), (p[3] - p[4]))))

    return volume
end

function calculate_tri_face_area(face_node_coordiantes)
    node_1_coords = face_node_coordiantes[1]
    node_2_coords = face_node_coordiantes[2]
    node_3_coords = face_node_coordiantes[3]

    total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

    return norm(total_area_vec)
end

function rebuild_fvm_geometry(
    cell_neighbors, cell_neighbors_node_ids::Vector{MVector{4, SVector{3, Int}}},
    all_cell_face_map, map_respective_node_ids::Vector{NTuple{3, Int}},
    node_coordinates, nodes_of_cells
)
    #The cells_data and connections map are still causing GC
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    n_cells = length(nodes_of_cells)

    cell_volumes = Vector{T}(undef, n_cells)
    cell_centroids = Vector{CoordType}(undef, n_cells)

    for cell_id in eachindex(nodes_of_cells)
        cell_nodes = nodes_of_cells[cell_id]

        n_nodes = length(cell_nodes)

        p = ntuple(4) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        vol = calculate_tetra_volume(p)

        cent = sum(p) / T(n_nodes)

        cell_volumes[cell_id] = vol
        cell_centroids[cell_id] = cent
    end

    n_cell_neighbors = length(cell_neighbors)

    connection_areas = Vector{MVector{4, T}}(undef, n_cell_neighbors)
    connection_normals = Vector{MVector{4, CoordType}}(undef, n_cell_neighbors)
    connection_distances = Vector{MVector{4, T}}(undef, n_cell_neighbors)

    #for geometry optimization in the future, we might want to store a separate version of cell_neighbors that doesn't contain duplicates
    #this would only require making normals negative for the neighboring cell and would reduce the math needed
    for this_cell_neighbors in cell_neighbors
        for (face_idx, neighbor_id) in enumerate(this_cell_neighbors)
            face_node_indices = cell_neighbors_node_ids[face_idx] #cell_neighbors_node_ids[face_idx] looks like (1, 4, 7) 
            node_1_coords = node_coordinates[face_node_indices[1]]
            node_2_coords = node_coordinates[face_node_indices[2]]
            node_3_coords = node_coordinates[face_node_indices[3]]

            #get area
            total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

            total_area = norm(total_area_vec)

            #get normal
            #make sure the cell's normal is actually pointing away
            vec_AB = cell_centroids[neighbor_id] - cell_centroids[cell_id]
            if dot(total_area_vec, vec_AB) < 0
                total_area_vec = -total_area_vec
            end

            cell_normal = normalize(total_area_vec)

            #get distance
            dist = norm(cell_centroids[cell_id] - cell_centroids[neighbor_id])

            connection_areas[cell_id][face_idx] = total_area
            connection_normals[cell_id][face_idx] = cell_normal
            connection_distances[cell_id][face_idx] = dist
        end
    end

    cell_face_areas = fill(zero(MVector{4, T}), n_cells)
    cell_face_normals = fill(zero(MVector{4, CoordType}), n_cells)

    for (i, (cell_id, face_idx)) in enumerate(all_cell_face_map)
        face_node_indices = map_respective_node_ids[i] #unconnected_map_respective_node_ids[i] looks like (1, 4, 7) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]

        total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

        total_area = norm(total_area_vec)

        cell_face_areas[cell_id][face_idx] = total_area

        vec_out = (node_1_coords + node_2_coords + node_3_coords) / 3 - cell_centroids[cell_id]

        cell_face_normals[cell_id][face_idx] = normalize(vec_out)
    end

    return cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances, cell_face_areas, cell_face_normals
end