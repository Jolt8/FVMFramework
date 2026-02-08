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
    cell_neighbor_map, neighbor_map_respective_node_ids::Vector{NTuple{3, Int}},
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

    n_connections = length(cell_neighbor_map)

    connection_areas = Vector{T}(undef, n_connections)
    connection_normals = Vector{CoordType}(undef, n_connections)
    #the above connection_normals causes some GC (4% of optimization runtime)
    connection_distances = Vector{T}(undef, n_connections)


    for (i, (cell_id, neighbor_id)) in enumerate(cell_neighbor_map)
        face_node_indices = neighbor_map_respective_node_ids[i] #map_respective_node_ids[i] looks like (1, 4, 7, 21) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]

        #get_area
        total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

        total_area = norm(total_area_vec)

        #make sure it's actuall ypointing away
        vec_AB = cell_centroids[neighbor_id] - cell_centroids[cell_id]
        if dot(total_area_vec, vec_AB) < 0
            total_area_vec = -total_area_vec
        end

        #get normal
        cell_normal = normalize(total_area_vec)

        #get distance 
        dist = norm(cell_centroids[cell_id] - cell_centroids[neighbor_id])

        connection_areas[i] = total_area
        connection_normals[i] = cell_normal
        connection_distances[i] = dist
    end

    cell_face_areas = fill(zero(MVector{4, T}), n_cells)
    cell_face_normals = fill(zero(MVector{4, CoordType}), n_cells)

    #= if the above breaks:
    areas = Vector{SVector{4, T}}(undef, n_cells)
    normals = Vector{SVector{4, CoordType}}(undef, n_cells)
    =#

    for (i, (cell_id, face_idx)) in enumerate(all_cell_face_map)
        #=
        curr_areas = MVector{4,T}(0.0, 0.0, 0.0, 0.0)
        curr_normals = MVector{4,CoordType}(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        )
        =#

        face_node_indices = map_respective_node_ids[i] #unconnected_map_respective_node_ids[i] looks like (1, 4, 7) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]

        #get_area
        total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

        total_area = norm(total_area_vec)

        cell_face_areas[cell_id][face_idx] = total_area

        vec_out = (node_1_coords + node_2_coords + node_3_coords) / 3 - cell_centroids[cell_id]

        cell_face_normals[cell_id][face_idx] = normalize(vec_out)

        #cell_face_areas[cell_id] = SVector(areas)
        #cell_face_normals[cell_id] = SVector(normals)
    end

    frozen_areas = Vector{SVector{4, T}}(undef, n_cells)
    frozen_normals = Vector{SVector{4, CoordType}}(undef, n_cells)

    #we freeze them back into SVectors
    for cell_id in eachindex(cell_face_areas)
        frozen_areas[cell_id] = SVector(cell_face_areas[cell_id])
        frozen_normals[cell_id] = SVector(cell_face_normals[cell_id])
    end

    return cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances, frozen_areas, frozen_normals
end