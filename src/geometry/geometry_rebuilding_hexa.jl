function calculate_hex_volume(p)
    c = sum(p) / 8.0
    
    # Faces of Hex8 (node indices for each face)
    faces = (
        (1, 4, 3, 2),  # bottom
        (1, 2, 6, 5),  # front
        (2, 3, 7, 6),  # right
        (3, 4, 8, 7),  # back
        (1, 5, 8, 4),  # left 
        (5, 6, 7, 8)   # top
    )
    
    total_vol = 0.0
    for face in faces
        node_1, node_2, node_3, node_4 = p[face[1]], p[face[2]], p[face[3]], p[face[4]]
        
        # Split quad face into 2 triangles, form tetrahedra with centroid
        total_vol += dot(node_1 - c, cross_product(node_2 - c, node_3 - c))
        total_vol += dot(node_1 - c, cross_product(node_3 - c, node_4 - c))
    end
    
    return abs(total_vol) / 6.0
end

function calculate_quad_face_area(face_node_coordiantes)
    node_1_coords = face_node_coordiantes[1]
    node_2_coords = face_node_coordiantes[2]
    node_3_coords = face_node_coordiantes[3]
    node_4_coords = face_node_coordiantes[4]

    cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
    cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
    
    area_vec_1 = 0.5 * cross_a
    area_vec_2 = 0.5 * cross_b

    total_area_vec = area_vec_1 + area_vec_2
    return norm(total_area_vec)
end

function rebuild_fvm_geometry(
        cell_neighbors, cell_neighbors_node_ids::Vector{MVector{6, SVector{4, Int}}}, 
        all_cell_face_map, map_respective_node_ids::Vector{NTuple{4, Int}}, 
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

        p = ntuple(8) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        vol = calculate_hex_volume(p)
        
        cent = sum(p) / T(n_nodes)

        cell_volumes[cell_id] = vol
        cell_centroids[cell_id] = cent
    end

    n_cell_neighbors = length(cell_neighbors)
    
    connection_areas = Vector{T}(undef, n_cell_neighbors)
    connection_normals = Vector{CoordType}(undef, n_cell_neighbors)
    #the above connection_normals causes some GC (4% of optimization runtime)
    connection_distances = Vector{T}(undef, n_cell_neighbors)
    

    for this_cell_neighbors in cell_neighbors
        for (face_idx, neighbor_id) in enumerate(this_cell_neighbors)
            face_node_indices = cell_neighbors_node_ids[face_idx] #cell_neighbors_node_ids[face_idx] looks like (1, 4, 7, 21) 
            node_1_coords = node_coordinates[face_node_indices[1]]
            node_2_coords = node_coordinates[face_node_indices[2]]
            node_3_coords = node_coordinates[face_node_indices[3]]
            node_4_coords = node_coordinates[face_node_indices[4]]

            #get_area
            cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
            cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
            
            area_vec_1 = 0.5 * cross_a
            area_vec_2 = 0.5 * cross_b

            total_area_vec = area_vec_1 + area_vec_2
            total_area = norm(total_area_vec)

            #make sure it's actually pointing away
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
    end
    
    cell_face_areas = fill(zero(MVector{6, T}), n_cells)
    cell_face_normals = fill(zero(MVector{6, CoordType}), n_cells)
    
    #if the above breaks:
    #unconnected_areas = Vector{SVector{6, T}}(undef, n_cells)
    #unconnected_normals = Vector{SVector{6, CoordType}}(undef, n_cells)
    
    
    for (i, (cell_id, face_idx)) in enumerate(all_cell_face_map)
        #=
        areas = MVector{6, T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        normals = MVector{6, CoordType}(
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0)
        )
        =#
        
        face_node_indices = map_respective_node_ids[i] #unconnected_map_respective_node_ids[i] looks like (1, 4, 7, 21) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]
        node_4_coords = node_coordinates[face_node_indices[4]]

        #get_area
        cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
        cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
        
        area_vec_1 = 0.5 * cross_a
        area_vec_2 = 0.5 * cross_b

        total_area_vec = area_vec_1 + area_vec_2
        total_area = norm(total_area_vec)

        cell_face_areas[cell_id][face_idx] = total_area 

        vec_out = (node_1_coords + node_2_coords + node_3_coords + node_4_coords) / 4 - cell_centroids[cell_id]

        cell_face_normals[cell_id][face_idx] = normalize(vec_out)

        #cell_face_areas[cell_id] = SVector(areas)
        #cell_face_normals[cell_id] = SVector(normals)
    end

    frozen_areas = Vector{SVector{6, T}}(undef, n_cells)
    frozen_normals = Vector{SVector{6, CoordType}}(undef, n_cells)

    #we freeze them back into SVectors
    for cell_id in eachindex(cell_face_areas)
        frozen_areas[cell_id] = SVector(cell_face_areas[cell_id])
        frozen_normals[cell_id] = SVector(cell_face_normals[cell_id])
    end
    
    return cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances, frozen_areas, frozen_normals
end