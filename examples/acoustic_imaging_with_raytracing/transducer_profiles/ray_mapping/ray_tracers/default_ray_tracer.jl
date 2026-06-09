using Ferrite
using LinearAlgebra
using FVMFramework
using SparseArrays

function find_traced_cells(origin_node_id, destination_node_id, grid, geo)
    origin_coordinates = grid.nodes[origin_node_id].x
    destination_coordinates = grid.nodes[destination_node_id].x
    
    ray_vector = destination_coordinates .- origin_coordinates
    ray_len = norm(ray_vector)
    
    if ray_len < 1e-12
        return Int[], Float64[], 0.0
    end
    
    ray_dir = ray_vector / ray_len
    
    intersected_cells = Tuple{Int, Float64, Float64}[]
    
    for cell_id in 1:getncells(grid)
        n_faces = nfacets(grid.cells[cell_id])
        t_entry = 0.0
        t_exit = ray_len
        intersect = true
        
        for face_idx in 1:n_faces
            # Outward normal of the face
            N_f = geo.cell_face_normals[cell_id][face_idx]
            
            # Get a point on the face by taking the first node of the face
            face_nodes = get_face_nodes(grid, cell_id, face_idx)
            C_f = grid.nodes[face_nodes[1]].x
            
            # Compute dot products
            D = dot(ray_dir, N_f)
            N = dot(C_f - origin_coordinates, N_f)
            
            if abs(D) < 1e-12
                # Parallel to the face plane
                if N < -1e-12
                    intersect = false
                    break
                end
            elseif D > 0.0
                # Pointing outwards relative to this face (exiting half-space)
                t_exit = min(t_exit, N / D)
            else
                # Pointing inwards relative to this face (entering half-space)
                t_entry = max(t_entry, N / D)
            end
        end
        
        if intersect && (t_exit - t_entry > 1e-10)
            push!(intersected_cells, (cell_id, t_exit - t_entry, t_entry)) #cell_id, intersection_distance, entry_distance_from_origin
        end
    end
    
    # Sort the intersected cells in the order the ray encounters them
    sort!(intersected_cells, by = x -> x[3])
    
    # Return a vector of (cell_id, intersection_distance)
    cell_ids = Int[c[1] for c in intersected_cells]
    corresponding_intersection_lengths = Float64[c[2] for c in intersected_cells]
    ray_length = sum(corresponding_intersection_lengths)
    return cell_ids, corresponding_intersection_lengths, ray_length
end