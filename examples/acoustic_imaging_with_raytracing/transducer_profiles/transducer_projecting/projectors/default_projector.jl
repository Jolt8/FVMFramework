using Ferrite
using LinearAlgebra
using FVMFramework

# Helper function to check if a point is inside the projection cone
function is_point_in_cone(P, source_transducer_coordinates, unit_ray_vector, projection_distance, cos_theta)
    v = P - source_transducer_coordinates
    d_axis = dot(v, unit_ray_vector)
    if d_axis < 0.0 || d_axis > projection_distance
        return false
    end
    dist_from_apex = norm(v)
    if dist_from_apex < 1e-12
        return true
    end
    return (d_axis / dist_from_apex) >= cos_theta
end


function get_cells_in_transducer_projection(
    source_transducer_node_id, opposing_transducer_node_id, grid, geo,
    most_consistent_timestamp_for_speed_of_sound,
    experimental_travel_time_interp_matrix, #used to find the speed of sound
    transducer_frequency,
    transducer_diameter
)
    source_transducer_coordinates = grid.nodes[source_transducer_node_id].x
    opposing_transducer_coordinates_raw = grid.nodes[opposing_transducer_node_id].x

    # If the opposing transducer is lower than the source transducer, we don't want to take that into account.
    # Instead, we just want where the source transducer is facing along the x and y axis.
    # Note: Ferrite Vecs are immutable, so we construct a new Vec rather than mutating.
    opposing_transducer_coordinates = Ferrite.Vec{3}((
        opposing_transducer_coordinates_raw[1],
        opposing_transducer_coordinates_raw[2],
        source_transducer_coordinates[3]
    ))
    
    ray_vector = opposing_transducer_coordinates .- source_transducer_coordinates
    unit_ray_vector = normalize(ray_vector)

    projection_distance = sqrt(
        (opposing_transducer_coordinates[1] - source_transducer_coordinates[1])^2 + 
        (opposing_transducer_coordinates[2] - source_transducer_coordinates[2])^2 + 
        (opposing_transducer_coordinates[3] - source_transducer_coordinates[3])^2
    )

    experimental_travel_time_at_consistent_timestamp = experimental_travel_time_interp_matrix[opposing_transducer_node_id, source_transducer_node_id](most_consistent_timestamp_for_speed_of_sound)
    speed_of_sound = projection_distance / experimental_travel_time_at_consistent_timestamp

    wavelength = speed_of_sound / transducer_frequency
    arg = 1.22 * wavelength / transducer_diameter
    if arg >= 1.0
        angle_of_projection = 90.0
    else
        angle_of_projection = asind(arg)
    end
    cos_theta = cosd(angle_of_projection)

    slant_distance = projection_distance / cos_theta

    projection_intersected_cells = Int[]
    projection_volume_in_cells = Float64[]
    projection_distances_from_source = Float64[]
    projection_unit_vectors_to_cell_centers = Vec{3, Float64}[]

    CellType = typeof(grid.cells[1])

    for cell_id in 1:length(grid.cells)
        cell_centroid = geo.cell_centroids[cell_id]
        cell_volume = geo.cell_volumes[cell_id]
        cell_nodes = grid.cells[cell_id].nodes
        n_nodes = length(cell_nodes)

        # Check all vertices and the centroid
        vertices_inside = 0
        p = ntuple(i -> grid.nodes[cell_nodes[i]].x, n_nodes)
        for i in 1:n_nodes
            if is_point_in_cone(p[i], source_transducer_coordinates, unit_ray_vector, projection_distance, cos_theta)
                vertices_inside += 1
            end
        end

        centroid_inside = is_point_in_cone(cell_centroid, source_transducer_coordinates, unit_ray_vector, projection_distance, cos_theta)

        if vertices_inside == n_nodes
            # Entirely inside the cone
            push!(projection_intersected_cells, cell_id)
            push!(projection_volume_in_cells, cell_volume)
            push!(projection_distances_from_source, norm(cell_centroid - source_transducer_coordinates))
        elseif vertices_inside > 0 || centroid_inside
            # Partially inside the cone, estimate volume via sub-cell sampling
            n_inside = 0
            n_total = 0
            if CellType <: Hexahedron
                n_total = 27
                for r in (0.25, 0.5, 0.75), s in (0.25, 0.5, 0.75), t in (0.25, 0.5, 0.75)
                    pt = (1 - r) * (1 - s) * (1 - t) * p[1] +
                         r * (1 - s) * (1 - t) * p[2] +
                         r * s * (1 - t) * p[3] +
                         (1 - r) * s * (1 - t) * p[4] +
                         (1 - r) * (1 - s) * t * p[5] +
                         r * (1 - s) * t * p[6] +
                         r * s * t * p[7] +
                         (1 - r) * s * t * p[8]
                    if is_point_in_cone(pt, source_transducer_coordinates, unit_ray_vector, projection_distance, cos_theta)
                        n_inside += 1
                    end
                end
            elseif CellType <: Tetrahedron
                for r in 0.1:0.2:0.9, s in 0.1:0.2:0.9, t in 0.1:0.2:0.9
                    if r + s + t <= 0.95
                        n_total += 1
                        pt = (1 - r - s - t) * p[1] + r * p[2] + s * p[3] + t * p[4]
                        if is_point_in_cone(pt, source_transducer_coordinates, unit_ray_vector, projection_distance, cos_theta)
                            n_inside += 1
                        end
                    end
                end
            else
                @error "Cell type $(CellType) not supported"
            end

            if n_total > 0 
                f = n_inside / n_total
            else 
                f = 0.0
            end

            if f > 0.0
                push!(projection_intersected_cells, cell_id)
                push!(projection_volume_in_cells, f * cell_volume)
                push!(projection_distances_from_source, norm(cell_centroid - source_transducer_coordinates))
            end
        else
            # For narrow cones, check if the cone axis (ray) itself intersects the cell.
            n_faces = nfacets(grid.cells[cell_id])
            t_ray_entry = 0.0
            t_ray_exit = projection_distance
            ray_intersect = true
            for face_idx in 1:n_faces
                N_f = geo.cell_face_normals[cell_id][face_idx]
                face_nodes = get_face_nodes(grid, cell_id, face_idx)
                C_f = grid.nodes[face_nodes[1]].x
                
                D = dot(unit_ray_vector, N_f)
                N = dot(C_f - source_transducer_coordinates, N_f)
                
                if abs(D) < 1e-12
                    if N < -1e-12
                        ray_intersect = false
                        break
                    end
                elseif D > 0.0
                    t_ray_exit = min(t_ray_exit, N / D)
                else
                    t_ray_entry = max(t_ray_entry, N / D)
                end
            end
            
            if ray_intersect && (t_ray_exit - t_ray_entry > 1e-10)
                # Calculate frustum volume of the narrow cone segment
                v_frustum = (pi * (tand(angle_of_projection)^2) / 3.0) * (t_ray_exit^3 - t_ray_entry^3)
                intersect_volume = min(cell_volume, v_frustum)
                if intersect_volume > 0.0
                    push!(projection_intersected_cells, cell_id)
                    push!(projection_volume_in_cells, intersect_volume)
                    push!(projection_distances_from_source, norm(cell_centroid - source_transducer_coordinates))
                end
            end
        end
    end

    for cell_id in projection_intersected_cells
        cell_displacement_vector = source_transducer_coordinates .- geo.cell_centroids[cell_id]
        push!(projection_unit_vectors_to_cell_centers, normalize(cell_displacement_vector))
    end

    # Sort results by distance from source transducer
    if !isempty(projection_distances_from_source)
        perm = sortperm(projection_distances_from_source)
        projection_intersected_cells = projection_intersected_cells[perm]
        projection_volume_in_cells = projection_volume_in_cells[perm]
        projection_unit_vectors_to_cell_centers = projection_unit_vectors_to_cell_centers[perm]
    end

    return projection_intersected_cells, projection_volume_in_cells, projection_distances_from_source, projection_distance, slant_distance, projection_unit_vectors_to_cell_centers
end
