using Ferrite
using LinearAlgebra
using FVMFramework

struct BeamProfile
    angle_of_projection::Float64      # degrees (far-field divergence half-angle)
    near_field_distance::Float64      # same units as grid coordinates (Fresnel zone length)
    transducer_radius::Float64        # same units as grid coordinates
end

# Check if a point is inside the beam envelope (near-field cylinder + far-field cone)
function is_point_in_beam(P, source_transducer_coordinates, unit_ray_vector, projection_distance, beam::BeamProfile)
    v = P - source_transducer_coordinates
    d_axis = dot(v, unit_ray_vector)
    if d_axis < 0.0 || d_axis > projection_distance
        return false
    end
    
    # Compute perpendicular distance from beam axis
    v_parallel = d_axis * unit_ray_vector
    dist_from_axis = norm(v - v_parallel)
    
    if d_axis <= beam.near_field_distance
        # Near-field (Fresnel zone): cylindrical beam
        return dist_from_axis <= beam.transducer_radius
    else
        # Far-field (Fraunhofer zone): diverging cone from the near-field boundary
        far_field_distance = d_axis - beam.near_field_distance
        max_radius = beam.transducer_radius + far_field_distance * tand(beam.angle_of_projection)
        return dist_from_axis <= max_radius
    end
end

# Compute the beam volume for a segment between d_entry and d_exit along the axis
function beam_segment_volume(d_entry, d_exit, beam::BeamProfile)
    N = beam.near_field_distance
    R = beam.transducer_radius
    tan_theta = tand(beam.angle_of_projection)
    
    vol = 0.0
    
    # Near-field portion (cylinder)
    nf_start = max(d_entry, 0.0)
    nf_end = min(d_exit, N)
    if nf_end > nf_start
        vol += pi * R^2 * (nf_end - nf_start)
    end
    
    # Far-field portion (expanding frustum)
    ff_start = max(d_entry, N)
    ff_end = d_exit
    if ff_end > ff_start
        r1 = R + (ff_start - N) * tan_theta
        r2 = R + (ff_end - N) * tan_theta
        vol += pi * (ff_end - ff_start) / 3.0 * (r1^2 + r1 * r2 + r2^2)
    end
    
    return vol
end


function get_cells_in_transducer_projection(
    source_transducer_node_id, opposing_transducer_node_id, grid, geo,
    beam::BeamProfile
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

    projection_distance = norm(ray_vector)

    slant_distance = projection_distance / cosd(beam.angle_of_projection)

    projection_intersected_cells = Int[]
    projection_volume_in_cells = Float64[]
    projection_distances_from_source = Float64[]
    projection_unit_vectors_to_cell_centers = Vector{Float64}[]

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
            if is_point_in_beam(p[i], source_transducer_coordinates, unit_ray_vector, projection_distance, beam)
                vertices_inside += 1
            end
        end

        centroid_inside = is_point_in_beam(cell_centroid, source_transducer_coordinates, unit_ray_vector, projection_distance, beam)

        if vertices_inside == n_nodes
            # Entirely inside the beam
            push!(projection_intersected_cells, cell_id)
            push!(projection_volume_in_cells, cell_volume)
            push!(projection_distances_from_source, norm(cell_centroid - source_transducer_coordinates))
        elseif vertices_inside > 0 || centroid_inside
            # Partially inside the beam, estimate volume via sub-cell sampling
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
                    if is_point_in_beam(pt, source_transducer_coordinates, unit_ray_vector, projection_distance, beam)
                        n_inside += 1
                    end
                end
            elseif CellType <: Tetrahedron
                for r in 0.1:0.2:0.9, s in 0.1:0.2:0.9, t in 0.1:0.2:0.9
                    if r + s + t <= 0.95
                        n_total += 1
                        pt = (1 - r - s - t) * p[1] + r * p[2] + s * p[3] + t * p[4]
                        if is_point_in_beam(pt, source_transducer_coordinates, unit_ray_vector, projection_distance, beam)
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
            # For narrow beams, check if the beam axis (ray) itself intersects the cell.
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
                # Calculate beam volume of the segment using the two-zone model
                intersect_volume = min(cell_volume, beam_segment_volume(t_ray_entry, t_ray_exit, beam))
                if intersect_volume > 0.0
                    push!(projection_intersected_cells, cell_id)
                    push!(projection_volume_in_cells, intersect_volume)
                    push!(projection_distances_from_source, norm(cell_centroid - source_transducer_coordinates))
                end
            end
        end
    end

    for (transducer_cell_idx, cell_id) in enumerate(projection_intersected_cells)
        cell_displacement_vector = geo.cell_centroids[cell_id] .- source_transducer_coordinates
        push!(projection_unit_vectors_to_cell_centers, normalize(cell_displacement_vector))
    end

    # Sort results by distance from source transducer
    if !isempty(projection_distances_from_source)
        perm = sortperm(projection_distances_from_source)
        projection_distances_from_source = projection_distances_from_source[perm]
        projection_intersected_cells = projection_intersected_cells[perm]
        projection_volume_in_cells = projection_volume_in_cells[perm]
        projection_unit_vectors_to_cell_centers = projection_unit_vectors_to_cell_centers[perm]
    end

    return projection_intersected_cells, projection_volume_in_cells, projection_distances_from_source, projection_distance, slant_distance, projection_unit_vectors_to_cell_centers
end
