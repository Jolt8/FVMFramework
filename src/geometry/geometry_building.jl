abstract type FVMGeometry end

struct FVMGeometryTetra{T,CoordType} <: FVMGeometry
    cell_neighbor_map::Vector{Tuple{Int, Int}}
    cell_volumes::Vector{T}
    cell_centroids::Vector{CoordType}
    connection_areas::Vector{T}
    connection_normals::Vector{CoordType}
    connection_distances::Vector{T}
    unconnected_areas::Vector{SVector{4, T}}
    unconnected_normals::Vector{SVector{4, CoordType}}
end

struct FVMGeometryHexa{T,CoordType} <: FVMGeometry
    cell_neighbor_map::Vector{Tuple{Int, Int}}
    cell_volumes::Vector{T}
    cell_centroids::Vector{CoordType}
    connection_areas::Vector{T}
    connection_normals::Vector{CoordType}
    connection_distances::Vector{T}
    unconnected_areas::Vector{SVector{6, T}}
    unconnected_normals::Vector{SVector{6, CoordType}}
end

function build_fvm_geo_into_struct(grid)
    initial_node_coordinates = get_node_coordinates(grid)

    cell_neighbor_map, neighbor_map_respective_node_ids = get_neighbor_map(grid)

    unconnected_cell_face_map, unconnected_map_respective_node_ids = get_unconnected_map(grid)

    nodes_of_cells = get_nodes_of_cells(grid)

    cell_volumes, cell_centroids, #cell volumes and cell centroids are accessed at the id of the cell
    connection_areas, connection_normals, connection_distances, #connection areas, normals, and distances are simply accessed by their location in the list which corresponds to the respective connection in cell_neighbor_map
    unconnected_areas, unconnected_normals = rebuild_fvm_geometry(
        cell_neighbor_map, neighbor_map_respective_node_ids,
        unconnected_cell_face_map, unconnected_map_respective_node_ids,
        initial_node_coordinates, nodes_of_cells
    )

    if typeof(grid) == Grid{3, Tetrahedron, Float64}
        return FVMGeometryTetra(
            cell_neighbor_map,
            cell_volumes, cell_centroids,
            connection_areas, connection_normals, connection_distances,
            unconnected_areas, unconnected_normals
        )
    elseif typeof(grid) == Grid{3, Hexahedron, Float64}
        return FVMGeometryHexa(
            cell_neighbor_map,
            cell_volumes, cell_centroids,
            connection_areas, connection_normals, connection_distances,
            unconnected_areas, unconnected_normals
        )
    end
end
