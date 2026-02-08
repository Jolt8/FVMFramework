#MAJOR CHANGES ARE ABOUT TO BE MADE TO SWITCH TO CELL CENTRIC CONNECTIONS

#note for the future to differentiate between Tetras and Hexas:
#typeof(grid) == Grid{3, Tetrahedron, Float64}

#typeof(grid.cells[1]) == Tetrahedron 

function get_node_coordinates(grid)
    n_nodes = length(grid.nodes)
    
    node_coordinates = Vector{SVector{3, Float64}}(undef, n_nodes)
    
    visited_map = Set{Int}()

    for node_id in 1:length(grid.nodes)
        if !(node_id in visited_map)
            #add it to the set so its coordinates do not get re-added
            push!(visited_map, node_id)
            coordinates = get_node_coordinate(grid.nodes[node_id])
            
            node_coordinates[node_id] = SVector(coordinates[1], coordinates[2], coordinates[3])
        end
    end
    return node_coordinates
end

#I think we could safely get rid of this one
function get_cell_topology(grid)
    n_nodes = length(grid.nodes)

    topology = Vector{Vector{Int}}(undef, n_nodes)

    visited_map = Set()

    for cell in CellIterator(grid)
        cell_id = cellid(cell)
        for i in eachindex(cell.nodes)
            node_id = cell.nodes[i]
            if !(node_id in visited_map)
                push!(visited_map, node_id) #add it to the set so its coordinates do not get re-added

                topology[node_id] = [cell_id]
            else
                push!(topology[node_id], cell_id)
            end
        end
    end
    return topology
end

function get_nodes_of_cells(grid)
    n_cells = length(grid.cells)

    node_topology = Vector{Vector{Int}}(undef, n_cells)

    for cell in CellIterator(grid)
        cell_id = cellid(cell)
        node_topology[cell_id] = collect(grid.cells[cell_id].nodes)
    end
    
    return node_topology
end

function get_face_nodes(grid, cell_id, face_idx)
    cell = grid.cells[cell_id]
    face_nodes = collect(Ferrite.faces(cell)[face_idx])
    n_face_nodes = length(face_nodes)
    return SVector{n_face_nodes, Int}(face_nodes)
end

#TODO: implement a method to add a "preferrability" input vector of types that makes sure that a certain type always comes before another
#this would be for avoiding dynamic dispatch when getting rho for example by guaranteeing that fluid cells come before solid cells 
function get_neighbor_map(grid)
    #returns:
    #cell/neighbor pairs (ex. (1, 2) (cell_1 is connected with cell_2))
    #respective node idxs of cell and neighbor (ex. (2, 8, 44, 38))

    cell_neighbor_pairs = Tuple{Int, Int}[]
    
    n_facets_per_cell = length(grid.cells[1].facets)
    n_nodes_per_face = length(get_face_nodes(grid, 1, 1))
    neighbor_map_respective_node_ids = NTuple{n_nodes_per_face, Int}[]

    n_cells = length(grid.cells)

    top = ExclusiveTopology(grid) #this causes A LOT of GC and slowdown, we may want to access it once
    for cell_id in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[cell_id])
            neighbor_info = top.face_face_neighbor[cell_id, face_idx]

            if !isempty(neighbor_info)
                neighbor_id = collect(neighbor_info[1].idx)[1]
                
                if cell_id < neighbor_id
                    respective_nodes = get_face_nodes(grid, cell_id, face_idx)
                    
                    push!(cell_neighbor_pairs, (cell_id, neighbor_id))
                    push!(neighbor_map_respective_node_ids, respective_nodes)
                end
            end
        end
    end
    return cell_neighbor_pairs, neighbor_map_respective_node_ids
end

#just read some of Ferrite's source code and automatically defining the cellset 
#depending on whether the whole grid or just a cellset is passed in is very smart
#I'm just worried if nfacets works on a cellset
function get_cell_neighbors(grid::AbstractGrid, top; cellset = grid.cells) 
    #returns:
    #cells and their respective neighbor ids at each of their face_idxs
    #ex. [(0, 50, 123, 99)...] (cell 1 has no neighbor at face 1, neighbors cell 50 at face 2, etc...)

    n_cells = getncells(grid)

    n_facets_per_cell = length(grid.cells[1].facets)
    n_nodes_per_face = length(get_face_nodes(grid, 1, 1))

    cell_neighbors = Vector{MVector{n_facets_per_cell, Int}}(undef, n_cells)
    #[Int[0, 0, 0, 0]...] for tetras for example
    neighbor_map_respective_node_ids = Vector{MVector{n_facets_per_cell, SVector{n_nodes_per_face, Int}}}(undef, n_cells)
    #this is pretty complicated, each cell has to store a vector of vector of vectors for the cell_id, facet_idx, and node_id

    for cell_id in cellset
        for face_idx in 1:nfacets(grid.cells[cell_id])
            neighbor_info = top.face_face_neighbor[cell_id, face_idx]

            if !isempty(neighbor_info)
                neighbor_id = collect(neighbor_info[1].idx)[1]

                cell_neighbors[cell_id][face_idx] = neighbor_id
                neighbor_map_respective_node_ids[cell_id][face_idx] = get_face_nodes(grid, cell_id, face_idx)
            end
        end
    end
    return cell_neighbors, neighbor_map_respective_node_ids
end

function get_unconnected_map(grid, top; cellset = grid.cells)
    nodes_of_cells = get_nodes_of_cells(grid)
    #returns:
    #list of unconnected faces for each cell (ex. (1, 5) (cell_1's 5th face_idxs is not connected))
    #respective node idxs of cell face (ex. (2, 8, 44, 38))

    n_cells = getncells(grid)
    n_nodes_per_face = length(get_face_nodes(grid, 1, 1))

    unconnected_cell_face_map = NTuple{2, Int}[]
    unconnected_map_respective_node_ids = NTuple{n_nodes_per_face, Int}[]

    for cell_id in cellset
        for face_idx in 1:nfacets(grid.cells[cell_id])
            neighbor_info = top.face_face_neighbor[cell_id, face_idx]

            if isempty(neighbor_info) #note that were checking if it isempty here, not !isempty like above
                curr_cell_nodes = nodes_of_cells[cell_id]
                respective_nodes = get_face_nodes(grid, cell_id, face_idx)

                push!(unconnected_cell_face_map, (cell_id, face_idx))
                push!(unconnected_map_respective_node_ids, ntuple(i -> respective_nodes[i], n_nodes_per_face))
            end
        end
    end
    return unconnected_cell_face_map, unconnected_map_respective_node_ids
end


function get_cell_face_map(grid)
    #returns:
    #list of faces for each cell (ex. (1, 5) (cell_1's 5th face_idx))
    #respective node idxs of cell face (ex. (2, 8, 44, 38))

    n_cells = length(grid.cells)
    n_nodes_per_face = length(get_face_nodes(grid, 1, 1))

    cell_face_map = NTuple{2, Int}[]
    map_respective_node_ids = NTuple{n_nodes_per_face, Int}[]

    for cell_id in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[cell_id])
            respective_nodes = get_face_nodes(grid, cell_id, face_idx)

            push!(cell_face_map, (cell_id, face_idx))
            push!(map_respective_node_ids, ntuple(i -> respective_nodes[i], n_nodes_per_face))
        end
    end
    return cell_face_map, map_respective_node_ids
end


function cross_product(a, b)
    x = a[2]*b[3] - a[3]*b[2]
    y = a[3]*b[1] - a[1]*b[3]
    z = a[1]*b[2] - a[2]*b[1]
    return SVector(x, y, z)
end
