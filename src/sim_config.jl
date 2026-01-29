struct FVMGeometry{T, CoordType}
    cell_neighbor_map::Vector{Tuple{Int, Int}}
    cell_volumes::Vector{T}
    cell_centroids::Vector{CoordType}
    connection_areas::Vector{T}
    connection_normals::Vector{CoordType}
    connection_distances::Vector{T}
    unconnected_areas::Vector{SVector{6, T}}
    unconnected_normals::Vector{SVector{6, CoordType}}
end

mutable struct RegionSetupInfo
    name::String
    physics::AbstractPhysics
    initial_conditions::ComponentArray
    cell_set_idxs::Set{Int}
end

mutable struct BoundarySetupInfo
    name::String
    fixed_conditions::ComponentArray
    cell_set_idxs::Set{Int}
end

struct SimulationConfigInfo
    grid::Ferrite.Grid
    regions::Vector{RegionSetupInfo}
    boundaries::Vector{BoundarySetupInfo}
    u_proto::ComponentArray
    u_axes::Axis
    geo::FVMGeometry
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

    return FVMGeometry(
        cell_neighbor_map,
        cell_volumes, cell_centroids,
        connection_areas, connection_normals, connection_distances,
        unconnected_areas, unconnected_normals
    )
end

function create_fvm_config(grid, u_proto)
    u_axes = getaxes(u_proto)[1]

    geo = build_fvm_geo_into_struct(grid)

    return SimulationConfigInfo(
        grid,
        RegionSetupInfo[],
        BoundarySetupInfo[],
        u_proto,
        u_axes,
        geo
    )
end

#if you're wondering what the type of u_axes is:
#Axis{(vel_x = ViewAxis(1:1, Shaped1DAxis((1,))), vel_y = ViewAxis(2:2, Shaped1DAxis((1,))), vel_z = ViewAxis(3:3, Shaped1DAxis((1,))), 
#pressure = ViewAxis(4:4, Shaped1DAxis((1,))), mass_fractions = ViewAxis(5:9, ShapedAxis((5, 1))), temp = ViewAxis(10:10, Shaped1DAxis((1,))))}

function add_region!(
    config, name;
    physics,
    initial_conditions
    )

    cell_set_idxs = Set(getcellset(config.grid, name))

    for cell_id in cell_set_idxs
        config.u_proto.vel_x[cell_id] = initial_conditions.vel_x
        config.u_proto.vel_y[cell_id] = initial_conditions.vel_y
        config.u_proto.vel_z[cell_id] = initial_conditions.vel_z
        config.u_proto.pressure[cell_id] = initial_conditions.pressure
        config.u_proto.temp[cell_id] = initial_conditions.temp
        config.u_proto.mass_fractions[:, cell_id] = initial_conditions.mass_fractions
    end

    push!(config.regions, RegionSetupInfo(name, physics, initial_conditions, cell_set_idxs))
end

function add_boundary!(
    config, name;
    fixed_conditions
    )
    cell_set_idxs = Set(getcellset(config.grid, name))

    push!(config.boundaries, BoundarySetupInfo(name, fixed_conditions, cell_set_idxs))
end

struct FVMSystem
    phys::Vector{AbstractPhysics}
    cell_phys_id_map::Vector{Int}
    fixed_idxs_and_vals_map::ComponentArray
    u_axes::Axis
end

function finish_fvm_config(config)
    n_cells = length(config.geo.cell_volumes)

    cell_phys_id_map = zeros(Int, n_cells)

    phys = AbstractPhysics[]

    for (i, region) in enumerate(config.regions)
        push!(phys, region.physics)

        for cell_id in region.cell_set_idxs
            #we could probably do something like if the same type of physics has already been used, 
            #assign it to the previous i, but I think that's an unncessary optimizaiton for now
            cell_phys_id_map[cell_id] = i
        end
    end
    
    property_names = propertynames(config.u_proto)
    n_props = length(property_names)
    containers = [Tuple{Int, Float64}[] for _ in 1:n_props]

    fixed_idxs_and_vals_map = ComponentArray(NamedTuple{property_names}(containers))

    for (i, boundary) in enumerate(config.boundaries)
        for cell_id in boundary.cell_set_idxs
            for field in propertynames(boundary.fixed_conditions)
                if field in property_names
                    push!(fixed_idxs_and_vals_map[field], (cell_id, boundary.fixed_conditions[field]))
                end
            end
        end
    end
    #we might have to flatten fixed_idxs_and_vals_map

    u0 = Vector(config.u_proto)
    du0 = u0 .* 0.0

    return du0, u0, config.geo, FVMSystem(
        phys,
        cell_phys_id_map,
        fixed_idxs_and_vals_map,
        config.u_axes,
    )
end