mutable struct RegionSetupInfo#{T<:Type{AbstractPhysics}} #how the heck can we accept any type that's a subset of AbstractPhysics
    name::String
    type::AbstractPhysics #this is used to define connections
    region_function::Function
    region_cells::Vector{Int}
end

mutable struct ControllerSetupInfo #this must be defined before SimulationConfigInfo
    id::Int
    name::String
    monitored_name::String
    affected_name::String
    controller_function::Function
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

struct SimulationConfigInfo
    grid::Ferrite.Grid
    regions::Vector{RegionSetupInfo}
    controllers::Vector{ControllerSetupInfo}
    geo::FVMGeometry
end

function create_fvm_config(grid)
    geo = build_fvm_geo_into_struct(grid)

    return SimulationConfigInfo(
        grid,
        RegionSetupInfo[],
        ControllerSetupInfo[],
        geo
    )
end

#if you're wondering what the type of u_axes is:
#Axis{(vel_x = ViewAxis(1:1, Shaped1DAxis((1,))), vel_y = ViewAxis(2:2, Shaped1DAxis((1,))), vel_z = ViewAxis(3:3, Shaped1DAxis((1,))), 
#pressure = ViewAxis(4:4, Shaped1DAxis((1,))), mass_fractions = ViewAxis(5:9, ShapedAxis((5, 1))), temp = ViewAxis(10:10, Shaped1DAxis((1,))))}

#VERY IMPORTANT!!!!
#= For future reference when getting properties using u[field]:
    u_cv.temp[cell_id] = val
        works!
    u_cv[field][cell_id] = val
        fails :(
    view(u_cv, field)[cell_id] = val
        works!
    getproperty(u_cv, field)[cell_id] = val
        works!
=#

function add_region!(
    config, name;
    type,
    region_function
)
    region_cells = collect(getcellset(config.grid, name))

    push!(config.regions, RegionSetupInfo(name, type, region_function, region_cells))
end

#this could probably also be handled by dynamic dispatch for facets, but it helps the user know a different routine is happening
#this needs to be updated, it currently doesn't work
function add_facet_region!(
    config, name;
    region_function,
)

    region_cells = get_cell_ids_in_facet_set(config.grid, name)

    push!(config.regions, RegionSetupInfo(name, type, region_function, region_cells))
end

function add_controller!(
    config, name;
    controller_function,
    monitored_cellset,
    affected_cellset
)
    monitored_cells = collect(getcellset(config.grid, monitored_cellset))
    affected_cells = collect(getcellset(config.grid, affected_cellset))

    controller_id = length(config.controllers) + 1

    push!(config.controllers, ControllerSetupInfo(controller_id, name, monitored_cellset, affected_cellset, controller_function, monitored_cells, affected_cells))
end

#the names here are only for the tracer 

struct ConnectionGroup{F<:Function}
    name_a::String
    name_b::String
    flux_function!::F
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}} #we could probably use ::Vector{Connection} here
end

#=
struct Connection
    idx_a::Int
    idx_b_face_idx_vec::Vector{Tuple{Int, Int}}
end
=#

struct ControllerGroup{F <:Function}
    id::Int
    name::String
    monitored_name::String
    affected_name::String
    controller_function!::F
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

struct RegionGroup{F<:Function}
    name::String
    region_function!::F
    region_cells::Vector{Int}
end


struct FVMSystem
    connection_groups::Vector{ConnectionGroup}
    controller_groups::Vector{ControllerGroup}
    region_groups::Vector{RegionGroup}
end

function finish_fvm_config(config, connection_map_function)
    n_cells = length(config.geo.cell_volumes)
    
    cell_phys_map = Vector{AbstractPhysics}(undef, n_cells)
    cell_region_name_map = Vector{String}(undef, n_cells)

    for (i, region) in enumerate(config.regions)
        @batch for cell_id in region.region_cells
            cell_phys_map[cell_id] = region.type
            cell_region_name_map[cell_id] = region.name
        end
    end

    connection_groups = ConnectionGroup[]

    unique_celltype_pairs = Vector{Tuple{AbstractPhysics, AbstractPhysics}}()
    unique_region_connection_pairs = Vector{Tuple{String, String}}()

    for (idx_a, idx_a_neighbors) in config.geo.cell_neighbors
        for (idx_b, face_idx) in idx_a_neighbors
            if idx_b <= 0
                continue
            end
            phys_a = cell_phys_map[idx_a]
            phys_b = cell_phys_map[idx_b]

            region_name_a = cell_region_name_map[idx_a]
            region_name_b = cell_region_name_map[idx_b]

            flux_function! = connection_map_function(typeof(phys_a), typeof(phys_b))

            connection_group_id = findfirst(item -> item == (region_name_a, region_name_b), unique_region_connection_pairs)

            if !((region_name_a, region_name_b) in unique_region_connection_pairs)
                push!(unique_region_connection_pairs, (region_name_a, region_name_b))
                connection_group_id = findfirst(item -> item == (region_name_a, region_name_b), unique_region_connection_pairs)
                push!(connection_groups, ConnectionGroup(
                    region_name_a,
                    region_name_b,
                    flux_function!,
                    [(idx_a, Tuple{Int,Int}[]) for idx_a in 1:n_cells]
                )
                )
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            elseif (region_name_a, region_name_b) in unique_region_connection_pairs && isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                println("I'm unsure if this actually ever happens, check to line 221 in sim config if it does")
                connection_group_id = findfirst(item -> item == (region_name_a, region_name_b), unique_region_connection_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a], (
                    (idx_a, Vector{Tuple{Int,Int}}((idx_b, face_idx)))
                )
                )
            elseif !isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                connection_group_id = findfirst(item -> item == (region_name_a, region_name_b), unique_region_connection_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            end
        end
    end

    #= This should probably be used for the actual main groups because it creates fewer unique groups, but it's bad for the tracer
    for (idx_a, idx_a_neighbors) in config.geo.cell_neighbors
        for (idx_b, face_idx) in idx_a_neighbors
            if idx_b <= 0
                continue
            end
            phys_a = cell_phys_map[idx_a]
            phys_b = cell_phys_map[idx_b]

            flux_function! = conneciton_map_function(typeof(phys_a), typeof(phys_b))

            connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)

            if !((phys_a, phys_b) in unique_celltype_pairs)
                push!(unique_celltype_pairs, (phys_a, phys_b))
                connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)
                push!(connection_groups, ConnectionGroup(
                    cell_region_name_map[idx_a],
                    cell_region_name_map[idx_b],
                    flux_function!,
                    [(idx_a, Tuple{Int,Int}[]) for idx_a in 1:n_cells]
                )
                )
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            elseif (phys_a, phys_b) in unique_celltype_pairs && isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                println("I'm unsure if this actually ever happens, check to line 221 in sim config if it does")
                connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a], (
                    (idx_a, Vector{Tuple{Int,Int}}((idx_b, face_idx)))
                )
                )
            elseif !isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            end
        end
    end
    =#
    

    #remove empty tuples from cell_neighbors
    
    for CG in connection_groups
        for (idx_a, idx_a_neighbors) in CG.cell_neighbors
            filter!(conn -> !(isempty(conn)), CG.cell_neighbors)
        end
    end
    
    region_groups = RegionGroup[]

    for region in config.regions
        push!(region_groups, RegionGroup(region.name, region.region_function, region.region_cells))
    end

    controller_groups = ControllerGroup[]

    for (controller_id, controller) in enumerate(config.controllers)
        push!(controller_groups, ControllerGroup(
                controller_id,
                controller.name,
                controller.monitored_name,
                controller.affected_name,
                controller.controller_function,
                controller.monitored_cells, controller.affected_cells
            )
        )
    end

    return config.geo, FVMSystem(connection_groups, controller_groups, region_groups)
end