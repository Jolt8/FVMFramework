mutable struct RegionSetupInfo#{T<:Type{AbstractPhysics}} #how the heck can we accept any type that's a subset of AbstractPhysics
    name::String
    type::AbstractPhysics #this is used to define connections
    initial_conditions::ComponentArray
    properties::ComponentArray
    region_function::Function
    region_cells::Vector{Int}
end

mutable struct ControllerSetupInfo #this must be defined before SimulationConfigInfo
    controller::ComponentArray
    monitored_cellset::String
    affected_cellset::String
    controller_function::Function
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

struct SimulationConfigInfo
    grid::Ferrite.Grid
    regions::Vector{RegionSetupInfo}
    controllers::Vector{ControllerSetupInfo}
    u_proto::ComponentArray
    geo::FVMGeometry
end

function create_fvm_config(grid, u_proto)
    geo = build_fvm_geo_into_struct(grid)

    return SimulationConfigInfo(
        grid,
        RegionSetupInfo[],
        ControllerSetupInfo[],
        u_proto,
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
    initial_conditions,
    properties,
    region_function
)
    region_cells = collect(getcellset(config.grid, name))

    for cell_id in region_cells
        for field in propertynames(initial_conditions)
            if length(propertynames(initial_conditions[field])) > 1
                #FIXME: we need to overload this so u.u_proto[field] = something actually works 
                for field_2 in propertynames(initial_conditions[field])
                    store = config.u_proto[field][field_2][cell_id]
                    #println(config.u_proto[field][field_2][cell_id])
                    #we should probably overload the basic getproperty / getindex with something that actually returns a value that mutates the original
                    #config.u_proto[field][field_2][cell_id] = initial_conditions[field][field_2] 
                    getproperty(getproperty(config.u_proto, field), field_2)[cell_id] = initial_conditions[field][field_2]
                    #println(config.u_proto[field][field_2][cell_id])
                    #println(store == config.u_proto[field][field_2][cell_id])
                end
            else
                getproperty(config.u_proto, field)[cell_id] = initial_conditions[field]
            end
        end
    end

    push!(config.regions, RegionSetupInfo(name, type, initial_conditions, properties, region_function, region_cells))
end

#this could probably also be handled by dynamic dispatch for facets, but it helps the user know a different routine is happening
#this needs to be updated, it currently doesn't work
function add_facet_region!(
    config, name;
    initial_conditions,
    region_physics,
    region_function,
)

    region_cells = get_cell_ids_in_facet_set(config.grid, name)

    for cell_id in region_cells
        for field in propertynames(initial_conditions)
            if ndims(getproperty(config.u_proto, field)) > 1 #for mass fractions
                getproperty(config.u_proto, field)[:, cell_id] = initial_conditions[field]
            else
                getproperty(config.u_proto, field)[cell_id] = initial_conditions[field]
            end
        end
    end

    push!(config.regions, RegionSetupInfo(name, type, initial_conditions, region_physics, region_function, region_cells))
end

function add_controller!(
    config;
    controller_properties,
    controller_function,
    monitored_cellset,
    affected_cellset
)
    monitored_cells = collect(getcellset(config.grid, monitored_cellset))
    affected_cells = collect(getcellset(config.grid, affected_cellset))

    push!(config.controllers, ControllerSetupInfo(controller_properties, monitored_cellset, affected_cellset, controller_function, monitored_cells, affected_cells))
end

struct ConnectionGroup{F<:Function}
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
    controller_function!::F
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

struct RegionGroup{F<:Function}
    region_function!::F
    region_cells::Vector{Int}
end


struct FVMSystem
    connection_groups::Vector{ConnectionGroup}
    controller_groups::Vector{ControllerGroup}
    region_groups::Vector{RegionGroup}
    u_merged_axes::Axis
end

function finish_fvm_config(config, conneciton_map_function)
    n_cells = length(config.geo.cell_volumes)
    
    cell_phys_map = Vector{AbstractPhysics}(undef, n_cells)

    for (i, region) in enumerate(config.regions)
        @batch for cell_id in region.region_cells
            cell_phys_map[cell_id] = region.type
        end
    end

    connection_groups = ConnectionGroup[]

    unique_celltype_pairs = Vector{Tuple{AbstractPhysics, AbstractPhysics}}()

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
                    flux_function!,
                    [(idx_a, Tuple{Int,Int}[]) for idx_a in 1:n_cells]
                )
                )
                push!(connection_groups[connection_group_id].connections[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            elseif (phys_a, phys_b) in unique_celltype_pairs && isempty(connection_groups[connection_group_id].connections[idx_a])
                println("I'm unsure if this actually ever happens, check to line 221 in sim config if it does")
                connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)
                push!(connection_groups[connection_group_id].connections[idx_a], (
                    (idx_a, Vector{Tuple{Int,Int}}((idx_b, face_idx)))
                )
                )
            elseif !isempty(connection_groups[connection_group_id].connections[idx_a])
                connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)
                push!(connection_groups[connection_group_id].connections[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            end
        end
    end

    #remove empty tuples from cell_neighbors
    for CG in connection_groups
        for (idx_a, idx_a_neighbors) in CG.cell_neighbors
            cell_neighbors[idx_a] = filter(!isempty, cell_neighbors[idx_a])
        end
    end

    region_groups = RegionGroup[]

    for region in config.regions
        push!(region_groups, RegionGroup(region.region_function, region.region_cells))
    end

    controller_groups = ControllerGroup[]

    for (controller_id, controller) in enumerate(config.controllers)
        push!(controller_groups, ControllerGroup(
                controller_id,
                controller.controller_function,
                controller.monitored_cells, controller.affected_cells
            )
        )
    end

    #this code is to merge the state variables such as velocity, temperature or pressure with non-state variables such as integral_error
    #integral error is indexed at its respective controller id
    n_controllers = length(config.controllers)
    u_non_state_proto = ComponentArray(
        integral_error=zeros(n_controllers)
    )
    #unlike u_proto, u_non_state_proto is accessed at controller_id

    u_merged = ComponentVector(config.u_proto; NamedTuple(u_non_state_proto)...)
    u_merged_axes = getaxes(u_merged)[1]

    u0 = Vector(u_merged)
    du0 = u0 .* 0.0

    return du0, u0, config.geo, FVMSystem(
        connection_groups, controller_groups, region_groups,
        u_merged_axes,
    )
end