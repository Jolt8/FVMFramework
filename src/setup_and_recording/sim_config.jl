mutable struct RegionSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    initial_conditions::ComponentArray
    region_properties::ComponentArray
    region_function::Function
    region_cells::Vector{Int}
end

abstract type AbstractController end

mutable struct ControllerSetupInfo #this must be defined before SimulationConfigInfo
    controller::AbstractController
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
    initial_conditions,
    region_properties,
    region_function,
)

    region_cells = collect(getcellset(config.grid, name))

    for cell_id in region_cells
        for field in propertynames(initial_conditions)
            if ndims(getproperty(config.u_proto, field)) > 1 #for mass fractions
                getproperty(config.u_proto, field)[:, cell_id] = initial_conditions[field]
            else
                getproperty(config.u_proto, field)[cell_id] = initial_conditions[field]
            end
        end
    end

    push!(config.regions, RegionSetupInfo(name, initial_conditions, region_properties, region_function, region_cells))
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

    push!(config.regions, RegionSetupInfo(name, initial_conditions, region_physics, region_function, region_cells))
end

function add_controller!(
    config;
    controller,
    monitored_cellset,
    affected_cellset,
    controller_function
)
    monitored_cells = collect(getcellset(config.grid, monitored_cellset))
    affected_cells = collect(getcellset(config.grid, affected_cellset))

    push!(config.controllers, ControllerSetupInfo(controller, monitored_cellset, affected_cellset, controller_function, monitored_cells, affected_cells))
end


#TODO: allow for face boundary conditions, boundary conditions applied on faces only apply to the faces right now

#we should probably find a way to make the creation of these structs automatic 
struct RegionGroup{T <: AbstractPhysics, F <: Function}
    phys::T
    region_function!::F
    region_cells::Vector{Int}
end

struct ControllerGroup{T <: AbstractController, F <: Function}
    id::Int
    controller::T
    controller_function!::F
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

#=
struct Connection
    idx_a::Int
    idx_b_face_idx_vec::Vector{Tuple{Int, Int}}
end
=#

struct ConnectionGroup{TA <: AbstractPhysics, TB <: AbstractPhysics, F <: Function}
    phys_a::TA
    phys_b::TB
    flux_function!::F
    connections::Vector{Tuple{Int, Tuple{Int, Int}}} #we could probably use ::Vector{Connection} here
end


struct FVMSystem
    phys::Vector{AbstractPhysics}
    cell_phys_id_map::Vector{Int}
    connection_groups::Vector{ConnectionGroup}
    controller_groups::Vector{ControllerGroup}
    region_groups::Vector{RegionGroup}
    u_merged_axes::Axis
end

function finish_fvm_config(config, conneciton_map_function)
    n_cells = length(config.geo.cell_volumes)

    phys = AbstractPhysics[]
    cell_phys_id_map = zeros(Int, n_cells)
    #the only reason this exists is because we still need to access the physics of the cell in the connections loop
    #we're doubling up on physics, but it makes everything cleaner later

    region_groups = RegionGroup[]

    for (i, region) in enumerate(config.regions)
        push!(region_structs, RegionGroup(region.region_physics, region.region_function, region.region_cells))
        push!(phys, region.region_physics)

        @batch for cell_id in region.region_cells
            cell_phys_id_map[cell_id] = i
        end
    end

    controller_groups = ControllerGroup[]

    for (controller_id, controller) in enumerate(config.controllers)
        push!(controller_groups, ControllerGroup(
            controller_id,
            controller.controller,
            controller.controller_function,
            controller.monitored_cells, controller.affected_cells
        )
        )
    end

    connection_groups = ConnectionGroup[]

    unique_celltype_pairs = Vector{Tuple{AbstractPhysics, AbstractPhysics}}()

    for (idx_a, idx_a_neighbors) in config.geo.cell_neighbors
        for (idx_b, face_idx) in idx_a_neighbors
            if idx_b <= 0
                continue
            end
            phys_a = phys[cell_phys_id_map[idx_a]]
            phys_b = phys[cell_phys_id_map[idx_b]]

            flux_function! = conneciton_map_function(phys_a, phys_b)

            if !((phys_a, phys_b) in unique_celltype_pairs)
                push!(unique_celltype_pairs, (phys_a, phys_b))
            end

            connection_group_id = findfirst(item -> item == (phys_a, phys_b), unique_celltype_pairs)

            if !((phys_a, phys_b) in unique_celltype_pairs)
                push!(connection_groups, ConnectionGroup(
                    phys_a, phys_b,
                    flux_function!,
                    [(idx_a, Tuple{Int,Int}[]) for idx_a in 1:n_cells]
                )
                )
                push!(connection_groups[connection_group_id].connections[idx_a][2], (
                    (idx_b, facet_idx)
                )
                )
            elseif (phys_a, phys_b) in unique_celltype_pairs && isempty([connection_group_id].connections[idx_a])
                println("I'm unsure if this actually ever happens, check to line 221 in sim config if it does")
                push!(connection_groups[connection_group_id].connections[idx_a], (
                    (idx_a, Vector{Tuple{Int,Int}}((idx_b, facet_idx)))
                )
                )
            elseif !isempty([connection_group_id].connections[idx_a])
                push!(connection_groups[connection_group_id].connections[idx_a][2], (
                    (idx_b, facet_idx)
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
        connection_groups, phys, cell_phys_id_map,
        regions_phys_func_cells,
        controllers_controller_id_controller_func_monitored_affected,
        u_merged_axes,
    )
end