mutable struct RegionSetupInfo{P <: AbstractPhysics} #this must be defined before SimulationConfigInfo
    name::String
    type::P
    initial_conditions::NamedTuple
    properties::NamedTuple
    state_syms::Vector{Symbol}
    cache_syms::Vector{Symbol}
    region_function::Function
    region_cells::Vector{Int}
end

mutable struct ControllerSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    controller::NamedTuple
    monitored_cellset::String
    affected_cellset::String
    controller_function::Function
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

struct SimulationConfigInfo
    grid::Ferrite.Grid
    geo::FVMGeometry
    regions::Vector{RegionSetupInfo}
    controllers::Vector{ControllerSetupInfo}
    u_proto::NamedTuple
end

function create_fvm_config(grid, u_proto)
    geo = build_fvm_geo_into_struct(grid)

    return SimulationConfigInfo(
        grid,
        geo,
        RegionSetupInfo[],
        ControllerSetupInfo[],
        u_proto
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
#=
function _drill_down_and_fill_region_properties!(merged_properties, initial_conditions)
    for (property_name, value) in pairs(initial_conditions)
        if merged_properties[property_name] isa NamedTuple
            _drill_down_and_fill_region_properties!(merged_properties[property_name], initial_conditions)
        elseif merged_properties[property_name] isa AbstractArray || merged_properties[property_name] isa Number
            merged_properties[property_name] = initial_conditions[property_name]
        end
    end
end
=#

function add_region!(
    config, name;
    type,
    initial_conditions,
    properties,
    state_syms,
    cache_syms,
    region_function
)

    region_cells = collect(getcellset(config.grid, name))

    for cell_id in region_cells
        for field in propertynames(initial_conditions)
            var = config.u_proto[field]
            initial_condition = initial_conditions[field]
            if var isa NamedTuple 
                for sub_name in propertynames(var)
                    var[sub_name][cell_id] = initial_condition[sub_name]
                end
            else
                var[cell_id] = initial_condition
            end
        end
    end

    push!(config.regions, RegionSetupInfo(name, type, initial_conditions, properties, state_syms, cache_syms, region_function, region_cells))
end

#this could probably also be handled by dynamic dispatch for facets, but it helps the user know a different routine is happening
#this needs to be updated, it currently doesn't work
function add_facet_region!(
    config, name;
    type,
    initial_conditions,
    region_physics,
    region_function,
)

    region_cells = get_cell_ids_in_facet_set(config.grid, name)

    for cell_id in region_cells
        for field in keys(initial_conditions)
            var = config.u_proto[field][cell_id]
            initial_condition = initial_conditions[field]
            if var isa NamedTuple 
                for sub_name in keys(var)
                    var[sub_name][cell_id] = initial_condition[sub_name]
                end
            else
                var[cell_id] = initial_condition
            end
        end
    end

    push!(config.regions, RegionSetupInfo(name, type, initial_conditions, region_physics, state_syms, cache_syms, region_function, region_cells))
end

function add_controller!(
    config, name;
    controller,
    monitored_cellset,
    affected_cellset,
    controller_function
)
    monitored_cells = collect(getcellset(config.grid, monitored_cellset))
    affected_cells = collect(getcellset(config.grid, affected_cellset))

    push!(config.controllers, ControllerSetupInfo(name, controller, monitored_cellset, affected_cellset, controller_function, monitored_cells, affected_cells))
end

#this...
function _build_blank_dict!(current_dict, region_properties, n_cells)
    for (property_name, value) in pairs(region_properties)
        if value isa NamedTuple
            if !haskey(current_dict, property_name)
                current_dict[property_name] = Dict{Symbol, Any}()
            end
            _build_blank_dict!(current_dict[property_name], value, n_cells)
        else
            if !haskey(current_dict, property_name)
                current_dict[property_name] = zeros(typeof(value), n_cells) 
            end
        end
    end
end

#... and this were AI generated because dealing with immutable NamedTuples with this level of nesting is fucked
#I do understand the code, however, so I don't think it that big of a deal
function _dict_to_namedtuple(d::Dict)
    # Recursively convert Dicts to NamedTuples
    named_tuple_pairs = Pair{Symbol, Any}[]
    for (property_name, value) in pairs(d)
        if value isa Dict
            push!(named_tuple_pairs, property_name => _dict_to_namedtuple(value))
        else
            push!(named_tuple_pairs, property_name => value)
        end
    end
    return (; named_tuple_pairs...)
end

function _drill_down_and_fill_properties!(merged_properties, region_properties, region_cells)
    for (property_name, value) in pairs(region_properties)
        if merged_properties[property_name] isa NamedTuple
            _drill_down_and_fill_properties!(merged_properties[property_name], region_properties[property_name], region_cells)
        else
            for cell_id in region_cells
                if merged_properties[property_name] isa AbstractArray
                    merged_properties[property_name][cell_id] = region_properties[property_name]
                elseif merged_properties[property_name] isa Number
                    merged_properties[property_name] = region_properties[property_name]
                end
            end
        end
    end
end

function merge_region_properties(config)
    n_cells = length(config.geo.cell_volumes)

    prop_dict = Dict{Symbol, Any}()

    for region in config.regions
        _build_blank_dict!(prop_dict, region.properties, n_cells)
    end
    
    merged_properties = _dict_to_namedtuple(prop_dict)

    for region in config.regions
        _drill_down_and_fill_properties!(merged_properties, region.properties, region.region_cells)
    end

    return merged_properties
end


function _drill_down_and_fill_caches!(merged_caches, region_cache_syms, special_caches)
    for property_name in region_cache_syms
        if property_name in keys(special_caches)
            
        elseif merged_caches[property_name] isa NamedTuple
            _drill_down_and_fill_caches!(merged_caches[property_name], region_cache_syms, (_ = 0.0,)) #we'll never recurse for special caches
        elseif merged_caches[property_name] isa AbstractArray
            merged_caches[property_name] .= 0.0
        else
            error("The cache $property_name was not handled")
        end
    end
end

function merge_region_caches(config, special_caches)
    n_cells = length(config.geo.cell_volumes)

    cache_dict = Dict{Symbol, Any}()

    for region in config.regions
        for cache_name in region.cache_syms
            if !haskey(cache_dict, cache_name)
                if cache_name in keys(special_caches)
                    cache_dict[cache_name] = special_caches[cache_name]
                else
                    cache_dict[cache_name] = zeros(Float64, n_cells)
                end
            end
        end
    end
    
    merged_caches = _dict_to_namedtuple(cache_dict)

    for region in config.regions
        _drill_down_and_fill_caches!(merged_caches, region.cache_syms, special_caches)
    end

    return merged_caches
end


#TODO: allow for face boundary conditions, boundary conditions applied on faces only apply to the faces right now

#we should probably find a way to make the creation of these structs automatic 
struct RegionGroup{P <: NamedTuple, F <: Function}
    name::String
    properties::P
    state_syms::Vector{Symbol}
    cache_syms::Vector{Symbol}
    region_function!::F
    region_cells::Vector{Int}
end

struct ControllerGroup{T <: NamedTuple, F <: Function}
    name::String
    id::Int
    controller::T
    controller_function!::F
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

#=
struct CellNeighbors
    idx_a::Int
    idx_b_face_idx_vec::Vector{Tuple{Int, Int}}
end
=#

struct ConnectionGroup{RA <: RegionGroup, RB <: RegionGroup, F <: Function}
    region_a::RA
    region_b::RB
    flux_function!::F
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}} #we could probably use ::Vector{CellNeighbors} here
end

struct FVMSystem
    connection_groups::Vector{ConnectionGroup}
    region_groups::Vector{RegionGroup}
    controller_groups::Vector{ControllerGroup}
    merged_properties::NamedTuple
    du_diff_cache_vec::DiffCache
    u_diff_cache_vec::DiffCache
    du_proto_axes::Axis
    u_proto_axes::Axis
    du_cache_axes::Axis
    u_cache_axes::Axis
end

function finish_fvm_config(config, connection_map_function, special_caches)
    n_cells = length(config.geo.cell_volumes)

    region_groups = RegionGroup[]
    cell_region_types_map = Vector{AbstractPhysics}(undef, n_cells) 
    #although this could be a part of RegionGroup, I'd rather not contaminate it with information not required in the simulation

    for region in config.regions
        push!(region_groups, RegionGroup(region.name, region.properties, region.state_syms, region.cache_syms, region.region_function, region.region_cells))
        
        for cell_id in region.region_cells
            cell_region_types_map[cell_id] = region.type
        end
    end

    cell_properties_map = Vector{NamedTuple}(undef, n_cells)
    cell_regions_map = Vector{RegionGroup}(undef, n_cells)

    for region in region_groups
        for cell_id in region.region_cells
            cell_properties_map[cell_id] = region.properties
            cell_regions_map[cell_id] = region
        end
    end

    controller_groups = ControllerGroup[]

    for (controller_id, controller) in enumerate(config.controllers)
        push!(controller_groups, ControllerGroup(
            controller.name,
            controller_id,
            controller.controller,
            controller.controller_function,
            controller.monitored_cells, controller.affected_cells
        )
        )
    end

    connection_groups = ConnectionGroup[]

    unique_region_connection_pairs = Vector{Tuple{String, String}}()
    #we specifically use strings here because checking if (region_a, region_b) == (region_b, region_a) was fragile for some reason

    for (idx_a, idx_a_neighbors) in config.geo.cell_neighbors
        for (idx_b, face_idx) in idx_a_neighbors
            if idx_b <= 0
                continue
            end

            region_a = cell_regions_map[idx_a]
            region_b = cell_regions_map[idx_b]

            flux_function! = connection_map_function(cell_region_types_map[idx_a], cell_region_types_map[idx_b])

            connection_group_id = findfirst(item -> item == (region_a.name, region_b.name), unique_region_connection_pairs)

            #I feel bad for anyone who has to read this
            if !((region_a.name, region_b.name) in unique_region_connection_pairs)
                push!(unique_region_connection_pairs, (region_a.name, region_b.name))
                new_connection_group_id = findfirst(item -> item == (region_a.name, region_b.name), unique_region_connection_pairs)
                push!(connection_groups, ConnectionGroup(
                    region_a,
                    region_b,
                    flux_function!,
                    [(idx_a, Tuple{Int, Int}[]) for idx_a in 1:n_cells]
                )
                )
                push!(connection_groups[new_connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            elseif (region_a.name, region_b.name) in unique_region_connection_pairs && isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                println("I'm unsure if this actually ever happens, check to line 237 in sim config if it does")
                connection_group_id = findfirst(item -> item == (region_a.name, region_b.name), unique_region_connection_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a], (
                    (idx_a, Vector{Tuple{Int, Int}}((idx_b, face_idx)))
                )
                )
            elseif !isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                connection_group_id = findfirst(item -> item == (region_a.name, region_b.name), unique_region_connection_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            end
        end
    end

    #If you want one that organizes based on unique groups, check the Tracer branch in sim_config and you can find it commented out there

    #remove empty tuples from cell_neighbors
    for CG in connection_groups
        filter!(conn -> !(isempty(conn[2])), CG.cell_neighbors)
    end

    #this code is to merge the state variables such as velocity, temperature or pressure with non-state variables such as integral_error
    #integral error is indexed at its respective controller id
    n_controllers = length(config.controllers)
    u_non_state_proto = ()
    if n_controllers != 0.0
        u_non_state_proto = (
            integral_error = zeros(n_controllers),
        )
    end
    #unlike u_proto, u_non_state_proto is accessed at controller_id

    u_merged = (; config.u_proto..., u_non_state_proto...)

    #deepcopy is needed for copying NamedTuples

    du_proto_nt = deepcopy(u_merged)
    #it might be worth it to convert it back into a ComponentArray so we can just broadcast .= 0.0 on everything 
    for (field, value) in pairs(du_proto_nt)
        if du_proto_nt[field] isa NamedTuple
            for (sub_field, value) in pairs(du_proto_nt[field])
                du_proto_nt[field][sub_field] .= 0.0
            end
        else
            du_proto_nt[field] .= 0.0
        end
    end
    #du_proto_nt = (; du_proto...)

    u_proto_nt = deepcopy(u_merged)
    #u_proto_nt = (; u_proto...)

    merged_caches = merge_region_caches(config, special_caches)

    du_cache_nt = deepcopy(merged_caches)
    #du_cache_nt = (; du_cache...)

    u_cache_nt = deepcopy(merged_caches)
    #u_cache_nt = (; u_cache...)

    du0_vec = Vector(ComponentArray(; du_proto_nt...))
    u0_vec = Vector(ComponentArray(; u_proto_nt...))
    du_cache_vec = Vector(ComponentArray(; du_cache_nt...))
    u_cache_vec = Vector(ComponentArray(; u_cache_nt...))

    N::Int = ForwardDiff.pickchunksize(length(u0_vec))

    du_proto_axes = getaxes(ComponentArray(; du_proto_nt...))[1]
    u_proto_axes = getaxes(ComponentArray(; u_proto_nt...))[1]
    du_cache_axes = getaxes(ComponentArray(; du_cache_nt...))[1]
    u_cache_axes = getaxes(ComponentArray(; u_cache_nt...))[1]

    du_diff_cache_vec = DiffCache(du_cache_vec, N)
    u_diff_cache_vec = DiffCache(u_cache_vec, N)

    merged_properties = merge_region_properties(config)

    return du0_vec, u0_vec, config.geo, FVMSystem(
        connection_groups, region_groups, controller_groups, 
        merged_properties, 
        du_diff_cache_vec, u_diff_cache_vec,
        du_proto_axes, u_proto_axes,
        du_cache_axes, u_cache_axes
    )
end