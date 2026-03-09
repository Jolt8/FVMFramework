mutable struct RegionSetupInfo{P <: AbstractPhysics} #this must be defined before SimulationConfigInfo
    name::String
    type::P
    initial_conditions::NamedTuple
    properties::NamedTuple
    cache_syms::Vector{Symbol}
    region_function::Function
    region_cells::Vector{Int}
end

mutable struct PatchSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    properties::NamedTuple
    patch_function::Function
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}}
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
    patches::Vector{PatchSetupInfo}
    controllers::Vector{ControllerSetupInfo}
    optimized_parameters::Dict{Symbol, Any}
    u_proto::NamedTuple
end

function create_fvm_config(grid, u_proto)
    geo = build_fvm_geo_into_struct(grid)

    return SimulationConfigInfo(
        grid,
        geo,
        RegionSetupInfo[],
        PatchSetupInfo[],
        ControllerSetupInfo[],
        Dict{Symbol, Any}(), 
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
#=
This was a scrapped pattern where you put optimized_ in front of properties that you wanted to optimized
I decided to go with the explicit list of provided optimized paramters just to prevent typos from messing everything up


function check_if_optimized!(properties, field, config)
    if startswith(String(field), "optimized_") #these will never be named tuples
        param_name = Symbol(String(field)[(length("optimized_")+1):end])
        push!(config.optimized_parameters, OptimizedParameterTracker(param_name, properties[field]))
    end
end

function scan_properties_for_optimized!(properties, config)
    for field in propertynames(properties)
        if properties[field] isa NamedTuple
            scan_properties_for_optimized!(properties[field], config)
        else
            check_if_optimized!(properties, field, config)
        end
    end
end
=#

function add_region!(
    config, name;
    type,
    initial_conditions,
    properties,
    optimized_syms,
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

    for field in optimized_syms
        config.optimized_parameters[field] = properties[field]
    end

    append!(cache_syms, optimized_syms) #this is so that the optimized parameters are initialized as a DiffCache

    push!(config.regions, RegionSetupInfo(name, type, initial_conditions, properties, cache_syms, region_function, region_cells))
end

#this could probably also be handled by dynamic dispatch for facets, but it helps the user know a different routine is happening
#this needs to be updated, it currently doesn't work
function get_neighboring_cell_and_face_idx_from_face_idx(cell_id, face_idx, top)
    neighbor_info = top.face_face_neighbor[cell_id, face_idx]

    if !isempty(neighbor_info)
        face_idx_neighbor = first(neighbor_info)
        neighbor_id = face_idx_neighbor.idx[1]
        neighbor_face_idx = face_idx_neighbor.idx[2]
    else
        return nothing, nothing
    end

    return neighbor_id, neighbor_face_idx
end

function add_patch!(
    config, name;
    properties,
    optimized_syms,
    patch_function,
)
    cell_ids_and_face_idxs = [cell_id_facet_idx.idx for cell_id_facet_idx in keys(config.grid.facetsets[name].dict)]
    #[(cell_id, face_idx), (cell_id, face_idx), ...]

    n_cells = length(config.grid.cells)
    
    cell_neighbors = [(cell_id, Vector{Tuple{Int, Int}}()) for cell_id in 1:n_cells]

    top = ExclusiveTopology(config.grid)

    for (cell_id, face_idx) in cell_ids_and_face_idxs
        neighbor_id, neighbor_face_idx = get_neighboring_cell_and_face_idx_from_face_idx(cell_id, face_idx, top)
        if !isnothing(neighbor_id)
            push!(cell_neighbors[cell_id][2], (neighbor_id, face_idx))
            #push!(cell_neighbors[neighbor_id][2], (cell_id, neighbor_face_idx)) #I don't think this is necessary
        end
    end

    
    filter!(conn -> !(isempty(conn[2])), cell_neighbors)
    #get rid of empty connections

    for field in optimized_syms
        config.optimized_parameters[field] = properties[field]
    end

    append!(config.regions[1].cache_syms, optimized_syms) #this is so that the optimized parameters are initialized as a DiffCache

    push!(config.patches, PatchSetupInfo(name, properties, patch_function, cell_neighbors))
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
function _build_blank_dict!(current_dict, properties, n_cells)
    for (property_name, value) in pairs(properties)
        if value isa NamedTuple
            if !haskey(current_dict, property_name)
                current_dict[property_name] = Dict{Symbol,Any}()
            end
            _build_blank_dict!(current_dict[property_name], value, n_cells)
        else
            if !haskey(current_dict, property_name)

                current_dict[property_name] = zeros(Float64, n_cells)
            end
        end
    end
end

#... and this were AI generated because dealing with immutable NamedTuples with this level of nesting is fucked
#I do understand the code, so I don't think it's that big of a deal
function _dict_to_namedtuple(d::Dict)
    # Recursively convert Dicts to NamedTuples
    named_tuple_pairs = Pair{Symbol,Any}[]
    for (property_name, value) in pairs(d)
        if value isa Dict
            push!(named_tuple_pairs, property_name => _dict_to_namedtuple(value))
        else
            push!(named_tuple_pairs, property_name => value)
        end
    end
    return (; named_tuple_pairs...)
end

function _drill_down_and_fill_properties!(merged_properties, properties, cells)
    for (property_name, value) in pairs(properties)
        if merged_properties[property_name] isa NamedTuple
            _drill_down_and_fill_properties!(merged_properties[property_name], properties[property_name], cells)
        else
            for cell_id in cells
                if merged_properties[property_name] isa AbstractArray
                    merged_properties[property_name][cell_id] = properties[property_name]
                elseif merged_properties[property_name] isa Number
                    merged_properties[property_name] = properties[property_name]
                end
            end
        end
    end
end

function _drill_down_and_fill_patch_properties!(merged_properties, properties, cells)
    for (property_name, value) in pairs(properties)
        if merged_properties[property_name] isa NamedTuple
            _drill_down_and_fill_patch_properties!(merged_properties[property_name], properties[property_name], cells)
        else
            merged_properties[property_name][1] = properties[property_name]
        end
    end
end

function merge_region_properties(config)
    n_cells = length(config.geo.cell_volumes)

    prop_dict = Dict{Symbol,Any}()

    for region in config.regions
        _build_blank_dict!(prop_dict, region.properties, n_cells)
    end

    for patch in config.patches
        _build_blank_dict!(prop_dict, patch.properties, n_cells)
    end

    merged_properties = _dict_to_namedtuple(prop_dict)

    for region in config.regions
        _drill_down_and_fill_properties!(merged_properties, region.properties, region.region_cells)
    end

    for patch in config.patches
        _drill_down_and_fill_patch_properties!(merged_properties, patch.properties, patch.cell_neighbors)
    end

    return merged_properties
end


function _drill_down_and_fill_caches!(merged_caches, region_cache_syms, special_caches, merged_properties)
    for property_name in region_cache_syms
        if property_name in keys(special_caches)

        elseif merged_caches[property_name] isa NamedTuple
            _drill_down_and_fill_caches!(merged_caches[property_name], region_cache_syms, (_ = 0.0,), merged_properties)
        elseif merged_caches[property_name] isa AbstractArray
            #since we use the properties as the initial value for the cache, we might not even need properties and could just preemptively merge them
            if hasproperty(merged_properties, property_name) && merged_properties[property_name] isa AbstractArray
                merged_caches[property_name] .= merged_properties[property_name]
            else
                merged_caches[property_name] .= 0.0
            end
        else
            error("The cache $property_name was not handled")
        end
    end
end

function merge_region_caches(config, special_caches, merged_properties)
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
        _drill_down_and_fill_caches!(merged_caches, region.cache_syms, special_caches, merged_properties)
    end

    return merged_caches
end


#TODO: allow for face boundary conditions, boundary conditions applied on faces only apply to the faces right now

#we should probably find a way to make the creation of these structs automatic 
struct RegionGroup{P <: NamedTuple, F <: Function}
    name::String
    properties::P
    cache_syms::Vector{Symbol}
    region_function!::F
    region_cells::Vector{Int}
end

struct PatchGroup{P <: NamedTuple, F <: Function}
    name::String
    properties::P
    patch_function!::F
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}}
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
    controller_groups::Vector{ControllerGroup}
    patch_groups::Vector{PatchGroup}
    region_groups::Vector{RegionGroup}
    p_vec::Vector{Float64}
    p_axes::NamedTuple
    merged_properties::NamedTuple
    du_diff_cache_vec::DiffCache
    u_diff_cache_vec::DiffCache
    du_proto_axes::NamedTuple
    u_proto_axes::NamedTuple
    du_cache_axes::NamedTuple
    u_cache_axes::NamedTuple
end

function finish_fvm_config(config, connection_map_function, special_caches)
    n_cells = length(config.geo.cell_volumes)

    #in order from when they're executed in the solver 
    #1. Connections
    connection_groups = ConnectionGroup[]
    
    #2. Controllers
    controller_groups = ControllerGroup[]

    #3. Patches
    patch_groups = PatchGroup[]

    #4. Regions
    region_groups = RegionGroup[]

    #Regions
    cell_region_types_map = Vector{AbstractPhysics}(undef, n_cells)
    #although this could be a part of RegionGroup, I'd rather not contaminate it with information not required in the simulation

    for region in config.regions
        push!(region_groups, RegionGroup(region.name, region.properties, region.cache_syms, region.region_function, region.region_cells))

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

    #Patches
    for patch in config.patches
        push!(patch_groups, PatchGroup(patch.name, patch.properties, patch.patch_function, patch.cell_neighbors))
    end

    #Controllers
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

    #Connections
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

    merged_properties = merge_region_properties(config)

    merged_caches = merge_region_caches(config, special_caches, merged_properties)

    du_cache_nt = deepcopy(merged_caches)
    #du_cache_nt = (; du_cache...)

    u_cache_nt = deepcopy(merged_caches)
    #u_cache_nt = (; u_cache...)

    du0_vec = Vector(ComponentArray(; du_proto_nt...))
    u0_vec = Vector(ComponentArray(; u_proto_nt...))
    du_cache_vec = Vector(ComponentArray(; du_cache_nt...))
    u_cache_vec = Vector(ComponentArray(; u_cache_nt...))

    N::Int = ForwardDiff.pickchunksize(length(u0_vec))

    du_diff_cache_vec = DiffCache(du_cache_vec, N)
    u_diff_cache_vec = DiffCache(u_cache_vec, N)
    
    du_proto_axes = create_axes(du_proto_nt, n_cells)
    u_proto_axes = create_axes(u_proto_nt, n_cells)
    du_cache_axes = create_axes(du_cache_nt, n_cells)
    u_cache_axes = create_axes(u_cache_nt, n_cells)

    optimized_parameters_keys = keys(config.optimized_parameters)
    optimized_parameters_1_element_vectors = [[config.optimized_parameters[field]] for field in keys(config.optimized_parameters)]
    #we have to make it into a 1 element vector because create_axes requires a vector for each name

    optimized_parameter_nt = (; zip(optimized_parameters_keys, optimized_parameters_1_element_vectors)...)
    p_vec = Vector(ComponentVector(; optimized_parameter_nt...))

    p_axes = create_axes(optimized_parameter_nt, n_cells)

    return du0_vec, u0_vec, config.geo, FVMSystem(
        connection_groups, controller_groups, patch_groups, region_groups,
        p_vec, p_axes,
        merged_properties,
        du_diff_cache_vec, u_diff_cache_vec,
        du_proto_axes, u_proto_axes,
        du_cache_axes, u_cache_axes
    )
end

