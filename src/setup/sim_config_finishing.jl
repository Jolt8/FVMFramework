
struct RegionGroup{P <: NamedTuple, F <: Function}
    name::String
    properties::P
    cache_syms_and_units::NamedTuple
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
    p_vec::Vector{Number}
    p_axes::NamedTuple
    merged_properties::NamedTuple
    du_diff_cache_vec::DiffCache
    u_diff_cache_vec::DiffCache
    du_proto_axes::NamedTuple
    u_proto_axes::NamedTuple
    du_cache_axes::NamedTuple
    u_cache_axes::NamedTuple
end

function finish_fvm_config(config, connection_map_function, special_caches; check_units::Bool)
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
        push!(region_groups, RegionGroup(region.name, region.properties, region.cache_syms_and_units, region.region_function, region.region_cells))

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
                #after using this framework for a while, I've never actually seen this trigger
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
                du_proto_nt[field][sub_field] .= 0.0 * unit(du_proto_nt[field][sub_field][1])
            end
        else
            du_proto_nt[field] .= 0.0 * unit(du_proto_nt[field][1])
        end
    end
    #du_proto_nt = (; du_proto...)

    u_proto_nt = deepcopy(u_merged)
    #u_proto_nt = (; u_proto...)

    merged_properties = merge_region_properties(config)

    merged_caches = merge_region_caches(config, special_caches, merged_properties)

    du_unitful_cache_nt = deepcopy(merged_caches)
    #du_cache_nt = (; du_cache...)

    u_unitful_cache_nt = deepcopy(merged_caches)
    #u_cache_nt = (; u_cache...)

    du0_vec_units = Vector(ComponentArray(; du_proto_nt...))
    u0_vec_units = Vector(ComponentArray(; u_proto_nt...))

    du0_vec = ustrip.(upreferred.(deepcopy(du0_vec_units)))
    u0_vec = ustrip.(upreferred.(deepcopy(u0_vec_units)))

    N::Int = ForwardDiff.pickchunksize(length(u0_vec))

    du_unitful_cache_vec = Vector(ComponentArray(; du_unitful_cache_nt...))
    u_unitful_cache_vec = Vector(ComponentArray(; u_unitful_cache_nt...))

    du_diff_cache_vec = DiffCache(ustrip.(upreferred.(deepcopy(du_unitful_cache_vec))), N)
    u_diff_cache_vec = DiffCache(ustrip.(upreferred.(deepcopy(u_unitful_cache_vec))), N)
    
    du_proto_axes = create_axes(du_proto_nt, n_cells)
    u_proto_axes = create_axes(u_proto_nt, n_cells)
    du_cache_axes = create_axes(du_unitful_cache_nt, n_cells)
    u_cache_axes = create_axes(u_unitful_cache_nt, n_cells)

    optimized_parameters_keys = keys(config.optimized_parameters)
    optimized_parameters_1_element_vectors = [[config.optimized_parameters[field]] for field in keys(config.optimized_parameters)]
    #we have to make it into a 1 element vector because create_axes requires a vector for each name

    optimized_parameter_nt = (; zip(optimized_parameters_keys, optimized_parameters_1_element_vectors)...)
    p_vec = Vector(ComponentVector(; optimized_parameter_nt...))

    p_axes = create_axes(optimized_parameter_nt, n_cells)


    #we have to split up properties to strip it of units
    properties_vec_units = Vector(ComponentArray(; merged_properties...))
    properties_vec = ustrip.(upreferred.(deepcopy(properties_vec_units)))

    properties_axes = create_axes(merged_properties, n_cells)

    if check_units == true
        system = FVMSystem(
            connection_groups, controller_groups, patch_groups, region_groups,
            p_vec, p_axes,
            create_views_inline(properties_vec_units, properties_axes),
            du_diff_cache_vec, u_diff_cache_vec,
            du_proto_axes, u_proto_axes,
            du_cache_axes, u_cache_axes
        )
        du_nt_units, u_nt_units = run_and_check_units(du0_vec_units, u0_vec_units, config.geo, system, du_unitful_cache_vec, u_unitful_cache_vec)
        return du_nt_units, u_nt_units, 0, 0
    end

    system = FVMSystem(
        connection_groups, controller_groups, patch_groups, region_groups,
        p_vec, p_axes,
        create_views_inline(properties_vec, properties_axes),
        du_diff_cache_vec, u_diff_cache_vec,
        du_proto_axes, u_proto_axes,
        du_cache_axes, u_cache_axes
    )

    return du0_vec, u0_vec, config.geo, system
end