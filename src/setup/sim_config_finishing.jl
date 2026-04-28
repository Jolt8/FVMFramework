
struct RegionGroup{P <: ComponentVector, F <: Function}
    name::String
    properties::P
    region_function!::F
    region_cells::Vector{Int}
end

struct PatchGroup{P <: ComponentVector, F <: Function}
    name::String
    properties::P
    patch_function!::F
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}}
end

struct ControllerGroup{T <: ComponentVector, F <: Function}
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
    du_virtual_axes::NamedTuple
    u_virtual_axes::NamedTuple
    du_diff_cache::DiffCache
    u_diff_cache::DiffCache
    properties_vec::Vector{Float64}
    properties_axes::Tuple
    p_vec::Vector{Float64}
    p_axes::Tuple
end

function finish_fvm_config(config, connection_map_function; check_units::Bool)
    n_cells = length(config.geo.cell_volumes)

    connection_groups = ConnectionGroup[]
    controller_groups = ControllerGroup[]
    patch_groups = PatchGroup[]
    region_groups = RegionGroup[]

    #Regions
    cell_region_phys_map = Vector{AbstractPhysics}(undef, n_cells)
    #although this could be a part of RegionGroup, I'd rather not contaminate it with information not required in the simulation

    for region in config.regions
        push!(region_groups, RegionGroup(region.name, region.properties, region.region_function, region.region_cells))

        for cell_id in region.region_cells
            cell_region_phys_map[cell_id] = region.type
        end
    end

    cell_regions_map = Vector{RegionGroup}(undef, n_cells)

    for region in region_groups
        for cell_id in region.region_cells
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
    #we specifically use strings here because checking if (region_a, region_b) == (region_b, region_a) was fragile

    for (idx_a, idx_a_neighbors) in config.geo.cell_neighbors
        for (idx_b, face_idx) in idx_a_neighbors
            if idx_b <= 0
                continue
            end

            region_a = cell_regions_map[idx_a]
            region_b = cell_regions_map[idx_b]

            flux_function! = connection_map_function(cell_region_phys_map[idx_a], cell_region_phys_map[idx_b])

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
            #=elseif (region_a.name, region_b.name) in unique_region_connection_pairs && isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                println("I'm unsure if this actually ever happens, check to line 237 in sim config if it does")
                #after using this framework for a while, I've never actually seen this trigger
                connection_group_id = findfirst(item -> item == (region_a.name, region_b.name), unique_region_connection_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a], (
                    (idx_a, Vector{Tuple{Int, Int}}((idx_b, face_idx)))
                )
                )=#
            elseif !isempty(connection_groups[connection_group_id].cell_neighbors[idx_a])
                connection_group_id = findfirst(item -> item == (region_a.name, region_b.name), unique_region_connection_pairs)
                push!(connection_groups[connection_group_id].cell_neighbors[idx_a][2], (
                    (idx_b, face_idx)
                )
                )
            end
        end
    end

    #remove empty tuples from cell_neighbors
    for CG in connection_groups
        filter!(conn -> !(isempty(conn[2])), CG.cell_neighbors)
    end

    merged_properties = merge_region_properties(config)
    properties_axes = getaxes(merged_properties)

    merged_caches = merge_region_caches(config, merged_properties)

    merged_caches = merge_in_second_order_caches(config, merged_caches, merged_properties)

    state_axes = getaxes(config.u_proto)
    
    du0_vec_units = Vector(deepcopy(upreferred.(config.u_proto)))
    du0_vec_units .*= 0.0
    u0_vec_units = Vector(deepcopy(upreferred.(config.u_proto)))

    du0_vec = ustrip.(upreferred.(deepcopy(du0_vec_units)))
    du0_vec .= 0.0
    u0_vec = ustrip.(upreferred.(deepcopy(u0_vec_units)))

    N::Int = ForwardDiff.pickchunksize(length(u0_vec))
    
    du_unitful_cache_vec = Vector(ComponentArray(; deepcopy(merged_caches)...))
    u_unitful_cache_vec = Vector(ComponentArray(; deepcopy(merged_caches)...))

    du_diff_cache = DiffCache(ustrip.(upreferred.(deepcopy(du_unitful_cache_vec))), N)
    u_diff_cache = DiffCache(ustrip.(upreferred.(deepcopy(u_unitful_cache_vec))), N)

    pre_lengthened_optimized_parameters_dict = Dict{Symbol, Any}()

    for name in propertynames(config.optimized_parameters)
        pre_lengthened_optimized_parameters_dict[name] = getproperty(config.optimized_parameters, name)
    end

    lengthened_optimized_parameters = ComponentVector(_dict_to_namedtuple(pre_lengthened_optimized_parameters_dict))

    p_vec_units = Vector(ComponentVector(lengthened_optimized_parameters))
    p_vec = ustrip.(upreferred.(Vector(ComponentVector(lengthened_optimized_parameters))))
    p_axes = getaxes(ComponentVector(lengthened_optimized_parameters))

    #we have to split up properties to strip it of units
    properties_vec_units = Vector(merged_properties)
    properties_vec = ustrip.(upreferred.(deepcopy(properties_vec_units)))

    #virtual_merge_axes takes in a tuple of ComponentArrays
    du_virtual_axes = virtual_merge_axes((ComponentVector(config.u_proto), ComponentVector(merged_caches)))
    if ComponentVector(config.optimized_parameters) == Union{}[] #if the user doesn't provide any optimized parameters, we ignore them
        u_virtual_axes = virtual_merge_axes((ComponentVector(config.u_proto), ComponentVector(merged_caches), ComponentVector(merged_properties)))
    else
        u_virtual_axes = virtual_merge_axes((ComponentVector(config.u_proto), ComponentVector(merged_caches), ComponentVector(merged_properties), ComponentVector(lengthened_optimized_parameters)))
    end

    if check_units == true
        system = FVMSystem(
            connection_groups, controller_groups, patch_groups, region_groups,
            du_virtual_axes, u_virtual_axes,
            du_diff_cache, u_diff_cache,
            properties_vec, properties_axes,
            p_vec, p_axes
        )
        du_units, u_units = run_and_check_units(du0_vec_units, u0_vec_units, config.geo, system, du_unitful_cache_vec, u_unitful_cache_vec, properties_vec_units, p_vec_units)
        return du_units, u_units, state_axes, 0, 0
    end

    system = FVMSystem(
        connection_groups, controller_groups, patch_groups, region_groups,
        du_virtual_axes, u_virtual_axes,
        du_diff_cache, u_diff_cache,
        properties_vec, properties_axes,
        p_vec, p_axes
    )

    return du0_vec, u0_vec, state_axes, config.geo, system
end