mutable struct RegionSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    initial_conditions::ComponentArray
    region_physics::AbstractPhysics
    region_function::Function
    region_cells::Vector{Int}
end

struct SimulationConfigInfo
    grid::Ferrite.Grid
    regions::Vector{RegionSetupInfo}
    u_proto::ComponentArray
    u_axes::Axis
    geo::FVMGeometry
end

function create_fvm_config(grid, u_proto)
    u_axes = getaxes(u_proto)[1]

    geo = build_fvm_geo_into_struct(grid)

    return SimulationConfigInfo(
        grid,
        RegionSetupInfo[],
        u_proto,
        u_axes,
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
    region_physics,
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

    push!(config.regions, RegionSetupInfo(name, initial_conditions, region_physics, region_function, region_cells))
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

#TODO: allow for face boundary conditions, boundary conditions applied on faces only apply to the faces right now

#we should probably find a way to make the creation of these structs automatic 
abstract type AbstractConnectionGroup end
finalize_connection_groups(groups::AbstractConnectionGroup) = groups

struct FVMSystem
    connection_groups::AbstractConnectionGroup
    phys::Vector{AbstractPhysics}
    cell_phys_id_map::Vector{Int}
    regions_phys_func_cells::Vector{Tuple{AbstractPhysics, Function, Vector{Int}}}
    u_axes::Axis
end

function finish_fvm_config(config, connection_catagorizer!, connection_groups::AbstractConnectionGroup)
    n_cells = length(config.geo.cell_volumes)

    phys = AbstractPhysics[]
    cell_phys_id_map = zeros(Int, n_cells)
    #the only reason this exists is because we still need to access the physics of the cell in the connections loop
    #we're doubling up on physics, but it makes everything cleaner later

    regions_phys_func_cells = Tuple{AbstractPhysics, Function, Vector{Int}}[]

    for (i, region) in enumerate(config.regions)
        push!(regions_phys_func_cells, (region.region_physics, region.region_function, region.region_cells))
        push!(phys, region.region_physics)

        @batch for cell_id in region.region_cells
            cell_phys_id_map[cell_id] = i
        end
    end

    @batch for (idx_a, idx_a_neighbors) in config.geo.cell_neighbors
        current_connection = (idx_a, Int[]) #(cell_id, Vector{neighbor_ids})
        for (idx_b, face_idx) in idx_a_neighbors
            if idx_b <= 0
                continue
            end
            phys_a = phys[cell_phys_id_map[idx_a]]
            phys_b = phys[cell_phys_id_map[idx_b]]

            push!(current_connection[2], idx_b)
            connection_catagorizer!(connection_groups, current_connection, idx_a, idx_b, face_idx, typeof(phys_a), typeof(phys_b))
        end
    end

    u0 = Vector(config.u_proto)
    du0 = u0 .* 0.0

    return du0, u0, config.geo, FVMSystem(
        finalize_connection_groups(connection_groups), phys, cell_phys_id_map,
        regions_phys_func_cells,
        config.u_axes,
    )
end