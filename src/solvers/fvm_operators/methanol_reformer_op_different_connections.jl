function methanol_reformer_f_test!(
    du_vec, u_vec, p, t, 

    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    connection_groups, controller_groups, region_groups, patch_groups,
    
    properties,

    du_diff_cache_vec, u_diff_cache_vec,
    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes
)
    du_cache_vec = get_tmp(du_diff_cache_vec, u_vec)
    u_cache_vec = get_tmp(u_diff_cache_vec, u_vec)
    
    #do you know what is fucking crazy?
    #the only reason this worked in the past was because I zeroed out the caches
    #otherwise, it just returns #undef for everything
    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = create_views_inline(du_cache_vec, du_cache_axes)
    u_cache_nt = create_views_inline(u_cache_vec, u_cache_axes)

    du_vec .= 0.0
    du_nt = create_views_inline(du_vec, du_proto_axes)
    u_nt = create_views_inline(u_vec, u_proto_axes)

    du = (; du_nt..., du_cache_nt...)
    #u = (; u_nt..., u_cache_nt..., properties...)
    u = (; properties..., u_nt..., u_cache_nt...)
    #remember, the right-most fields overwrite other fields with the same name

    u.rho .= (u_cache_nt.rho .+ properties.rho)
    #u = (; properties..., u_nt..., u_cache_nt...)
    #while the above would seem to be better, for some reason it doesn't work
    #the issue is that because we have to set caches to .= 0.0 to prevent them from being undefined, the properties inside the cache vec are overwritten

    #Check heat_transfer_minimal_allocs for more info on comparing fetching NamedTuples
    #I FOUND THE PROBLEM!
    #This pattern doesn't work for mass fractions or other variables that have nested structures such as mass fractions
    #that's also where the SubArrays were coming from
    #thank fucking god!

    #=
    du_cache_ca = ComponentArray(du_cache_vec, du_cache_axes)
    du_cache_nt = NamedTuple{propertynames(du_cache_ca)}(Tuple(getproperty(du_cache_ca, p) for p in propertynames(du_cache_ca)))

    u_cache_ca = ComponentArray(u_cache_vec, u_cache_axes)
    u_cache_nt = NamedTuple{propertynames(u_cache_ca)}(Tuple(getproperty(u_cache_ca, p) for p in propertynames(u_cache_ca)))

    du_vec .= 0.0
    du_ca = ComponentArray(du_vec, du_proto_axes)
    du_nt = NamedTuple{propertynames(du_ca)}(Tuple(getproperty(du_ca, p) for p in propertynames(du_ca)))

    u_ca = ComponentArray(u_vec, u_proto_axes)
    u_nt = NamedTuple{propertynames(u_ca)}(Tuple(getproperty(u_ca, p) for p in propertynames(u_ca)))

    du = (; du_nt..., du_cache_nt...)
    u = (; u_nt..., u_cache_nt..., properties...)
    =#

    #Connection Loops
    for conn in connection_groups
        solve_connection_group!(
            du, u,
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        )
    end

    #=
    #Controller Loops
    for cont in controller_groups
        solve_controller_group!(
            du, u, cont.id, cont.controller, cont.controller_function!, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
    end
    =#

    for patch in patch_groups
        solve_patch_group!(
            du, u,
            patch.patch_function!, patch.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            cell_volumes
        )
    end

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops
    for reg in region_groups
        solve_region_group!(
            du, u,
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end

    #=
    for reg in region_groups
        if reg.name == "surrounding_tissue"
            total_methylene_blue_dm_dt = 0.0
            total_water_dm_dt = 0.0
            total_reservoir_mass = 0.0
            total_reservoir_volume = 0.0

            for cell_id in reg.region_cells
                total_methylene_blue_dm_dt += du.mass_fractions.methylene_blue[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]
                total_water_dm_dt += du.mass_fractions.water[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]

                total_reservoir_mass += u.rho[cell_id] * cell_volumes[cell_id]
                total_reservoir_volume += cell_volumes[cell_id]
            end

            well_mixed_methylene_blue_dt = total_methylene_blue_dm_dt / total_reservoir_mass
            well_mixed_water_dt = total_water_dm_dt / total_reservoir_mass
            
            for cell_id in reg.region_cells
                du.mass_fractions.methylene_blue[cell_id] = well_mixed_methylene_blue_dt
                du.mass_fractions.water[cell_id] = well_mixed_water_dt
            end
        end
        if reg.name == "implant_interior"
            total_methylene_blue_dm_dt = 0.0
            total_water_dm_dt = 0.0
            total_reservoir_mass = 0.0
            total_reservoir_volume = 0.0

            for cell_id in reg.region_cells
                total_methylene_blue_dm_dt += du.mass_fractions.methylene_blue[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]
                total_water_dm_dt += du.mass_fractions.water[cell_id] * u.rho[cell_id] * cell_volumes[cell_id]

                total_reservoir_mass += u.rho[cell_id] * cell_volumes[cell_id]
                total_reservoir_volume += cell_volumes[cell_id]
            end

            well_mixed_methylene_blue_dt = total_methylene_blue_dm_dt / total_reservoir_mass
            well_mixed_water_dt = total_water_dm_dt / total_reservoir_mass
            
            for cell_id in reg.region_cells
                du.mass_fractions.methylene_blue[cell_id] = well_mixed_methylene_blue_dt
                du.mass_fractions.water[cell_id] = well_mixed_water_dt
            end
        end
    end
    =#


    #du_vec[1:length(du.mass_fractions.methylene_blue)] = Vector(du.mass_fractions.methylene_blue)
    #du_vec[length(du.mass_fractions.methylene_blue)+1:length(du.mass_fractions.methylene_blue)+length(du.mass_fractions.water)] = Vector(du.mass_fractions.water)

    #=
    pressure = []
    mass_fractions = []
    temp = []
    heat = []
    rho = []
    mw_avg = []

    region_of_interest = region_groups[3]

    for cell_id in region_of_interest.region_cells
        if isnan(du.pressure[cell_id])
            push!(pressure, du.pressure[cell_id])
        end
        if isnan(du.mass_fractions.methanol[cell_id])
            push!(mass_fractions, du.mass_fractions.methanol[cell_id])
        end
        if isnan(du.temp[cell_id])
            push!(temp, du.temp[cell_id])
            push!(heat, du.heat[cell_id])
        end
        if isnan(du.rho[cell_id]) || u.rho[cell_id] == 0.0
            push!(rho, du.rho[cell_id])
        end
        if u.mw_avg[cell_id] == 0.0
            push!(mw_avg, du.mw_avg[cell_id])
        end
    end

    println(region_of_interest.name)
    println("pressure: ", pressure)
    println("mass_fractions: ", mass_fractions)
    println("temp: ", temp)
    println("heat: ", heat)
    println("rho: ", rho)
    println("mw_avg: ", mw_avg)



    #=
    for region in region_groups
        region_encountered_nan_symbols = Symbol[]
        for key in keys(u)
            println(length(u[key]))
            if !(key in region_encountered_nan_symbols)
                for cell_id in region.region_cells
                    if u[key] isa SubArray
                        if isnan(u[key][1])
                            push!(region_encountered_nan_symbols, key)
                            break
                        end
                    elseif length(u[key]) >= length(region.region_cells)
                        if length(u[key]) == length(region.region_cells)
                            if isnan(u[key][cell_id])
                                push!(region_encountered_nan_symbols, key)
                                break
                            end
                        else
                            if isnan(u[key][1][1])
                                push!(region_encountered_nan_symbols, key)
                                break
                            end
                        end
                    end
                end
            end
        end

        println("encountered NaN in: $(region.name) : $(region_encountered_nan_symbols)")
    end

    =#


    for region in region_groups
        for cell_id in region.region_cells
            #println(u.mw_avg[cell_id])
            if isnan(u.rho[cell_id])
                println(region.name)
                println("mw_avg: ", u.mw_avg[cell_id])
                println("rho: ", u.rho[cell_id])
                println("pressure: ", u.pressure[cell_id])
                println("temp: ", u.temp[cell_id])
            end
        end
    end

    for reg in region_groups
        debug_region!(du_nt, u_nt, reg)
    end
    =#
end

#VERY IMPORTANT!!!!
#= For future reference when writing to properties using u[field]:
    u_cv.temp[cell_id] = val
        works!
    u_cv[field][cell_id] = val
        fails :(
    view(u_cv, field)[cell_id] = val
        works!
    getproperty(u_cv, field)[cell_id] = val
        works!
=#

#OH SHIT!, mass_fractions[species_name] only works for ComponentVectors it's a scalar, not a vector
#so when mass_fractions[species_name] is a vector of n_cells we have to use view(mass_fractions, species_name)[cell_id]

#instead of updating rho in each flux and internal physics, we could do an initial property retrieval loop
#=
@batch for cell_id in eachindex(cell_volumes)
    mw_avg!(u, cell_id, molecular_weights, mw_avg_cache)
    rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)
end
=#