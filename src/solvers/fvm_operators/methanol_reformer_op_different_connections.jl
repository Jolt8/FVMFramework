struct MethanolReformerPhysics <: AbstractFluidPhysics
    k::Float64
    cp::Float64
    mu::Float64
    permeability::Float64
    species_diffusion_coeffs::Vector{Float64}
    species_molecular_weights::Vector{Float64}
    chemical_reactions::Vector{AbstractReaction}
    cell_kg_cat_per_m3_for_each_reaction::Vector{Float64}
end

struct WallPhysics <: AbstractSolidPhysics
    k::Float64
    rho::Float64
    cp::Float64
end

struct MethanolReformerConnectionGroups <: AbstractConnectionGroup
    fluid_fluid::Vector{Tuple{Int, Vector{Tuple{Int,Int}}}}
    solid_solid::Vector{Tuple{Int, Vector{Tuple{Int,Int}}}}
    fluid_solid::Vector{Tuple{Int, Vector{Tuple{Int,Int}}}}
    solid_fluid::Vector{Tuple{Int, Vector{Tuple{Int,Int}}}}
end

struct MethanolReformerConnectionGroupsBuilder <: AbstractConnectionGroup
    fluid_fluid::Vector{Vector{Tuple{Int, Int}}}
    solid_solid::Vector{Vector{Tuple{Int, Int}}}
    fluid_solid::Vector{Vector{Tuple{Int, Int}}}
    solid_fluid::Vector{Vector{Tuple{Int, Int}}}
end

function methanol_reformer_init_conn_groups(grid)
    n_cells = length(grid.cells)
    return MethanolReformerConnectionGroupsBuilder(
        [Vector{Tuple{Int, Int}}() for _ in 1:n_cells],
        [Vector{Tuple{Int, Int}}() for _ in 1:n_cells],
        [Vector{Tuple{Int, Int}}() for _ in 1:n_cells],
        [Vector{Tuple{Int, Int}}() for _ in 1:n_cells]
    )
end

#this is to remove cells that have no connecitons in each connection group 
function finalize_connection_groups(builder::MethanolReformerConnectionGroupsBuilder)
    return MethanolReformerConnectionGroups(
        [(idx_a, builder.fluid_fluid[idx_a]) for idx_a in 1:length(builder.fluid_fluid) if !isempty(builder.fluid_fluid[idx_a])],
        [(idx_a, builder.solid_solid[idx_a]) for idx_a in 1:length(builder.solid_solid) if !isempty(builder.solid_solid[idx_a])],
        [(idx_a, builder.fluid_solid[idx_a]) for idx_a in 1:length(builder.fluid_solid) if !isempty(builder.fluid_solid[idx_a])],
        [(idx_a, builder.solid_fluid[idx_a]) for idx_a in 1:length(builder.solid_fluid) if !isempty(builder.solid_fluid[idx_a])]
    )
end

function connection_catagorizer!(builder::MethanolReformerConnectionGroupsBuilder, current_connection, idx_a, idx_b, face_idx, type_a, type_b)
    if type_a <: AbstractFluidPhysics && type_b <: AbstractFluidPhysics
        push!(builder.fluid_fluid[idx_a], (idx_b, face_idx))
    elseif type_a <: AbstractSolidPhysics && type_b <: AbstractSolidPhysics
        push!(builder.solid_solid[idx_a], (idx_b, face_idx))
    elseif (type_a <: AbstractFluidPhysics && type_b <: AbstractSolidPhysics)
        push!(builder.fluid_solid[idx_a], (idx_b, face_idx))
    elseif (type_a <: AbstractSolidPhysics && type_b <: AbstractFluidPhysics)
        push!(builder.solid_fluid[idx_a], (idx_b, face_idx))
    end 
end

function fluid_fluid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    phys_a, phys_b,
    rho_cache, mw_avg_cache,
    change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
)
    area = cell_neighbor_areas[idx_a][face_idx]
    dist = cell_neighbor_distances[idx_a][face_idx]
    norm = cell_neighbor_normals[idx_a][face_idx]

    rho_a = rho_ideal!(u, idx_a, rho_cache, mw_avg_cache)
    rho_b = rho_ideal!(u, idx_b, rho_cache, mw_avg_cache)

    #mutating-ish, it mutates du.pressure for a
    face_m_dot = continuity_and_momentum_darcy(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        rho_a, rho_b,
        phys_a, phys_b
    )

    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b
    )

    all_species_advection!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b,
        face_m_dot
    )

    enthalpy_advection!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b,
        face_m_dot
    )

    diffusion_mass_fraction_exchange!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        rho_a, rho_b,
        phys_a, phys_b
    )
end

function solid_solid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    phys_a, phys_b,
    rho_cache, mw_avg_cache,
    change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
)
    area = cell_neighbor_areas[idx_a][face_idx]
    dist = cell_neighbor_distances[idx_a][face_idx]
    norm = cell_neighbor_normals[idx_a][face_idx]

    rho_a = phys_a.rho
    rho_b = phys_b.rho

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b
    )
end

function fluid_solid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    phys_a, phys_b,
    rho_cache, mw_avg_cache,
    change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
)
    area = cell_neighbor_areas[idx_a][face_idx]
    dist = cell_neighbor_distances[idx_a][face_idx]
    norm = cell_neighbor_normals[idx_a][face_idx]

    rho_a = rho_ideal!(u, idx_a, rho_cache, mw_avg_cache)
    rho_b = phys_b.rho

    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b
    )
end

function solid_fluid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    phys_a, phys_b,
    rho_cache, mw_avg_cache,
    change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
)
    area = cell_neighbor_areas[idx_a][face_idx]
    dist = cell_neighbor_distances[idx_a][face_idx]
    norm = cell_neighbor_normals[idx_a][face_idx]

    rho_a = phys_a.rho
    rho_b = rho_ideal!(u, idx_b, rho_cache, mw_avg_cache)

    diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b
    )
end

function methanol_reformer_f_test!(
        du, u, p, t,
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        unconnected_cell_face_map, cell_face_areas, cell_face_normals,
        connection_groups::MethanolReformerConnectionGroups,
        phys::Vector{AbstractPhysics}, cell_phys_id_map::Vector{Int},
        regions_phys_func_cells::Vector{Tuple{AbstractPhysics, Function, Vector{Int}}},
        ax,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    du .= 0.0

    rho_cache = get_tmp(rho_cache, u)
    rho_cache .= 0.0

    mw_avg_cache = get_tmp(mw_avg_cache, u)
    mw_avg_cache .= 0.0

    change_in_molar_concentrations_cache = get_tmp(change_in_molar_concentrations_cache, u)
    change_in_molar_concentrations_cache .= 0.0

    molar_concentrations_cache = get_tmp(molar_concentrations_cache, u) #just using mass fractions for cell 1, this may cause some issues later!
    molar_concentrations_cache .= 0.0

    net_rates_cache = get_tmp(net_rates_cache, u)
    net_rates_cache .= 0.0
    #even though react_cell!() already sets change_in_molar_concentrations_cache to .= 0.0, we're just doing .*= 0.0 to all to make sure



    #= 
    #idea for how we could easily handle fluxes  
    #this would allow us to easily create "hidden" connections for specific types like MethanolReformerPhysics instead of AbstractFluidPhysics
    for (two_phys, connection_function!, group_cells) in connection_groups_two_phys_func_cells
        for (idx_a, idx_a_neighbors) in enumerate(group_cells)
            for (face_idx, idx_b) in enumerate(idx_a_neighbors)
                #the only thing that I don't like about doing this is that whenever an error occurs, the error message just tells you
                #connection_function! errored rather than whatever specific connection_function errored
                #further, profview can't tell you what specific conneciton_function is causing issues
                connecton_function!(
                    du, u, idx_a, idx_b,
                    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                    phys[1], phys[2],
                    rho_cache, mw_avg_cache,
                    change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
                )
            end
        end
    end
    =#

    #instead of updating rho in each flux and internal physics, we could do an initial property retrieval loop
    #the only disadvantage is that once again, solids can't do rho_ideal!() and just have phys.rho
    #=
    @batch for cell_id in eachindex(cell_volumes)
        mw_avg!(u, cell_id, molecular_weights, mw_avg_cache)
        rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)
    end
    =#

    #connections loops
    #look into storing this data into a CSR in the future
    @batch for (idx_a, neighbors) in connection_groups.fluid_fluid
        phys_a = phys[cell_phys_id_map[idx_a]]
        for (idx_b, face_idx) in neighbors
            phys_b = phys[cell_phys_id_map[idx_b]]

            fluid_fluid_flux!(
                du, u, idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                phys_a, phys_b,
                rho_cache, mw_avg_cache,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
            )
        end
    end

    @batch for (idx_a, neighbors) in connection_groups.solid_solid
        phys_a = phys[cell_phys_id_map[idx_a]]
        for (idx_b, face_idx) in neighbors
            phys_b = phys[cell_phys_id_map[idx_b]]

            solid_solid_flux!(
                du, u, idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                phys_a, phys_b,
                rho_cache, mw_avg_cache,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
            )
        end
    end

    @batch for (idx_a, neighbors) in connection_groups.fluid_solid
        phys_a = phys[cell_phys_id_map[idx_a]]
        for (idx_b, face_idx) in neighbors
            phys_b = phys[cell_phys_id_map[idx_b]]

            fluid_solid_flux!(
                du, u, idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                phys_a, phys_b,
                rho_cache, mw_avg_cache,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
            )
        end
    end

    @batch for (idx_a, neighbors) in connection_groups.solid_fluid
        phys_a = phys[cell_phys_id_map[idx_a]]
        for (idx_b, face_idx) in neighbors
            phys_b = phys[cell_phys_id_map[idx_b]]

            solid_fluid_flux!(
                du, u, idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                phys_a, phys_b,
                rho_cache, mw_avg_cache,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
            )
        end
    end

    # ---- Internal Physics, Sources, Boundary Conditions, and Capacities ----
    @batch for (region_phys, region_function!, region_cells) in regions_phys_func_cells
        for cell_id in region_cells
            region_function!(
                du, u, cell_id,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
                rho_cache, mw_avg_cache,
                cell_volumes[cell_id], region_phys
            )
        end
    end
end

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

#= saving this for later if we ever stop using functional BCs
    for field in propertynames(fixed_idxs_and_vals_map)
        #these views are the only way I can seem to get 
        u_field = view(u, field)
        du_field = view(du, field)

        for (cell_id, value) in fixed_idxs_and_vals_map[field]
            if ndims(u_field) > 1 # for mass fractions
                u_field[:, cell_id] .= value
                du_field[:, cell_id] .= 0.0
            else
                u_field[cell_id] = value
                du_field[cell_id] = 0.0
            end
        end
    end
=#

#= also saving previous implementation of Sources
    for cell_id in cell_groups.mass_fraction_sources
        cell_phys_id = cell_phys_id_map[cell_id]

        vol = cell_volumes[cell_id]

        S = phys[cell_phys_id].mass_fraction_vol_source_term * vol

        du.mass_fractions[:, cell_id] += S
    end

    for cell_id in cell_groups.pressure_sources
        cell_phys_id = cell_phys_id_map[cell_id]

        vol = cell_volumes[cell_id]

        S = phys[cell_phys_id].pressure_vol_source_term * vol

        du.pressure[cell_id] += S
    end

    for cell_id in cell_groups.heat_sources
        cell_phys_id = cell_phys_id_map[cell_id]

        vol = cell_volumes[cell_id]

        S = phys[cell_phys_id].heat_vol_source_term * vol

        du.temp[cell_id] += S
    end
=#

#= legacy scaling loops 
    for cell_id in cell_groups.pressure_scaling
        du.pressure[cell_id] /= (1 / (R_gas * u.temp[cell_id]))
    end

    for cell_id in cell_groups.heat_scaling
        cell_phys_id = cell_phys_id_map[cell_id]

        vol = cell_volumes[cell_id]

        rho = get_cell_rho(u, phys[cell_phys_id], cell_id)
        cp = phys[cell_phys_id].cp

        cell_mass = rho * vol
        cap = cp * cell_mass

        du.temp[cell_id] /= cap
    end
=#
