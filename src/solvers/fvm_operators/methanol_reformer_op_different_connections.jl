struct MethanolReformerPhysics <: AbstractFluidPhysics
    k::Float64
    cp::Float64
    mu::Float64
    permeability::Float64
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
    fluid_fluid::Vector{Tuple{Int,Tuple{Int,Int}}} #connection idx, (cell idx a, cell idx b)
    solid_solid::Vector{Tuple{Int,Tuple{Int,Int}}}
    fluid_solid::Vector{Tuple{Int,Tuple{Int,Int}}}
end

function methanol_reformer_init_conn_groups()
    return MethanolReformerConnectionGroups(
        Tuple{Int,Tuple{Int,Int}}[],
        Tuple{Int,Tuple{Int,Int}}[],
        Tuple{Int,Tuple{Int,Int}}[]
    )
end

#this fucking sucks, but I can't think of anything better, there's got to be a way to leverage 
#dynamic dispatch, but I can't think of anything
#perhaps a function like apply_flux!(args) that takes in all the things needed to do a 
#flux calculation, but dynamically dispatches depending on the type of a_phys and b_phys
#the only issue with this is that instead of allowing the CPU to perform the same operations for a set like 
#fluid_fluid, the methods would randomly change each time (bad for SIMD I think)
#another method would be to create a function that acts like connection catagorizer below, but pushes depending on the types
#ex. :
#=
function connection_catagorizer(connection_groups::MethanolReformerConnectionGroups, conn_idx, idx_a, idx_b, a_phys_type::Type{<:AbstractFluidPhysics}, b_phys_type::Type{<:AbstractFluidPhysics})
    push!(connection_groups.fluid_fluid, (conn_idx, (idx_a, idx_b)))
end

function connection_catagorizer(connection_groups::MethanolReformerConnectionGroups, conn_idx, idx_a, idx_b, a_phys_type::Type{<:AbstractSolidPhysics}, b_phys_type::Type{<:AbstractSolidPhysics})
    push!(connection_groups.solid_solid, (conn_idx, (idx_a, idx_b)))
end

function connection_catagorizer(connection_groups::MethanolReformerConnectionGroups, conn_idx, idx_a, idx_b, a_phys_type::Type{<:AbstractFluidPhysics}, b_phys_type::Type{<:AbstractSolidPhysics})
    push!(connection_groups.fluid_solid, (conn_idx, (idx_a, idx_b)))
end

function connection_catagorizer(connection_groups::MethanolReformerConnectionGroups, conn_idx, idx_a, idx_b, a_phys_type::Type{<:AbstractSolidPhysics}, b_phys_type::Type{<:AbstractFluidPhysics})
    push!(connection_groups.solid_fluid, (conn_idx, (idx_a, idx_b)))
end
=#
# I honestly think the above is worse

function connection_catagorizer!(connection_groups::MethanolReformerConnectionGroups, conn_idx, idx_a, idx_b, type_a, type_b)
    if type_a <: AbstractFluidPhysics && type_b <: AbstractFluidPhysics
        push!(connection_groups.fluid_fluid, (conn_idx, (idx_a, idx_b)))
    elseif type_a <: AbstractSolidPhysics && type_b <: AbstractSolidPhysics
        push!(connection_groups.solid_solid, (conn_idx, (idx_a, idx_b)))
    elseif (type_a <: AbstractSolidPhysics && type_b <: AbstractFluidPhysics) ||
           (type_a <: AbstractFluidPhysics && type_b <: AbstractSolidPhysics)
        push!(connection_groups.fluid_solid, (conn_idx, (idx_a, idx_b)))
    end
end

function methanol_reformer_f_test!(
        du, u, p, t,
        cell_neighbor_map,
        cell_volumes, cell_centroids,
        connection_areas, connection_normals, connection_distances,
        #cell volumes and cell centroids are accessed at the id of the cell
        unconnected_areas,
        #connection areas, normals, and distances are simply accessed by their location in the 
        #list which corresponds to the respective connection in cell_neighbor_map

        connection_groups::MethanolReformerConnectionGroups, phys::Vector{AbstractPhysics}, cell_phys_id_map::Vector{Int},
        regions_phys_func_cells::Vector{Tuple{AbstractPhysics, Function, Vector{Int}}},
        ax,
    )

    #A_Ea_pairs = eachcol(reshape(p, :, n_reactions))
    #unflattened_p would be [[reaction_1_kf_A, reaction_1_kf_Ea], [reaction_2_kf_A, reaction_2_kf_Ea], etc..] 

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)

    du .= 0.0

    n_cells = length(cell_volumes)

    #connections loop
    for (conn_idx, (idx_a, idx_b)) in connection_groups.fluid_fluid
        phys_a_id = cell_phys_id_map[idx_a]
        phys_b_id = cell_phys_id_map[idx_b]
        #in the future to maybe reduce cache misses, we could access phys by phys[cell_id] instead

        area = connection_areas[conn_idx]
        dist = connection_distances[conn_idx]
        norm = connection_normals[conn_idx]

        rho_a = get_cell_rho(u, phys[phys_a_id], idx_a)
        rho_b = get_cell_rho(u, phys[phys_b_id], idx_b)

        #mutating-ish, it mutates du.pressure for a and b
        face_m_dot = continuity_and_momentum_darcy(
            du, u, #if we need views later view(du, :)
            idx_a, idx_b,
            area, norm, dist,
            rho_a, rho_b,
            phys[phys_a_id], phys[phys_b_id],
        )

        diffusion_temp_exchange!(
            du, u,
            idx_a, idx_b,
            area, dist,
            phys[phys_a_id], phys[phys_b_id],
        )

        all_species_advection!(
            du, u,
            idx_a, idx_b,
            face_m_dot,
            area, norm, dist
        )

        enthalpy_advection!(
            du, u,
            idx_a, idx_b,
            face_m_dot,
            area, norm, dist,
            phys[phys_a_id], phys[phys_b_id],
        )
    end

    for (conn_idx, (idx_a, idx_b)) in connection_groups.solid_solid
        phys_a_id = cell_phys_id_map[idx_a]
        phys_b_id = cell_phys_id_map[idx_b]

        area = connection_areas[conn_idx]
        dist = connection_distances[conn_idx]
        norm = connection_normals[conn_idx]

        diffusion_temp_exchange!(
            du, u,
            idx_a, idx_b,
            area, dist,
            phys[phys_a_id], phys[phys_b_id],
        )
    end

    for (conn_idx, (idx_a, idx_b)) in connection_groups.fluid_solid
        phys_a_id = cell_phys_id_map[idx_a]
        phys_b_id = cell_phys_id_map[idx_b]

        area = connection_areas[conn_idx]
        dist = connection_distances[conn_idx]
        norm = connection_normals[conn_idx]

        diffusion_temp_exchange!(
            du, u,
            idx_a, idx_b,
            area, dist,
            phys[phys_a_id], phys[phys_b_id],
        )
    end

    #these are necessary because we can't pass in caches very easily so this will have to do for now. 
    #this is not ideal and should be changed in the future
    change_in_molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1]))
    molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1])) #just using mass fractions for cell 1, this may cause some issues later!
    net_rates_cache = zeros(eltype(u.mass_fractions), length(phys[1].chemical_reactions))
    
    # ---- Internal Physics, Sources, Boundary Conditions, and Capacities ----
    for (region_phys, region_function!, region_cells) in regions_phys_func_cells
        for cell_id in region_cells
            vol = cell_volumes[cell_id]
            rho = get_cell_rho(u, region_phys, cell_id)

            region_function!(
                du, u, cell_id,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
                vol, rho, region_phys
            )
        end
    end
end

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
