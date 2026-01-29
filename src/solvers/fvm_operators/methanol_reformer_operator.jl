struct MethanolReformerPhysics <: AbstractPhysics
    k::Float64
    rho::Float64
    mu::Float64
    cp::Float64
    permeability::Float64
    species_molecular_weights::Vector{Float64}
    chemical_reactions::Vector{AbstractReaction}
    cell_kg_cat_per_m3_for_each_reaction::Vector{Float64}
    chemical_vol_source_term::Vector{Float64} #chemical addition
    heat_vol_source_term::Float64 #volumetric heating
end

function methanol_reformer_f!(
    du, u, p, t,
    cell_neighbor_map,
    cell_volumes, cell_centroids,
    connection_areas, connection_normals, connection_distances,
    #cell volumes and cell centroids are accessed at the id of the cell
    unconnected_areas,
    #connection areas, normals, and distances are simply accessed by their location in the 
    #list which corresponds to the respective connection in cell_neighbor_map
    phys::Vector{AbstractPhysics}, cell_phys_id_map, fixed_idxs_and_vals_map::ComponentArray, 
    ax, 
    )

    #A_Ea_pairs = eachcol(reshape(p, :, n_reactions))
    #unflattened_p would be [[reaction_1_kf_A, reaction_1_kf_Ea], [reaction_2_kf_A, reaction_2_kf_Ea], etc..] 

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)

    du .= 0.0

    n_cells = length(cell_volumes)

    #connections loop
    for (i, (idx_a, idx_b)) in enumerate(cell_neighbor_map)
        phys_a_id = cell_phys_id_map[idx_a]
        phys_b_id = cell_phys_id_map[idx_b]

        k_a = phys[phys_a_id].k #maybe use heat_phys later
        k_b = phys[phys_b_id].k
        #would it be better to index phys[] by [idx_a] instead? I think this would reduce cache misses

        species_molecular_weights_a = phys[phys_a_id].species_molecular_weights
        species_molecular_weights_b = phys[phys_b_id].species_molecular_weights

        area = connection_areas[i]
        dist = connection_distances[i]
        norm = connection_normals[i]

        du_mass_fractions_a = view(du.mass_fractions, :, idx_a)
        du_mass_fractions_b = view(du.mass_fractions, :, idx_b)
        mass_fractions_a = view(u.mass_fractions, :, idx_a)
        mass_fractions_b = view(u.mass_fractions, :, idx_b)

        mw_avg_a = get_mw_avg(mass_fractions_a, species_molecular_weights_a)
        mw_avg_b = get_mw_avg(mass_fractions_b, species_molecular_weights_b)

        rho_a = cell_rho_ideal(u.pressure[idx_a], u.temp[idx_a], mw_avg_a)
        rho_b = cell_rho_ideal(u.pressure[idx_b], u.temp[idx_b], mw_avg_b)

        mu_a = phys[phys_a_id].mu
        mu_b = phys[phys_b_id].mu

        #mutating-ish, it mutates du.pressure for a and b
        face_m_dot = continuity_and_momentum_darcy(
            @view(du.pressure[idx_a]), @view(du.pressure[idx_b]),
            u.pressure[idx_a], u.pressure[idx_b],
            area, norm, dist,
            cell_volumes[idx_a], cell_volumes[idx_b],
            rho_a, rho_b,
            mu_a, mu_b,
            phys[phys_a_id].permeability
        )

        diffusion_temp_exchange!(
            @view(du.temp[idx_a]), @view(du.temp[idx_b]),
            u.temp[idx_a], u.temp[idx_b],
            area, dist,
            k_a, k_b,
        )

        for i in eachindex(mass_fractions_a)
            species_advection!(
                @view(du_mass_fractions_a[i]), @view(du_mass_fractions_b[i]),
                face_m_dot,
                mass_fractions_a[i], mass_fractions_b[i],
                area, norm, dist
            )
        end

        enthalpy_advection!(
            @view(du.temp[idx_a]), @view(du.temp[idx_b]),
            face_m_dot,
            u.temp[idx_a], u.temp[idx_b],
            area, norm, dist,
            phys[phys_a_id].cp, phys[phys_b_id].cp
        )
    end

    #these are necessary because we can't pass in caches very easily so this will have to do for now. 
    #this is not ideal and should be changed in the future
    change_in_molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1]))
    molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1])) #just using mass fractions for cell 1, this may cause some issues later!
    net_rates_cache = zeros(eltype(u.mass_fractions), length(phys[1].chemical_reactions))

    # Source and Capacity Loop
    for cell_id in 1:n_cells
        #basic phys
        vol = cell_volumes[cell_id]
        cell_phys_id = cell_phys_id_map[cell_id]

        species_molecular_weights = phys[cell_phys_id].species_molecular_weights

        mw_avg_cell = get_mw_avg(u.mass_fractions[:, cell_id], species_molecular_weights)

        k = phys[cell_phys_id].k
        rho = cell_rho_ideal(u.pressure[cell_id], u.temp[cell_id], mw_avg_cell)
        cp = phys[cell_phys_id].cp
        cell_kg_cat_per_m3_for_each_reaction = phys[cell_phys_id].cell_kg_cat_per_m3_for_each_reaction
        reactions = phys[cell_phys_id].chemical_reactions

        #chemical reactions loop
        #mass_fractions = u.mass_fractions[:, cell_id]
        cell_temp = u.temp[cell_id]
        du_mass_fractions  = view(du.mass_fractions, :, cell_id)
        mass_fractions = view(u.mass_fractions, :, cell_id) #we should maybe use views here, probably does matter that much

        #react_cell! also adds the heat of reaction to the cell 
        react_cell!(
            du_mass_fractions, @view(du.temp[cell_id:cell_id]),
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
            mass_fractions, cell_temp,
            vol,
            rho,
            species_molecular_weights, reactions, cell_kg_cat_per_m3_for_each_reaction
        )

        # ----- Source Loop ----- 
        S = phys[cell_phys_id].heat_vol_source_term * vol

        du.temp[cell_id] += S

        # ----- Capacity Loop -----
        cell_mass = rho * vol
        cap = cp * cell_mass

        du.temp[cell_id] /= cap
        du.pressure[cell_id] /= (1 / (R_gas * u.temp[cell_id]))
    end

    # ----- Dirichlet Loop ----- 
    for field in propertynames(fixed_idxs_and_vals_map)
        for (cell_id, value) in fixed_idxs_and_vals_map[field]
            du[field][cell_id] = value
        end
    end
end