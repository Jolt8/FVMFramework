
struct PowerLawReaction <: AbstractReaction
    heat_of_reaction::Float64
    delta_gibbs_free_energy::Float64
    K_gibbs_free_ref_temp::Float64
    kf_A::Float64 # Pre-exponential factor for forward reaction
    kf_Ea::Float64 # Activation energy for forward reaction
    reactants::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    reactant_stoich_coeffs::Vector{Int} #(stoich_coeff_for_reactant_1, stoich_coeff_for_reactant_2, etc...)
    products::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    product_stoich_coeffs::Vector{Int} #[stoich_coeff_for_product_1, stoich_coeff_for_product_2, etc...]
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoich_for_chemical_id_2, etc...) #-1 = reactant, 1 = product
end

function net_reaction_rate(chemical_reaction::PowerLawReaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    kf = (kf_A * exp(-kf_Ea / (R_gas * T)))
    kr = (kr_A * exp(-kr_Ea / (R_gas * T)))

    forward_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.reactants)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.reactant_stoich_coeffs[i]
        forward_term *= concentration^stoich_coeff
    end

    reverse_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.products)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.product_stoich_coeffs[i]
        reverse_term *= concentration^stoich_coeff
    end

    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))

    return net_reaction_rate
end

function react_cell!(
        #input idxs 
        cell_id,
        #mutated vars
        du,
        #u data
        u,
        #caches (also mutated)
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache,
        #geometry data
        vol,
        #props
        rho,
        #phys data
        phys
    )
    change_in_molar_concentrations_cache .= 0.0

    for (species_id, species_mass_fraction) in enumerate(u.mass_fractions[:, cell_id])
        molar_concentrations_cache[species_id] = (rho * species_mass_fraction) / phys.species_molecular_weights[species_id]
    end

    for (reaction_id, reaction) in enumerate(phys.chemical_reactions)
        kf_A = reaction.kf_A
        kf_Ea = reaction.kf_Ea

        #find reverse pre exponential_factor
        K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, u.temp[cell_id], reaction.delta_gibbs_free_energy, reaction.heat_of_reaction)

        kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (R_gas * u.temp[cell_id]))

        #find reverse Ea
        kr_Ea = kf_Ea - reaction.heat_of_reaction

        net_rates_cache[reaction_id] = net_reaction_rate(reaction, molar_concentrations_cache, u.temp[cell_id], kf_A, kf_Ea, kr_A, kr_Ea) * phys.cell_kg_cat_per_m3_for_each_reaction[reaction_id]
        #rate returned by net_reactions_rate is in mol / (kg_cat * s), so we times by the cell's kg_cat / m3
        #rate is now in mol / (m3 * s)
    end

    for (species_id, species_molar_concentration) in enumerate(molar_concentrations_cache)
        for (reaction_id, reaction) in enumerate(phys.chemical_reactions)
            stoich = reaction.all_stoich_coeffs[species_id]
            change_in_molar_concentrations_cache[species_id] += net_rates_cache[reaction_id] * stoich
        end

        du.mass_fractions[species_id, cell_id] += (change_in_molar_concentrations_cache[species_id] * phys.species_molecular_weights[species_id]) / rho
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
    end

    for (reaction_id, reaction) in enumerate(phys.chemical_reactions)
        du.temp[cell_id] += net_rates_cache[reaction_id] * (-reaction.heat_of_reaction) * vol
        # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
    end
end