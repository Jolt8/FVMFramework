
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

#I think I know what is causing the slowdown here, the reaciton dynamic dispatch seems to be causing some issues
#whenever a reactin property is grabbed using reaciton.something, it causes GC

function react_cell!(
        du, u, cell_id,
        vol, 
        phys
    )
    du.molar_concentrations .= 0.0

    for species_id in eachindex(u.molar_concentrations)
        u.molar_concentrations[species_id] = (u.rho[cell_id] * u.mass_fractions[species_id, cell_id]) / u.species_molecular_weights[species_id]
    end

    for (reaction_id, reaction) in enumerate(phys.chemical_reactions)
        kf_A = reaction.kf_A
        kf_Ea = reaction.kf_Ea

        #find reverse pre exponential_factor
        K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, u.temp[cell_id], reaction.delta_gibbs_free_energy, reaction.heat_of_reaction) #this causes GC

        kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (R_gas * u.temp[cell_id])) #this causes GC

        #find reverse Ea
        kr_Ea = kf_Ea - reaction.heat_of_reaction #this causes GC

        u.net_rates[reaction_id] = net_reaction_rate(reaction, u.molar_concentrations, u.temp[cell_id], kf_A, kf_Ea, kr_A, kr_Ea) * phys.cell_kg_cat_per_m3_for_each_reaction[reaction_id] #this causes GC
        #rate returned by net_reactions_rate is in mol / (kg_cat * s), so we times by the cell's kg_cat / m3
        #rate is now in mol / (m3 * s)
    end

    for species_id in eachindex(u.molar_concentrations)
        for (reaction_id, reaction) in enumerate(phys.chemical_reactions)
            du.molar_concentrations[species_id] += u.net_rates[reaction_id] * reaction.all_stoich_coeffs[species_id] #this causes GC
        end
        du.mass_fractions[species_id, cell_id] += (du.molar_concentrations[species_id] * phys.species_molecular_weights[species_id]) / u.rho[cell_id]
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
    end

    for (reaction_id, reaction) in enumerate(phys.chemical_reactions)
        du.temp[cell_id] += u.net_rates[reaction_id] * (-reaction.heat_of_reaction) * vol #this causes GC
        # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
    end
end