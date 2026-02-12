#we should probably move to using cell instead of cell_id inside these functions and cell_a cell_b 
#however in the main loop inside the op function, we should still use cell because cell can sometimes be used
function net_reaction_rate!(du, u, cell, reaction, kr_A, kr_Ea) #reaction is u.reactions[reaction_id]
    kf = (reaction.kf_A * exp(-reaction.kf_Ea / (R_gas * u.temp[cell])))
    kr = (kr_A * exp(-reaction.kr_Ea / (R_gas * u.temp[cell])))

    forward_term = 1.0
    for species_name in propertynames(reaction.all_stoich_coeffs)
        concentration = u.molar_concentrations[species_name][cell]
        stoich_coeff = reaction.reactant_stoich_coeffs[species_name]
        forward_term *= concentration^stoich_coeff
    end

    reverse_term = 1.0
    for species_name in propertynames(reaction.all_stoich_coeffs)
        concentration = u.molar_concentrations[species_name][cell]
        stoich_coeff = reaction.product_stoich_coeffs[species_name]
        reverse_term *= concentration^stoich_coeff
    end

    reaction.net_rate = ((kf * forward_term) - (kr * reverse_term))
end

#I think I know what is causing the slowdown here, the reaciton dynamic dispatch seems to be causing some issues
#whenever a reactin property is grabbed using reaciton.something, it causes GC


#after rewriting this, I think it would be better to have u.species = ComponentVector(molecular_weight = 10.0, mass_fraction = 0.4, molar_concentration = 0.3, diffusion_coefficient = 1e-5, etc.)
#also, I'm absolutely loving writing physics functions while only referencing du and u
#it's absolutely amazing to write out the equations almost in their base form and then only have to keep track where each variable is stored 
#however in the future, I'm sure this will become more and more concrete and it will only become more and more intuitive
function power_law_react_cell!(du, u, cell, reaction, vol)
    du.molar_concentrations .= 0.0 
    #IMPORTANT:
    #oh wait, if we ever decide to make molar_concentrations a state variable (even though mass_fractions is almost always used),
    #we would have to get rid of this line
    #I guess this is one disadvantage as treating caches as pseudo-state variables

    #maybe we don't have to include species_name, we could just do species and index based on u.molar_concentrations[species][cell]
    #this is another instance of weird indexing, we could make it so that u contains species and we do:
    #for species in u.species
    #   species.molar_concentration[cell] = (u.rho[cell] * species.mass_fraction / species.mw)
    #end
    #this updates the molar_concentrations_cache for the cell
    molar_concentrations!(u, cell)

    #find reverse pre exponential_factor
    K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, u.temp[cell], reaction.delta_gibbs_free_energy, reaction.heat_of_reaction) 

    kr_A = (reaction.kf_A / K_ref) * exp(-reaction.heat_of_reaction / (R_gas * u.temp[cell])) 

    #find reverse Ea
    kr_Ea = reaction.kf_Ea - reaction.heat_of_reaction 

    #mutates reaction.net_rate
    net_reaction_rate!(du, u, cell, reaction, kr_A, kr_Ea) * u.kg_cat_per_m3_for_each_reaction[cell] 
    #rate returned by net_reactions_rate is in mol / (kg_cat * s), so we times by the cell's kg_cat / m3
    #rate is now in mol / (m3 * s)

    for species_name in propertynames(u.molar_concentrations)
        du.molar_concentrations[species_name][cell] += reaction.net_rates * reaction.all_stoich_coeffs[species_name] 
        du.mass_fractions[species_name][cell] += (du.molar_concentrations[species_name][cell] * u.species_molecular_weights[species_name]) / u.rho[cell]
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
    end

    du.temp[cell] += u.net_rates[reaction_id] * (-reaction.heat_of_reaction) * vol 
    # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
end