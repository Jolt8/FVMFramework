function PAM_reforming_react_cell!(du, u, cell_id, vol)
    mw_avg!(u, cell_id)
    rho_ideal!(u, cell_id)
    molar_concentrations!(u, cell_id)

    #we divide by 1e-5 because PAM parameters are typically in bar
    # methanol, water, carbon_monoxide, hydrogen, carbon_dioxide
    conversion_factor = ((R_gas * u.temp[cell_id]) * 1e-5)
    P_CH3OH = u.molar_concentrations.methanol[cell_id] * conversion_factor
    P_H2O = u.molar_concentrations.water[cell_id] * conversion_factor
    P_CO = u.molar_concentrations.carbon_monoxide[cell_id] * conversion_factor
    P_H2 = u.molar_concentrations.hydrogen[cell_id] * conversion_factor
    P_CO2 = u.molar_concentrations.carbon_dioxide[cell_id] * conversion_factor

    K_CH3O = van_t_hoff(
        u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_A.CH3O, 
        u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_dH.CH3O, 
        u.temp[cell_id]
    )
    K_HCOO = van_t_hoff(
        u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_A.HCOO, 
        u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_dH.HCOO, 
        u.temp[cell_id]
    )
    K_OH = van_t_hoff(
        u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_A.OH, 
        u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_dH.OH, 
        u.temp[cell_id]
    )

    term_CH3O = K_CH3O * (P_CH3OH / sqrt(P_H2))
    term_HCOO = K_HCOO * (P_CO2 * sqrt(P_H2))
    term_OH = K_OH * (P_H2O / sqrt(P_H2))

    DEN = 1.0 + term_CH3O + term_HCOO + term_OH

    for (reforming_reaction_name, reforming_reaction) in pairs(u.reactions.reforming_reactions)
        u.net_rates.reforming_reactions[reforming_reaction_name] = 0.0
    end

    #MSR
    k_MSR = u.reactions.reforming_reactions.MSR_rxn.kf_A * exp(-u.reactions.reforming_reactions.MSR_rxn.kf_Ea / (R_gas * u.temp[cell_id]))
    K_eq_MSR = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.MSR_rxn)

    driving_force = 1.0 - ((P_CO2 * P_H2^3) / (K_eq_MSR * P_CH3OH * P_H2O))

    println("MSR rate: ", (k_MSR * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2))

    u.net_rates.reforming_reactions.MSR_rxn = (k_MSR * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2)

    #MD
    k_MD = u.reactions.reforming_reactions.MD_rxn.kf_A * exp(-u.reactions.reforming_reactions.MD_rxn.kf_Ea / (R_gas * u.temp[cell_id]))
    K_eq_MD = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.MD_rxn)
    driving_force = 1.0 - ((P_CO * P_H2^2) / (K_eq_MD * P_CH3OH))

    println("MD rate: ", (k_MD * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2))

    u.net_rates.reforming_reactions.MD_rxn = (k_MD * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2)

    #WGS
    k_WGS = u.reactions.reforming_reactions.WGS_rxn.kf_A * exp(-u.reactions.reforming_reactions.WGS_rxn.kf_Ea / (R_gas * u.temp[cell_id]))
    K_eq_WGS = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.WGS_rxn)
    driving_force = 1.0 - ((P_CO2 * P_H2) / (K_eq_WGS * P_CO * P_H2O))

    println("WGS rate: ", (k_WGS * K_OH * ((P_CO * P_H2O) / sqrt(P_H2)) * driving_force) / (DEN^2))

    u.net_rates.reforming_reactions.WGS_rxn = (k_WGS * K_OH * ((P_CO * P_H2O) / sqrt(P_H2)) * driving_force) / (DEN^2)

    for (species_name, species_id) in pairs(u.species_ids)
        for (reforming_reaction, _) in pairs(u.reactions.reforming_reactions)
            getproperty(du.molar_concentrations, species_name)[cell_id] += u.net_rates.reforming_reactions[reforming_reaction] * u.reactions.reforming_reactions[reforming_reaction].stoich_coeffs[species_name] 
            #hold up, we need to decide if we should do u.net_rates.reactions.reforming_reactions[reaction_name] OR reforming_reaction.net_rate OR u.reactions.net_rates.reforming_reaction
            #the disadvantage of the second one is that we start mixing up cached values (net_rates) and fixed property values (u.reforming_reactions)
        end
        getproperty(du.mass_fractions, species_name)[cell_id] += getproperty(du.molar_concentrations, species_name)[cell_id] * u.molecular_weights[species_name] / u.rho[cell_id]
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
    end

    for (reforming_reaction, _) in pairs(u.reactions.reforming_reactions)
        du.heat[cell_id] += u.net_rates.reforming_reactions[reforming_reaction] * (-u.reactions.reforming_reactions[reforming_reaction].heat_of_reaction) * vol 
        # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
    end
end