using Clapeyron

function get_mole_fractions_vec(du, u, cell_id, vol)
    gas_mole_fractions_vec = [
        max(u.gas_mole_fractions.methanol[cell_id], 1e-15),
        max(u.gas_mole_fractions.water[cell_id], 1e-15)
    ]
    liquid_mole_fractions_vec = [
        max(u.liquid_mole_fractions.methanol[cell_id], 1e-15),
        max(u.liquid_mole_fractions.water[cell_id], 1e-15)
    ]
    # Normalize so they sum to 1
    gas_mole_fractions_vec ./= sum(gas_mole_fractions_vec)
    liquid_mole_fractions_vec ./= sum(liquid_mole_fractions_vec)
    return gas_mole_fractions_vec, liquid_mole_fractions_vec
end

function multi_species_update_eos_densities!(du, u, cell_id, vol, clapeyron_model)
    # Compute molar volumes (m^3/mol)
    #this will allocate, we'll try views for now, but it will probably error
    gas_mole_fractions_vec, liquid_mole_fractions_vec = get_mole_fractions_vec(du, u, cell_id, vol)
    
    V_l = Clapeyron.volume(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], liquid_mole_fractions_vec, phase = :liquid)
    V_g = Clapeyron.volume(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], gas_mole_fractions_vec, phase = :gas)

    # Update intrinsic densities (with units kg/m^3)
    u.eos_liquid_density[cell_id] = u.liquid_mw_avg[cell_id] / V_l
    u.eos_gas_density[cell_id] = u.gas_mw_avg[cell_id] / V_g
end

function multi_species_update_K_vle!(du, u, cell_id, vol, clapeyron_model)
    # Compute fugacity coefficients for both phases
    gas_mole_fractions_vec, liquid_mole_fractions_vec = get_mole_fractions_vec(du, u, cell_id, vol)
    
    phi_l = Clapeyron.fugacity_coefficient(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], liquid_mole_fractions_vec, phase = :liquid)
    phi_g = Clapeyron.fugacity_coefficient(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], gas_mole_fractions_vec, phase = :gas)

    # Update K-values (K = phi_l / phi_g)
    i = 1
    for_fields!(u.K_vle) do species, u_K_vle
        u_K_vle[species[cell_id]] = phi_l[i] / phi_g[i]
        i += 1
    end

    #u.K_vle.methanol[cell_id] = 0.5
    #u.K_vle.water[cell_id] = 0.5
end

function single_species_update_eos_densities!(du, u, cell_id, vol, clapeyron_model)
    # Compute molar volumes (m^3/mol)
    V_l = Clapeyron.volume(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], phase = :liquid)
    V_g = Clapeyron.volume(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], phase = :gas)

    # Update intrinsic densities (with units kg/m^3)
    u.eos_liquid_density[cell_id] = u.liquid_mw_avg[cell_id] / V_l
    u.eos_gas_density[cell_id] = u.gas_mw_avg[cell_id] / V_g
end

function single_species_update_K_vle!(du, u, cell_id, vol, clapeyron_model)
    # Compute fugacity coefficients for both phases
    phi_l = Clapeyron.fugacity_coefficient(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], phase = :liquid)
    phi_g = Clapeyron.fugacity_coefficient(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], phase = :gas)

    # Update K-values (K = phi_l / phi_g)
    i = 1
    for_fields!(u.K_vle) do species, u_K_vle
        u_K_vle[species[cell_id]] = phi_l[i] / phi_g[i]
        i += 1
    end

    #u.K_vle.methanol[cell_id] = 0.5
    #u.K_vle.water[cell_id] = 0.5
end