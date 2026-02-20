const R_gas = 8.314

function upwind(du, u, idx_a, idx_b, face_idx, var_a, var_b)
    if du.mass[idx_a] < 0.0
        return var_a
    else
        return var_b
    end
end

function harmonic_mean(val_a, val_b)
    2 * val_a * val_b / (val_a + val_b)
end

#we should probably decide whether or not to pass du into every function just for consistency
#function mw_avg!(du, u, cell)
function mw_avg!(u, cell)
    for species_name in propertynames(u.mass_fractions)
        u.mw_avg[cell] += u.mass_fractions[species_name][cell] / u.molecular_weights[species_name]
    end

    u.mw_avg[cell] = u.mw_avg[cell]^-1.0
end

function rho_ideal!(u, cell)
    u.rho[cell] = (u.pressure[cell] * u.mw_avg[cell]) / (R_gas * u.temp[cell])
end

function molar_concentrations!(u, cell)
    for species_name in propertynames(u.molar_concentrations)
        u.molar_concentrations[species_name][cell] = (u.rho[cell] * u.mass_fractions[species_name][cell]) / u.molecular_weights[species_name]
    end
end

function molar_fractions!(u, cell_id)
    for species_name in propertynames(u.molar_concentrations)
        u.molar_concentrations[species_name][cell_id] = u.mass_fractions[species_name][cell_id] / u.molecular_weights[species_name]
    end
end

#honestly, this could probably be removed and derrived only in functions that actually need it
function cp_avg!(u, cell)
    for species_name in propertynames(u.mass_fractions)
        u.cp_avg[cell] += u.mass_fractions[species_name][cell] * u.species_cps[species_name]
    end
end
#=
function partial_pressure(u, cell)
    return mole_fractions .* P_total_bar
end
=#

function van_t_hoff(A, dH, T)
    #K = A * exp(dH/RT)
    return A * exp(-dH / (R_gas * T))
end

function K_gibbs_free(u, cell, reaction)
    K_ref = exp(-reaction.ΔG_rxn_ref / (8.314e-3 * reaction.T_ref)) #R is in kJ

    ln_K_ratio = (-reaction.ΔH_rxn_ref / 8.314e-3) * (1 / u.temp[cell] - 1 / reaction.T_ref)

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

#probably unecessary, but whatever
function arrenhius_k(A, Ea, T)
    return A * exp(-Ea / (R_gas * T))
end
