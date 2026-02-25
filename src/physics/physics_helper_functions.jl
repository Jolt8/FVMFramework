const R_gas = 8.314

function upwind(
    du, u,
    idx_a, idx_b, face_idx,
    var_a, var_b
)
    if du.mass_face[idx_a][face_idx] < 0.0
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
function mw_avg!(u, cell_id)
    u.mw_avg[cell_id] = 0.0

    for (species_name, species_id) in pairs(u.species_ids)
        u.mw_avg[cell_id] += u.mass_fractions[species_name][cell_id] / u.molecular_weights[species_name][cell_id]
    end

    u.mw_avg[cell_id] = u.mw_avg[cell_id]^-1.0

    #println(u.mw_avg[cell_id])
end

function rho_ideal!(u, cell_id)
    #should we even make mw_avg a cache variable?
    u.rho[cell_id] = (u.pressure[cell_id] * u.mw_avg[cell_id]) / (R_gas * u.temp[cell_id])
end

function molar_concentrations!(u, cell_id)
    for (species_name, species_id) in pairs(u.species_ids)
        u.molar_concentrations[species_name][cell_id] = (u.rho[cell_id] * u.mass_fractions[species_name][cell_id]) / u.molecular_weights[species_name][cell_id]
    end
end

function molar_fractions!(u, cell_id)
    for (species_name, species_id) in pairs(u.species_ids)
        view(u.molar_fractions, species_name)[cell_id] = view(u.mass_fractions, species_name)[cell_id] / u.molecular_weights[species_name]
    end
end

#honestly, this could probably be removed and derrived only in functions that actually need it
function cp_avg!(u, cell_id)
    for (species_name, species_id) in pairs(u.species_ids)
        u.cp_avg[cell_id] += view(u.mass_fractions, species_name)[cell_id] * u.species_cps[species_name]
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

function K_gibbs_free(u, cell_id, reaction)
    K_ref = exp(-reaction.ref_delta_G / (R_gas * reaction.ref_temp))

    ln_K_ratio = (-reaction.heat_of_reaction / R_gas) * (1 / u.temp[cell_id] - 1 / reaction.ref_temp)

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

#probably unecessary, but whatever
function arrenhius_k(A, Ea, T)
    return A * exp(-Ea / (R_gas * T))
end