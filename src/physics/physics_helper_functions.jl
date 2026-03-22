function upwind(
    du, u,
    idx_a, idx_b, face_idx,
    var_a, var_b
)
    if ustrip(du.mass_face[idx_a, face_idx]) < 0.0
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
    u.mw_avg[cell_id] *= 0.0

    foreach_field_at!(u.mass_fractions, u.molecular_weights) do species, mass_fractions, molecular_weights
        u.mw_avg[cell_id] += mass_fractions[species[cell_id]] * molecular_weights[species[cell_id]]
    end
end

function rho_ideal!(u, cell_id)
    u.rho[cell_id] = (u.pressure[cell_id] * u.mw_avg[cell_id]) / (u.R_gas[cell_id] * u.temp[cell_id])
end

function molar_concentrations!(u, cell_id)
    foreach_field_at!(u.mass_fractions, u.molecular_weights, u.molar_concentrations) do species, mass_fractions, molecular_weights, molar_concentrations
        molar_concentrations[species[cell_id]] = (u.rho[cell_id] * mass_fractions[species[cell_id]]) / molecular_weights[species[cell_id]]
    end
end

function molar_fractions!(u, cell_id)
    foreach_field_at!(u.mass_fractions, u.molecular_weights, u.molar_fractions) do species, mass_fractions, molecular_weights, molar_fractions
        molar_fractions[species[cell_id]] = mass_fractions[species[cell_id]] / molecular_weights[species[cell_id]]
    end
end

#honestly, this could probably be removed and derrived only in functions that actually need it
function cp_avg!(u, cell_id)
    u.cp_avg[cell_id] *= 0.0
    foreach_field_at!(u.mass_fractions, u.species_cps, u.cp_avg) do species, mass_fractions, species_cps, cp_avg
        cp_avg[species[cell_id]] += mass_fractions[species[cell_id]] * species_cps[species[cell_id]]
    end
end
#=
function partial_pressure(u, cell)
    return mole_fractions .* P_total_bar
end
=#

function van_t_hoff(A, dH, T, R)
    #K = A * exp(dH/RT)
    return A * exp(-dH / (R * T))
end

function K_gibbs_free(u, cell_id, reaction)
    K_ref = exp(-reaction.ref_delta_G[cell_id] / (u.R_gas[cell_id] * reaction.ref_temp[cell_id]))

    ln_K_ratio = (-reaction.heat_of_reaction[cell_id] / u.R_gas[cell_id]) * (1 / u.temp[cell_id] - 1 / reaction.ref_temp[cell_id])

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

#probably unecessary, but whatever
function arrenhius_k(A, Ea, T, R)
    return A * exp(-Ea / (R * T))
end
