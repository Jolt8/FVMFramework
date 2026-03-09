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

    map(keys(u.mass_fractions)) do species_name
        #println("species_name: $(species_name)")
        #println("u.mass_fractions[species_name][cell_id] = $(u.mass_fractions[species_name][cell_id])")
        #println("u.molecular_weights[species_name][cell_id] = $(u.molecular_weights[species_name][cell_id])")
        u.mw_avg[cell_id] += u.mass_fractions[species_name][cell_id] / u.molecular_weights[species_name][cell_id]
    end

    u.mw_avg[cell_id] = u.mw_avg[cell_id]^-1.0
end

function rho_ideal!(u, cell_id)
    #println("u.pressure[cell_id] = $(u.pressure[cell_id])")
    #println("u.mw_avg[cell_id] = $(u.mw_avg[cell_id])")
    #println("u.temp[cell_id] = $(u.temp[cell_id])")
    #println("R_gas = $(R_gas)")
    u.rho[cell_id] = (u.pressure[cell_id] * u.mw_avg[cell_id]) / (R_gas * u.temp[cell_id])
end

function molar_concentrations!(u, cell_id)
    map(keys(u.mass_fractions)) do species_name
        u.molar_concentrations[species_name][cell_id] = (u.rho[cell_id] * u.mass_fractions[species_name][cell_id]) / u.molecular_weights[species_name][cell_id]
    end
end

function molar_fractions!(u, cell_id)
    map(keys(u.mass_fractions)) do species_name
        u.molar_fractions[species_name][cell_id] = u.mass_fractions[species_name][cell_id] / u.molecular_weights[species_name][cell_id]
    end
end

#honestly, this could probably be removed and derrived only in functions that actually need it
function cp_avg!(u, cell_id)
    map(keys(u.mass_fractions)) do species_name
        u.cp_avg[cell_id] += u.mass_fractions[species_name][cell_id] * u.species_cps[species_name][cell_id]
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
    K_ref = exp(-reaction.ref_delta_G[cell_id] / (R_gas * reaction.ref_temp[cell_id]))

    ln_K_ratio = (-reaction.heat_of_reaction[cell_id] / R_gas) * (1 / u.temp[cell_id] - 1 / reaction.ref_temp[cell_id])

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

#probably unecessary, but whatever
function arrenhius_k(A, Ea, T)
    return A * exp(-Ea / (R_gas * T))
end
