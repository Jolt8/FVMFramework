const R_gas = 8.314

function upwind(var_a, var_b, mass_flow_rate)
    if mass_flow_rate > 0.0
        return var_a
    else
        return var_b
    end
end

function harmonic_mean(val_a, val_b)
    2 * val_a * val_b / (val_a + val_b)
end

function mw_avg!(u, cell_id, molecular_weights, mw_avg_cache)
    for i in eachindex(molecular_weights)
        mw_avg_cache[cell_id] += u.mass_fractions[i, cell_id] / molecular_weights[i]
    end

    mw_avg_cache[cell_id] = mw_avg_cache[cell_id]^-1.0
end

function rho_ideal!(
        u, #u values
        cell_id,
        rho_cache,
        mw_avg_cache
    )
    rho_cache[cell_id] = (u.pressure[cell_id] * mw_avg_cache[cell_id]) / (R_gas * u.temp[cell_id])
end

function cell_rho!(u, phys::AbstractSolidPhysics, cell_id, rho_cache, mw_avg_cache)
    return rho_cache[cell_id] = phys.rho
end

function cell_rho!(u, phys::AbstractFluidPhysics, cell_id, rho_cache, mw_avg_cache)
    mw_avg!(u, cell_id, phys.species_molecular_weights, mw_avg_cache)
    rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)
end

function get_cell_cp(
    mass_fractions, #u values
    species_cps #other props
)
    cp_avg_cache_for_cell = 0.0

    for i in eachindex(mass_fractions)
        cp_avg_cache_for_cell += mass_fractions[i] * species_cps[i]
    end

    return cp_avg_cache_for_cell
end

function get_partial_pressures(
    mass_fractions, P_total_bar, #u values
    molecular_weights #other props
)
    moles = mass_fractions ./ molecular_weights
    total_moles = sum(moles)
    mole_fractions = moles ./ total_moles

    return mole_fractions .* P_total_bar
end

function van_t_hoff(A, dH, T)
    #K = A * exp(dH/RT)
    return A * exp(-dH / (R_gas * T))
end

function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    K_ref = exp(-ΔG_rxn_ref / (R_gas * T_ref))

    ln_K_ratio = (-ΔH_rxn_ref / R_gas) * (1 / T_actual - 1 / T_ref)

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

#probably unecessary, but whatever
function arrenhius_k(A, Ea, T)
    return A * exp(-Ea / (R_gas * T))
end
