function cap_evaporation_rate_to_phase_holdup!(du, u, cell_id, vol)
    du.liquid_holdup[cell_id] -= (du.mass_evaporated[cell_id] / u.liquid_rho[cell_id]) / vol
    du.gas_holdup[cell_id] += (du.mass_evaporated[cell_id] / u.gas_rho[cell_id]) / vol

    #doesn't this cause gas_holdup to exceed 1.0?
    #if so, how do we have the excess gas leave the cell?
end

function vaporization_model!(du, u, cell_id, vol)
    if u.liquid_holdup[cell_id] > 0.0 && u.temp[cell_id] > u.saturation_temp[cell_id] 
        du.mass_evaporated[cell_id] += u.mass_transfer_coeff_vap[cell_id] * u.packing_surface_area[cell_id] * (u.temp[cell_id] - u.saturation_temp[cell_id])
        
        du.mass[cell_id] += du.mass_evaporated[cell_id]
        #oh, I see, mass is conserved through the density decreasing
        #wait, but why is this -=, shouldn't it be += if the density decreases to compensate?
        #ok, yeah, setting this to -= makes the solver unable to solve the system


        du.heat[cell_id] -= du.mass_evaporated[cell_id] * u.heat_of_vaporization[cell_id]
        
        #why is this even necessary?
        
        for_fields!(du.mass_fractions) do species, du_mass_fractions
            # Here we assume the evaporated mass has the composition of the liquid feed
            # In a more complex model, you'd use the liquid-phase mass fractions directly
            du_mass_fractions[species[cell_id]] += (du.mass_evaporated[cell_id] * u.liquid_feed_mass_fractions[species[cell_id]]) / (u.rho[cell_id] * vol)
        end
        
    end
end


function rho_multiphase!(du, u, cell_id, vol)
    total_holdup = u.liquid_holdup[cell_id] + u.gas_holdup[cell_id]
    normalized_liquid_holdup = u.liquid_holdup[cell_id] / total_holdup
    normalized_gas_holdup = u.gas_holdup[cell_id] / total_holdup
    u.rho[cell_id] = normalized_liquid_holdup * u.liquid_rho[cell_id] + normalized_gas_holdup * u.gas_rho[cell_id]
end







