function cap_evaporation_rate_to_phase_holdup!(du, u, cell_id, vol)
    du.liquid_holdup[cell_id] -= (du.mass_evaporated[cell_id] / u.liquid_rho[cell_id]) / vol
    du.gas_holdup[cell_id] += (du.mass_evaporated[cell_id] / u.gas_rho[cell_id]) / vol
end

function vaporization_model!(du, u, cell_id, vol)
    if u.liquid_holdup[cell_id] > 0.0 && u.temp[cell_id] > u.saturation_temp[cell_id]
        du.mass_evaporated[cell_id] += u.mass_transfer_coeff_vap[cell_id] * u.packing_surface_area[cell_id] * (u.temp[cell_id] - u.saturation_temp[cell_id])
        #du.mass[cell_id] -= du.mass_evaporated[cell_id]

        du.heat[cell_id] -= du.mass_evaporated[cell_id] * u.heat_of_vaporization[cell_id]
        
        for_fields!(du.mass_fractions) do species, du_mass_fractions
            # Here we assume the evaporated mass has the composition of the liquid feed
            # In a more complex model, you'd use the liquid-phase mass fractions directly
            du_mass_fractions[species[cell_id]] += (du.mass_evaporated[cell_id] * u.liquid_feed_mass_fractions[species[cell_id]]) / (u.rho[cell_id] * vol)
        end
    end
end






