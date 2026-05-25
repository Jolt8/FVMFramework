
function drift_flux_mass_conservation!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    for_fields!(u.species_gas_mass_per_volume, du.volumetric_gas_mass_generation, u.gas_mass_fractions) do u_species_gas_mass_accumulation_per_volume, du_volumetric_gas_mass_generation, u_gas_mass_fractions

        change_in_gas_mass_per_volume = du_volumetric_gas_mass_generation[species[idx_a]] - 
        ((u.gas_holdup[idx_a] * u.gas_density[idx_a] * u.gas_velocity[idx_a] * u_gas_mass_fractions[species[idx_a]]) / dist)

        u_species_gas_mass_accumulation_per_volume[species[idx_a]] += change_in_gas_mass_per_volume
    end

    for_fields!(du.species_liquid_mass_per_volume, du.volumetric_liquid_mass_generation, u.liquid_mass_fractions) do du_species_liquid_mass_per_volume, du_volumetric_liquid_mass_generation, u_liquid_mass_fractions

        change_in_liquid_mass_per_volume = du_volumetric_liquid_mass_generation[species[idx_a]] - 
        ((u.liquid_holdup[idx_a] * u.liquid_density[idx_a] * u.liquid_velocity[idx_a] * u_liquid_mass_fractions[species[idx_a]]) / dist)

        du_species_liquid_mass_per_volume[species[idx_a]] += change_in_liquid_mass_per_volume

        du.overall_liquid_mass_per_volume[idx_a] += change_in_liquid_mass_per_volume
    end


end

function drift_flux_mixture_momentum!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
end


function drift_flux_energy_balance!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    for_fields!(du.heat, du.volumetric_enthalpy_generation, u.mass_fractions) do du_enthalpy, du_volumetric_enthalpy_generation, u_enthalpy
        du.heat[idx_a] += (du.gas_velocity[idx_a] * u.gas_mass_per_volume[idx_a] * u.gas_enthalpy[idx_a])
        
    end
end

#non flux functions tracking accumulation
du.species_gas_mass_accumulation_per_volume[species[cell_id]] += (gas_holdup[cell_id] * u.gas_density[cell_id] * du.mass[cell_id] / dist) + u.boiling_rate[cell_id]

du.fluid_volume = darcy_law()#we use darcy's law

gas_velocity = u.distribution_parameter[cell_id] * du.fluid_volume[cell_id] + u.local_gas_drift_velocity[cell_id]

liquid_velocity = (j - u.gas_holdup[cell_id] * u.gas_velocity[cell_id]) / (1 - u.gas_holdup[cell_id])
