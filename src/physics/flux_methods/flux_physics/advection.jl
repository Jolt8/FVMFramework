function all_species_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    for_fields!(u.mass_fractions, du.species_mass_flows) do species, mass_fractions, species_mass_flows
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, mass_fractions[species[idx_a]], mass_fractions[species[idx_b]])
        species_mass_flows[species[idx_a]] += (du.mass_face[idx_a, face_idx] * upwinded_mass_fraction)
    end
end

function enthalpy_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    cp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.cp[idx_a], u.cp[idx_b])

    temp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.temp[idx_a], u.temp[idx_b])

    energy_flux = du.mass_face[idx_a, face_idx] * cp_upwinded * temp_upwinded
    
    du.heat[idx_a] += energy_flux
end


