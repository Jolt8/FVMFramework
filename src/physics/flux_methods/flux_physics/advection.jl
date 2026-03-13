function all_species_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    map(keys(du.mass_fractions)) do species_name
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[species_name][idx_a], u.mass_fractions[species_name][idx_b])
        du.species_mass_flows[species_name][idx_a] += (du.mass_face[idx_a][face_idx] * upwinded_mass_fraction)
    end
end

function enthalpy_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    cp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.cp[idx_a], u.cp[idx_b])

    temp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.temp[idx_a], u.temp[idx_b])

    energy_flux = du.mass_face[idx_a][face_idx] * cp_upwinded * temp_upwinded
    
    du.heat[idx_a] += energy_flux
end


