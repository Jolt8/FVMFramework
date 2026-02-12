#use get_cell_rho from helper functions to get rho for cell a and b

function species_advection!(
        du, u, idx_a, idx_b, face_idx, species_name,
        face_m_dot,
        area, norm, dist,
    )
    upwinded_mass_fractions = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[species_name][idx_a], u.mass_fractions[species_name][idx_b])

    du.mass_fractions[species_idx, idx_a] -= face_m_dot * upwinded_mass_fractions
end

function all_species_advection!(
        du, u, idx_a, idx_b, face_idx,
        area, norm, dist
    )
    for species_name in propertynames(u.mass_fractions)
        species_advection!(du, u, idx_a, idx_b, face_idx, species_name, face_m_dot, area, norm, dist)
    end
end

#use get_cell_cp from helper functions to get cp for cell a and b

function enthalpy_advection!(
        du, u, idx_a, idx_b, face_idx,
        area, norm, dist
    )
    cp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.cp[idx_a], u.cp[idx_b])

    temp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.temp[idx_a], u.temp[idx_b])

    energy_flux = -du.mass[idx_a] * cp_upwinded * temp_upwinded
    #-du.mass[idx_a] because mass is being taken out of [idx_a] 

    du.temp[idx_a] -= energy_flux
end


