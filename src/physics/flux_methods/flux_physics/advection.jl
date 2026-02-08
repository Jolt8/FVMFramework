#use get_cell_rho from helper functions to get rho for cell a and b

function species_advection!(
        du, u, idx_a, idx_b, species_idx,
        face_m_dot,
        area, norm, dist,
    )
    upwinded_mass_fractions = upwind(u.mass_fractions[species_idx, idx_a], u.mass_fractions[species_idx, idx_b], face_m_dot)

    du.mass_fractions[species_idx, idx_a] -= face_m_dot * upwinded_mass_fractions
end

function all_species_advection!(
        du, u, idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b,
        face_m_dot
    )
    for species_idx in eachindex(view(u.mass_fractions, :, idx_a))
        species_advection!(du, u, idx_a, idx_b, species_idx, face_m_dot, area, norm, dist)
    end
end

#use get_cell_cp from helper functions to get cp for cell a and b

function enthalpy_advection!(
        du, u, idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b,
        face_m_dot
    )
    cp_upwinded = upwind(phys_a.cp, phys_b.cp, face_m_dot)

    temp_upwinded = upwind(u.temp[idx_a], u.temp[idx_b], face_m_dot)

    energy_flux = face_m_dot * cp_upwinded * temp_upwinded

    du.temp[idx_a] -= energy_flux
end


