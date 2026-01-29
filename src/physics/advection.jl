#use get_cell_rho from helper functions to get rho for cell a and b

function species_advection!(
    #mutated vars
    du_species_mass_fractions_a, du_species_mass_fractions_b,
    #u values
    face_m_dot,
    species_mass_fractions_a, species_mass_fractions_b,
    #geometry data
    area, norm, dist,
    #other props 
    )

    upwinded_mass_fractions = upwind(species_mass_fractions_a, species_mass_fractions_b, face_m_dot)

    du_species_mass_fractions_a[1] -= face_m_dot * upwinded_mass_fractions
    du_species_mass_fractions_b[1] += face_m_dot * upwinded_mass_fractions
end

#use get_cell_cp from helper functions to get cp for cell a and b

function enthalpy_advection!(
    #mutated vars
    du_temp_a, du_temp_b,
    #u values
    face_m_dot,
    temp_a, temp_b,
    #geometry data
    area, norm, dist,
    #other props 
    cp_a, cp_b
    )
    cp_upwinded = upwind(cp_a, cp_b, face_m_dot)

    temp_upwinded = upwind(temp_a, temp_b, face_m_dot)

    energy_flux = face_m_dot * cp_upwinded * temp_upwinded 

    du_temp_a[1] -= energy_flux
    du_temp_b[1] += energy_flux
end


