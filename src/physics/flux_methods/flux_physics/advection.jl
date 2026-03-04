#use get_cell_rho from helper functions to get rho for cell a and b

function all_species_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    #map(keys(u.mass_fractions)) do species_name
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[:methanol][idx_a], u.mass_fractions[:methanol][idx_b])
        du.mass_fractions[:methanol][idx_a] += du.mass_face[idx_a][face_idx] * upwinded_mass_fraction
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[:water][idx_a], u.mass_fractions[:water][idx_b])
        du.mass_fractions[:water][idx_a] += du.mass_face[idx_a][face_idx] * upwinded_mass_fraction
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[:carbon_monoxide][idx_a], u.mass_fractions[:carbon_monoxide][idx_b])
        du.mass_fractions[:carbon_monoxide][idx_a] += du.mass_face[idx_a][face_idx] * upwinded_mass_fraction
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[:hydrogen][idx_a], u.mass_fractions[:hydrogen][idx_b])
        du.mass_fractions[:hydrogen][idx_a] += du.mass_face[idx_a][face_idx] * upwinded_mass_fraction
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[:carbon_dioxide][idx_a], u.mass_fractions[:carbon_dioxide][idx_b])
        du.mass_fractions[:carbon_dioxide][idx_a] += du.mass_face[idx_a][face_idx] * upwinded_mass_fraction
    #end
end

#use get_cell_cp from helper functions to get cp for cell a and b

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


