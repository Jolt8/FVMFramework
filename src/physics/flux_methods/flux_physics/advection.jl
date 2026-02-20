#use get_cell_rho from helper functions to get rho for cell a and b

#=
function species_advection!(
        du, u, idx_a, idx_b, face_idx, species_name,
        area, norm, dist
    )
    upwinded_mass_fractions = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[species_name][idx_a], u.mass_fractions[species_name][idx_b])

    #it seems like the thing that is really slowing this down is accessing or writing to using for property_name in propertynames(something), this needs to be addressed
    du.mass_fractions[species_name][idx_a] -= -du.mass_face[idx_a][face_idx] * upwinded_mass_fractions
    #-du.mass_face[idx_a][face_idx] because species are being taken out of [idx_a] 
    #I'm still keeping it as -= even though both negative signs could be taken out, but I think always thinking that 
    #stuff is being taken out of [idx_a] is more intuitive 
end

function all_species_advection!(
        du, u, idx_a, idx_b, face_idx,
        area, norm, dist
    )
    for species_name in propertynames(u.mass_fractions)
        species_advection!(du, u, idx_a, idx_b, face_idx, species_name, area, norm, dist)
    end
end
=#


function all_species_advection!(
    du, u, idx_a, idx_b, face_idx,
    area, norm, dist
)
    for species_name in propertynames(u.mass_fractions)
        upwinded_mass_fractions = upwind(du, u, idx_a, idx_b, face_idx, u.mass_fractions[species_name][idx_a], u.mass_fractions[species_name][idx_b])

        du.mass_fractions[species_name][idx_a] -= -du.mass_face[idx_a][face_idx] * upwinded_mass_fractions
    end
end



#use get_cell_cp from helper functions to get cp for cell a and b

function enthalpy_advection!(
        du, u, idx_a, idx_b, face_idx,
        area, norm, dist
    )
    cp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.cp[idx_a], u.cp[idx_b])

    temp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.temp[idx_a], u.temp[idx_b])

    energy_flux = -du.mass_face[idx_a][face_idx] * cp_upwinded * temp_upwinded
    #-du.mass_face[idx_a][face_idx] because mass is being taken out of [idx_a] 

    du.heat[idx_a] -= energy_flux
end


