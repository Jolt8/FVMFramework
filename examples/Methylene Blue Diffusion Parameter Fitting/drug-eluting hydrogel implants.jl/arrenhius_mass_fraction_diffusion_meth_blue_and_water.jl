function diffusion_arrenhius_equation(du, u, cell_id, temp_avg)
    return u.diffusion_pre_exponential_factor[1] * exp(-u.diffusion_activation_energy[1] / (R_gas * temp_avg))
end

#the reason these look so weird is that optimizing a the diffusion_pre_exponential_factor and diffusion_activation_energy
#for both water and methylene blue would give 4 parameters to optimize, which is too many for the data we have.
#so instead, we are only optimizing the diffusion_pre_exponential_factor and diffusion_activation_energy
#for methylene blue, and assuming that each m3 of methylene_blue diffused causes an equal volume of water to diffuse
#in the opposite direction.

function arrenhius_mass_fraction_diffusion_meth_blue_and_water!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    temp_avg = 0.5 * (u.temp[idx_a] + u.temp[idx_b])

    meth_blue_species_diffusion_coefficient = diffusion_arrenhius_equation(du, u, idx_a, temp_avg)

    meth_blue_concentration_gradient = (u.mass_fractions.methylene_blue[idx_b] - u.mass_fractions.methylene_blue[idx_a]) / dist
    meth_blue_diffusion = -rho_avg * meth_blue_species_diffusion_coefficient * meth_blue_concentration_gradient * area

    #I should probably check the units here

    #we should probably divide by rho * vol in a capacity function instead later
    #if du.mass_fractions.methylene_blue[idx_a] > 1.0
        #println(du.mass_fractions.methylene_blue[idx_a])
    #end
    #if du.mass_fractions.wateridx_a] > 1.0
        #println(du.mass_fractions.water[idx_a])
    #end
    du.mass_fractions.methylene_blue[idx_a] -= meth_blue_diffusion / (u.rho[idx_a] * vol_a)
    du.mass_fractions.water[idx_a] += meth_blue_diffusion / (u.rho[idx_a] * vol_a)
    
    #=
    if idx_a == 64 
        println("64 as idx_a was applied to")
    elseif idx_a == 1008
        println("1008 as idx_a was applied to")
    end

    if idx_b == 64 
        println("64 as idx_b was applied to")
    elseif idx_b == 1008
        println("1008 as idx_b was applied to")
    end
    =#
end