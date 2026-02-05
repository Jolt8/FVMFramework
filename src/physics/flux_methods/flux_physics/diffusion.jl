function species_numerical_flux(rho_avg, diffusion_coeff, mass_fraction_a, mass_fraction_b, area, dist) #not sure what k_average would represent in respect to diffusion
    concentration_gradient = (mass_fraction_b - mass_fraction_a) / dist
    diffusion = -rho_avg * diffusion_coeff * concentration_gradient
    return diffusion * area
end

function diffusion_mass_fraction_exchange!(
        #mutated vars
        du, u, idx_a, idx_b,
        #geometry data
        area, norm, dist,
        #props
        rho_a, rho_b,
        #diffusion_coeff_a, diffusion_coeff_b, #not sure if these are even necessary
        species_diffusion_coeffs_a, species_diffusion_coeffs_b
        #other data
    )

    rho_avg = 0.5 * (rho_a + rho_b)

    for i in eachindex(u.mass_fractions[:, idx_a])
        diffusion_coeff_effective = harmonic_mean(species_diffusion_coeffs_a[i], species_diffusion_coeffs_b[i])
        numerical_flux = species_numerical_flux(rho_avg, diffusion_coeff_effective, u.mass_fractions[i, idx_a], u.mass_fractions[i, idx_b], area, dist)
        du.mass_fractions[i, idx_a] -= numerical_flux
        du.mass_fractions[i, idx_b] += numerical_flux
    end
end