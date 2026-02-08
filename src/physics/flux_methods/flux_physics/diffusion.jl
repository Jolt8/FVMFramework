function species_numerical_flux(rho_avg, diffusion_coeff, mass_fraction_a, mass_fraction_b, area, dist)
    concentration_gradient = (mass_fraction_b - mass_fraction_a) / dist
    diffusion = -rho_avg * diffusion_coeff * concentration_gradient
    return diffusion * area
end

function diffusion_mass_fraction_exchange!(
        du, u, idx_a, idx_b,
        area, norm, dist,
        rho_a, rho_b,
        phys_a, phys_b
    )

    rho_avg = 0.5 * (rho_a + rho_b)

    for i in eachindex(u.mass_fractions[:, idx_a])
        diffusion_coeff_effective = harmonic_mean(phys_a.species_diffusion_coeffs[i], phys_b.species_diffusion_coeffs[i])
        numerical_flux = species_numerical_flux(rho_avg, diffusion_coeff_effective, u.mass_fractions[i, idx_a], u.mass_fractions[i, idx_b], area, dist)
        du.mass_fractions[i, idx_a] -= numerical_flux
    end
end