function species_numerical_flux(rho_avg, diffusion_coeff, mass_fraction_a, mass_fraction_b, area, dist)
    concentration_gradient = (mass_fraction_b - mass_fraction_a) / dist
    diffusion = -rho_avg * diffusion_coeff * concentration_gradient
    return diffusion * area
end

function mass_fraction_diffusion!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
)
    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    for_fields!(u.mass_fractions, du.mass_fractions, u.diffusion_coefficients) do species, mass_fractions, du_mass_fractions, diffusion_coefficients
        diffusion_coeff_effective = harmonic_mean(diffusion_coefficients[species[idx_a]], diffusion_coefficients[species[idx_b]])
        numerical_flux = species_numerical_flux(rho_avg, diffusion_coeff_effective, mass_fractions[species[idx_a]], mass_fractions[species[idx_b]], area, dist)
        du_mass_fractions[species[idx_a]] -= numerical_flux
    end
end