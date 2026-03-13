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

    map(keys(u.mass_fractions)) do species_name
        diffusion_coeff_effective = harmonic_mean(u.diffusion_coefficients[species_name][idx_a], u.diffusion_coefficients[species_name][idx_b])
        numerical_flux = species_numerical_flux(rho_avg, diffusion_coeff_effective, u.mass_fractions[species_name][idx_a], u.mass_fractions[species_name][idx_b], area, dist)
        du.mass_fractions[species_name][idx_a] -= numerical_flux
    end
end