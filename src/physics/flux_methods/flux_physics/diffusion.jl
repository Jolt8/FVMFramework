function species_numerical_flux(rho_avg, diffusion_coeff, mass_fraction_a, mass_fraction_b, area, dist)
    concentration_gradient = (mass_fraction_b - mass_fraction_a) / dist
    diffusion = -rho_avg * diffusion_coeff * concentration_gradient
    return diffusion * area
end

function diffusion_mass_fraction_exchange!(
        du, u, idx_a, idx_b, face_idx,
        area, norm, dist,
    )

    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    for species_name in propertynames(u.mass_fractions)
        diffusion_coeff_effective = harmonic_mean(u.species_diffusion_coeffs[species_name][idx_a], u.species_diffusion_coeffs[species_name][idx_b])

        concentration_gradient = (u.mass_fractions[species_name][idx_b] - u.mass_fractions[species_name][idx_a]) / dist
        diffusion = -rho_avg * diffusion_coeff_effective * concentration_gradient
        du.mass_fractions[species_name][idx_a] -= diffusion * area
    end
end