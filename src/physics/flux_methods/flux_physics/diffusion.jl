function species_numerical_flux(rho_avg, diffusion_coeff, mass_fraction_a, mass_fraction_b, area, dist)
    concentration_gradient = (mass_fraction_b - mass_fraction_a) / dist
    diffusion = -rho_avg * diffusion_coeff * concentration_gradient
    return diffusion * area
end

function diffusion_mass_fraction_exchange!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
)

    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    for (species_name, species_id) in pairs(u.species_ids)
        diffusion_coeff_effective = harmonic_mean(u.diffusion_coefficients[species_name][idx_a], u.diffusion_coefficients[species_name][idx_b])
        
        numerical_flux = species_numerical_flux(rho_avg, diffusion_coeff_effective, u.mass_fractions[species_name][idx_a], u.mass_fractions[species_name][idx_b], area, dist)

        #=
        if u.mass_fractions[species_name][idx_a] !== u.mass_fractions[species_name][idx_b]
            println("rho_avg: ", rho_avg)
            println("diffusion_coeff_effective: ", diffusion_coeff_effective)
            println("mass_fraction_a: ", u.mass_fractions[species_name][idx_a])
            println("mass_fraction_b: ", u.mass_fractions[species_name][idx_b])
            println("area: ", area)
            println("dist: ", dist)
            println("numerical_flux: ", numerical_flux)
            println(typeof(du.mass_fractions))
            println("du_mass_fractions $(du.mass_fractions[species_name][idx_a])")
        end=#

        du.mass_fractions[species_name][idx_a] -= numerical_flux

        #=if u.mass_fractions[species_name][idx_a] !== u.mass_fractions[species_name][idx_b]
            println("du_mass_fractions $(du.mass_fractions[species_name][idx_a])")
        end=#
    end
end