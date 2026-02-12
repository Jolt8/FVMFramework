
function get_darcy_mass_flux(rho_avg, permeability, viscosity, pressure_a, pressure_b, area, dist) #not sure what k_average would represent in respect to diffusion
    p_grad = (pressure_b - pressure_a) / dist
    m_dot = -rho_avg * (permeability / viscosity) * p_grad * area
    return m_dot
end

function continuity_and_momentum_darcy!(
    du, u, idx_a, idx_b, face_idx,
    area, norm, dist,
    )

    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])
    mu_avg = 0.5 * (u.mu[idx_a] + u.mu[idx_b])
    permeability_avg = 0.5 * (u.permeability[idx_a] + u.permeability[idx_b])

    p_grad = (u.pressure[idx_b] - u.pressure[idx_a]) / dist
    face_m_dot = -rho_avg * (permeability_avg / mu_avg) * p_grad * area

    du.pressure[idx_a] -= face_m_dot
    du.mass[idx_a] -= face_m_dot
end

