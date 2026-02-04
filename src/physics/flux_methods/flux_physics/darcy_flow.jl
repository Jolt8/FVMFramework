
function get_darcy_mass_flux(rho_avg, permeability, viscosity, pressure_a, pressure_b, area, dist) #not sure what k_average would represent in respect to diffusion
    p_grad = (pressure_b - pressure_a) / dist
    m_dot = -rho_avg * (permeability / viscosity) * p_grad * area
    return m_dot
end

function continuity_and_momentum_darcy(
    #NOTE!!:
    #this also returns face_m_dot even though it's a f!() function
    du, u,
    idx_a, idx_b,
    #geometry data
    area, norm, dist,
    #other props 
    rho_a, rho_b, #kinda a u value because it changes with time but not explicitly tracked through u values
    phys_a, phys_b
    )

    pressure_a = u.pressure[idx_a]
    pressure_b = u.pressure[idx_b]

    rho_avg = 0.5 * (rho_a + rho_b)
    mu_avg = 0.5 * (phys_a.mu + phys_b.mu)
    permeability_avg = 0.5 * (phys_a.permeability + phys_b.permeability)

    face_m_dot = get_darcy_mass_flux(rho_avg, permeability_avg, mu_avg, pressure_a, pressure_b, area, dist)

    #perhaps we should save the vol_a and vol_b division for the capacity functions later
    du.pressure[idx_a] -= face_m_dot
    du.pressure[idx_b] += face_m_dot

    return face_m_dot
end

