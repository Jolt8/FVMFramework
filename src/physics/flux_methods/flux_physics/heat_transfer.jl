function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

function numerical_flux(k_avg, T_L, T_R, area, dist)
    grad_T = (T_R - T_L) / dist
    q = -k_avg * grad_T
    return q * area
end

function diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        connection_area, connection_distance,
        #props
        phys_a, phys_b
        #other data

    )
    k_effective = get_k_effective(phys_a.k, phys_b.k)

    F = numerical_flux(k_effective, u.temp[idx_a], u.temp[idx_b], connection_area, connection_distance)
    
    du.temp[idx_a] -= F
    du.temp[idx_b] += F
end