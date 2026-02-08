function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

function numerical_flux(k_avg, temp_a, temp_b, area, dist)
    grad_T = (temp_b - temp_a) / dist
    q = -k_avg * grad_T
    return q * area
end

function diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b,
        area, norm, dist,
        phys_a, phys_b
    )
    k_effective = get_k_effective(phys_a.k, phys_b.k)

    F = numerical_flux(k_effective, u.temp[idx_a], u.temp[idx_b], area, dist)
    
    du.temp[idx_a] -= F
end