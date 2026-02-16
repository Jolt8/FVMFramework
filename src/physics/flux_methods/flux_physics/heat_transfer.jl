function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

#we could either use cell_a, idx_a, or a here. I'm not sure which is best.
function numerical_flux(u, idx_a, idx_b, k_avg, area, dist)
    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist
    q = -k_avg * grad_T
    return q * area
end

function diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        area, norm, dist
    )
    k_effective = get_k_effective(u.k[idx_a], u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist
    q = -k_effective * grad_T
    F = q * area
    
    du.heat[idx_a] -= F
end