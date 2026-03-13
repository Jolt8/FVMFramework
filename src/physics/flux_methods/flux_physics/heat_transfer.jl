function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

function temp_diffusion!(
    du, u, 
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    k_effective = get_k_effective(u.k[idx_a], u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist
    F = -k_effective * grad_T * area

    du.heat[idx_a] -= F
end