function fluid_to_steel_pipe_convective_heat_transfer!(
    du, u, 
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    du.heat[idx_a] += u.fluid_to_steel_pipe_convective_heat_transfer_coefficient[idx_a] * (u.temp[idx_b] - u.temp[idx_a]) * area
end