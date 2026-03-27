function wall_heat_flux!(du, u, cell_id, vol)
    surface_area = pi * u.pipe_inside_diameter[cell_id] * u.pipe_length[cell_id]

    du.heat[cell_id] += u.overall_heat_transfer_coefficient[cell_id] * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area
end