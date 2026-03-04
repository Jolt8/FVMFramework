function ergun_pressure_drop!(du, u, cell_id, vol)
    #superficial_mass_velocity = du.mass_face[cell_id][3] / u.pipe_area[cell_id]

    term1 = 150 * (1 - u.bed_void_fraction[cell_id])^2 / u.bed_void_fraction[cell_id]^3 * u.viscosity[cell_id] * u.superficial_mass_velocity[cell_id] / (u.rho[cell_id] * u.catalyst_particle_diameter[cell_id]^2)
    term2 = 1.75 * (1 - u.bed_void_fraction[cell_id]) / u.bed_void_fraction[cell_id]^3 * u.superficial_mass_velocity[cell_id]^2 / (u.catalyst_particle_diameter[cell_id] * u.rho[cell_id])
    pressure_drop = (term1 + term2) * u.pipe_length[cell_id]
    
    du.pressure[cell_id] -= pressure_drop
end