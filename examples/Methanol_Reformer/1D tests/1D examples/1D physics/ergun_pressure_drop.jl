function ergun_pressure_drop!(du, u, cell_id, vol)
    term1 = 150 * (1 - u.bed_void_fraction[cell_id])^2 / u.bed_void_fraction[cell_id]^3 * u.viscosity[cell_id] * u.superficial_mass_velocity[cell_id] / (u.rho[cell_id] * u.catalyst_particle_diameter[cell_id]^2)
    term2 = 1.75 * (1 - u.bed_void_fraction[cell_id]) / u.bed_void_fraction[cell_id]^3 * u.superficial_mass_velocity[cell_id]^2 / (u.catalyst_particle_diameter[cell_id] * u.rho[cell_id])
    pressure_drop = (term1 + term2) * u.pipe_length[cell_id]
    
    du.pressure[cell_id] -= pressure_drop 
    #this is incorrect, ergun is spatial while this is temporal, getting it to be spatial is going to be quite trickly
end