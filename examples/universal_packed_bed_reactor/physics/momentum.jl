function ergun_momentum_friction!(du, u, cell_id, vol)
    # The Ergun equation represents the pressure gradient (friction) in the packed bed.
    # In an FVM framework where pressure is resolved, this acts as a resistive body force (f_ergun).
    # dP/dz = - (A * mu * u + B * rho * u^2)
    
    A = 150.0 * (1.0 - u.bed_void_fraction[cell_id])^2 / (u.particle_diameter[cell_id]^2 * u.bed_void_fraction[cell_id]^3)
    B = 1.75 * (1.0 - u.bed_void_fraction[cell_id]) / (u.particle_diameter[cell_id] * u.bed_void_fraction[cell_id]^3)

    # Calculate the pressure drop over the cell's length
    dp_dz = A * u.dynamic_viscosity[cell_id] * u.superficial_velocity[cell_id] + B * u.rho[cell_id] * abs(u.superficial_velocity[cell_id]) * u.superficial_velocity[cell_id]

    du.pressure[cell_id] -= (dp_dz * u.cell_lengths_along_pipe[cell_id]) 
    #is it more logical to use per_cell_pipe_length or cell_lengths_along_pipe? 
    #wait a minute, this doesn't make any sense because the ergun equation is spatial while this is temporal
    #this is just plain wrong
    #I guess we'll have to pre-calculate the pressure gradient 
end

function update_velocity!(du, u, cell_id, vol)
    # Here we would update velocity_superficial from mass_flow or current state
    # u_superficial = mass_flow / (rho * area)
    u.superficial_velocity[cell_id] = u.pipe_mass_flow[cell_id] / (u.rho[cell_id] * u.pipe_area[cell_id])
end

function new_pressure_driven_mass_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    face_area, face_normal, distance,
)
    # Darcy's Law / Linear Ergun Approximation for mass flow:
    # m_dot = (rho * K * Area / mu) * (dP / dL)
    
    # Permeability K from Ergun Bed Properties
    # We take the average of the two cells for property evaluation
    eps = (u.bed_void_fraction[idx_a] + u.bed_void_fraction[idx_b]) / 2.0
    dp = (u.particle_diameter[idx_a] + u.particle_diameter[idx_b]) / 2.0
    dynamic_viscosity = (u.dynamic_viscosity[idx_a] + u.dynamic_viscosity[idx_b]) / 2.0
    rho = (u.rho[idx_a] + u.rho[idx_b]) / 2.0
    
    K = (dp^2 * eps^3) / (150.0 * (1.0 - eps)^2)
    
    du.mass_face[idx_a, face_idx] -= (rho * K * face_area / dynamic_viscosity) * (u.pressure[idx_a] - u.pressure[idx_b]) / distance
end
