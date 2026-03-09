function Nu_packed_bed_Gnielinski!(du, u, cell_id, Re, Pr)
    #borrowed from ht.py (love that package)
    #Calculates Nusselt number of a fluid passing over a bed of particles
    
    Nu_lam = 0.664 * Re^0.5 * Pr^(1/3.)
    
    Nu_turb = 0.037 * Re^0.8 * Pr/(1 + 2.443 * Re^-0.1 * (Pr^(2/3.) - 1))
    
    Nu_sphere = 2 + (Nu_lam^2 + Nu_turb^2)^0.5
    
    fa = 1.0 + 1.5 * (1.0 - u.bed_void_fraction[cell_id])

    Nu = Nu_sphere * fa
    
    return Nu
end

function UA_prime(du, u, cell_id, vol)
    #println("u.rho[cell_id] = $(u.rho[cell_id])")
    #println("u.velocity[cell_id] = $(u.velocity[cell_id])")
    #println("u.pipe_hydraulic_diameter[cell_id] = $(u.pipe_hydraulic_diameter[cell_id])")
    #println("u.viscosity[cell_id] = $(u.viscosity[cell_id])")
    #println("u.cp[cell_id] = $(u.cp[cell_id])")
    #println("u.k[cell_id] = $(u.k[cell_id])")
    Re = (u.rho[cell_id] * u.velocity[cell_id] * u.pipe_hydraulic_diameter[cell_id]) / (u.viscosity[cell_id])
    #println("Re = $(Re)")
    Pr = (u.viscosity[cell_id] * u.cp[cell_id]) / (u.k[cell_id])
    #println("Pr = $(Pr)")

    pipe_nusselt_number = Nu_packed_bed_Gnielinski!(du, u, cell_id, Re, Pr)

    pipe_heat_transfer_coeff = (pipe_nusselt_number * u.cp[cell_id]) / u.pipe_hydraulic_diameter[cell_id]

    conv_resistance = 1 / (pipe_heat_transfer_coeff * u.pipe_hydraulic_diameter[cell_id])

    total_cond_resistance = u.pipe_thickness[cell_id] / (u.pipe_k[cell_id] * u.pipe_hydraulic_diameter[cell_id])

    UA_prime = 1 / (conv_resistance + total_cond_resistance)

    return UA_prime
end