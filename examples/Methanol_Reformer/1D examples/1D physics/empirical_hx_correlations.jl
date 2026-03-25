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

function overall_heat_transfer_coefficient(du, u, cell_id, vol)
    #@show u.rho[cell_id]
    #@show u.velocity[cell_id]
    #@show u.catalyst_particle_diameter[cell_id]
    #@show u.viscosity[cell_id]
    #@show u.cp[cell_id]
    #@show u.k[cell_id]
    
    Re = (u.rho[cell_id] * u.velocity[cell_id] * u.catalyst_particle_diameter[cell_id]) / (u.viscosity[cell_id])
    Pr = (u.viscosity[cell_id] * u.cp[cell_id]) / (u.k[cell_id])

    #@show Re
    #@show Pr


    if Re <= 0 || Pr <= 0
        @show u.rho[cell_id]
        @show u.velocity[cell_id]
        @show u.catalyst_particle_diameter[cell_id]
        @show u.viscosity[cell_id]
        @show u.cp[cell_id]
        @show u.k[cell_id]
        @show Re
        @show Pr
    end

    pipe_nusselt_number = Nu_packed_bed_Gnielinski!(du, u, cell_id, Re, Pr)

    pipe_convective_heat_transfer_coeff = (pipe_nusselt_number * u.k[cell_id]) / u.catalyst_particle_diameter[cell_id]

    conv_resistance = 1 / pipe_convective_heat_transfer_coeff

    Di = u.pipe_inside_diameter[cell_id]
    Do = u.pipe_outside_diameter[cell_id]

    internal_cond_resistance = (Di * log(Do / Di)) / (2 * u.pipe_k[cell_id])

    overall_heat_transfer_coefficient = 1 / (conv_resistance + internal_cond_resistance)

    return overall_heat_transfer_coefficient
end