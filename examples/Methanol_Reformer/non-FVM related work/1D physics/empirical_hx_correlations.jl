function Nu_packed_bed_Gnielinski!(du, u, cell_id, vol)
    #borrowed from ht.py
    #Calculates Nusselt number of a fluid passing over a bed of particles

    Re = u.rho * u.viscosity * u.catalyst_particle_diameter / u.mu / u.bed_void_fraction
    
    Nu_lam = 0.664 * Re^0.5 * u.Pr^(1/3.)
    
    Nu_turb = 0.037 * Re^0.8 * u.Pr/(1 + 2.443 * Re^-0.1 * (u.Pr^(2/3.) - 1))
    
    Nu_sphere = 2 + (Nu_lam^2 + Nu_turb^2)^0.5
    
    fa = 1.0 + 1.5 * (1.0 - u.bed_void_fraction)
    
    u.Nu = Nu_sphere * fa
end