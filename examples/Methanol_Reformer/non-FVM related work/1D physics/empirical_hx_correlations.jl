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

function calculate_UA_prime_system(du, u, cell_id, vol)
    Re_1 = (u.rho[cell_id] * u.velocity[cell_id] * u.channel_hydraulic_diameter[cell_id]) / (u.viscosity[cell_id])
    Pr_1 = (u.viscosity[cell_id] * u.cp[cell_id]) / (u.pipe_k[cell_id])

    Re2f_1 = u.Re2f_interpolation
    if Re_1 < 2200
        fd_1 = 64/Re_1 # Laminar approximation
    else
        fd_1 = Re2f_1(Re_1)
    end

    # -- Pipe 2 -- 
    fluid_2_thermal_conductivity = 0.034u"W/(m*K)"

    fluid_2_molar_density = GasProperties.get_mixture_molar_density(model_2_T, model_2_P)
    fluid_2_density = GasProperties.get_mixture_density(model_2.coolprop_components, model_2_common.z_vec, model_2_T, model_2_P)

    #println("model_1_T, ", model_1_T)
    println("model_2_T, ", model_2_T)
    #println("model_2_common.z_vec, ", model_2_common.z_vec)
    fluid_2_specific_heat = GasProperties.get_mixture_cp_mass(model_2.coolprop_components, model_2_common.z_vec, model_2_T)

    fluid_2_total_molar_flow = sum(model_2_molar_flows)

    fluid_2_total_volumetric_flow = fluid_2_total_molar_flow / fluid_2_molar_density

    fluid_2_velocity = fluid_2_total_volumetric_flow / model_2.pressure_model.channel_cross_sectional_area

    fluid_2_dynamic_viscosity = 3.5e-5u"Pa*s" #just using something close

    Re_2 = (fluid_2_density * fluid_2_velocity * model_2.pressure_model.channel_hydraulic_diameter) / (fluid_2_dynamic_viscosity)
    #println("fluid_2_dynamic_viscosity ", fluid_2_dynamic_viscosity)
    #println("fluid_2_specific_heat ", fluid_2_specific_heat)
    #println("fluid_2_thermal_conductivity ", fluid_2_thermal_conductivity)
    Pr_2 = (fluid_2_dynamic_viscosity * fluid_2_specific_heat) / (fluid_2_thermal_conductivity)
    #println("Pr_2 ", Pr_2)

    Re2f_2 = model_2.pressure_model.Re2f_interpolation
    if Re_2 < 2200
        fd_2 = 64/Re_2 # Laminar approximation
    else
        fd_2 = Re2f_2(Re_2)
    end

    phys_diameter = model_1.pressure_model.channel_hydraulic_diameter * (pi + 2) / pi

    return get_linear_heat_transfer_coefficient(
        z_position, 
        model_1.pressure_model.channel_hydraulic_diameter, fluid_1_thermal_conductivity, Re_1, Pr_1, fd_1, 
        model_2.pressure_model.channel_hydraulic_diameter, fluid_2_thermal_conductivity, Re_2, Pr_2, fd_2,
        phys_diameter,
        layers
    )
end