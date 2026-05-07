function return_common_reformer_properties(pipe_length, n_cells, cell_lengths_along_pipe)
    van_t_hoff_A = ComponentVector(CH3O = 1.7e-6u"s^-1", HCOO = 4.74e-13u"s^-1", OH = 3.32e-14u"s^-1")
    van_t_hoff_dH = ComponentVector(CH3O = -46800.0u"J/mol", HCOO = -115000.0u"J/mol", OH = -110000.0u"J/mol")

    MSR_rxn = ComponentVector(
        heat_of_reaction = 49500.0u"J/mol", 
        ref_delta_G = -3800.0u"J/mol", 
        ref_temp = 298.15u"K", 
        kf_A = 1.59e10u"s^-1", #sources online point to values around 1.25e7 mol / (kg * s * bar)
        kf_Ea = 103000.0u"J/mol",
        reactant_stoich_coeffs = (methanol = 1, water = 1),
        product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3), 
        stoich_coeffs = (methanol = -1, water = -1, carbon_dioxide = 1, hydrogen = 3), # CH3OH + H2O ⇋ CO2 + 3(H2)
        van_t_hoff_A = van_t_hoff_A, 
        van_t_hoff_dH = van_t_hoff_dH
    )

    MD_rxn = ComponentVector(
        heat_of_reaction = 90200.0u"J/mol", 
        ref_delta_G = 24800.0u"J/mol", 
        ref_temp = 298.15u"K", 
        kf_A = 1.46e13u"s^-1", #sources online point to values around 1.15e11 mol / (kg * s * bar)
        kf_Ea = 170000.0u"J/mol",
        reactant_stoich_coeffs = (methanol = 1,), 
        product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), 
        stoich_coeffs = (methanol = -1, carbon_monoxide = 1, hydrogen = 2), # CH3OH ⇋ CO + 2(H2)
        van_t_hoff_A = van_t_hoff_A,
        van_t_hoff_dH = van_t_hoff_dH
    )

    WGS_rxn = ComponentVector(
        heat_of_reaction = -41100.0u"J/mol", 
        ref_delta_G = -28600.0u"J/mol", 
        ref_temp = 298.15u"K", 
        kf_A = 4.63e10u"s^-1", #sources online point to values around 3.65e7 mol / (kg * s * bar)
        kf_Ea = 87500.0u"J/mol",
        reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), 
        product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), 
        stoich_coeffs = (carbon_monoxide = -1, water = -1, carbon_dioxide = 1, hydrogen = 1), # CO + H2O ⇋ CO2 + H2
        van_t_hoff_A = van_t_hoff_A, 
        van_t_hoff_dH = van_t_hoff_dH
    )

    empty_mass_fractions = ComponentVector(
        methanol = 1e-20u"kg/kg",
        water = 1e-20u"kg/kg",
        carbon_monoxide = 1e-20u"kg/kg",
        hydrogen = 1e-6u"kg/kg",
        carbon_dioxide = 1e-20u"kg/kg",
        air = 1.0u"kg/kg"
    )

    total_empty_mass_fractions = sum(empty_mass_fractions)
    empty_mass_fractions ./= total_empty_mass_fractions

    pipe_area = pi * (pipe_inside_diameter / 2)^2

    #cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)]
    #this only works because the base of the pipe is at z = 0.0

    TC1_position_along_reactor = 1.5u"inch"
    TC1_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC1_position_along_reactor))

    TC2_position_along_reactor = 3.0u"inch"
    TC2_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC2_position_along_reactor))

    TC3_position_along_reactor = 5.0u"inch"
    TC3_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC3_position_along_reactor))

    TC4_position_along_reactor = 7.5u"inch"
    TC4_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC4_position_along_reactor))

    TC5_position_along_reactor = 10.0u"inch"
    TC5_closest_cell_id = argmin(abs.(cell_lengths_along_pipe .- TC5_position_along_reactor))

    reforming_area_properties = ComponentVector(
        k = 0.026u"W/(m*K)", 
        cp = 4186u"J/(kg*K)",
        dynamic_viscosity = 1.81e-5u"Pa*s",
        rho = 998.2u"kg/m^3",
        R_gas = 8.314u"J/(mol*K)",

        pipe_mass_flow = 0.0u"g/minute",

        pipe_inside_diameter = pipe_inside_diameter,
        pipe_area = pipe_area,
        pipe_length = pipe_length,
        per_cell_pipe_length = pipe_length / n_cells,
        cell_lengths_along_pipe = cell_lengths_along_pipe,

        bed_void_fraction = 0.4,
        packing_surface_area = 100.0u"m^2/m^3",
        particle_diameter = 5.0u"mm",

        overall_heat_transfer_coefficient = 1000.0u"W/(m^2*K)",
        measured_room_temp = 18.0u"°C", #change for each experiment

        inlet_temperature = 80u"°C",
        outlet_temperature = 350.0u"°C",

        TC1_closest_cell_id = TC1_closest_cell_id,
        TC2_closest_cell_id = TC2_closest_cell_id,
        TC3_closest_cell_id = TC3_closest_cell_id,
        TC4_closest_cell_id = TC4_closest_cell_id,
        TC5_closest_cell_id = TC5_closest_cell_id,

        external_temp = 300.0u"°C",
        saturation_temp = 72.4u"°C",
        liquid_rho = 791.0u"kg/m^3",
        gas_rho = 1.225u"kg/m^3",
        mass_transfer_coeff_vap = 0.001u"kg/(m^2*s*K)",
        heat_of_vaporization = 1.5u"kJ/g",
        liquid_feed_mass_fractions = empty_mass_fractions,

        diffusion_coefficients = (
            methanol = 1e-5u"m^2/s",
            water = 1e-5u"m^2/s",
            carbon_monoxide = 1e-5u"m^2/s",
            hydrogen = 1e-5u"m^2/s",
            carbon_dioxide = 1e-5u"m^2/s",
            air = 1e-5u"m^2/s"
        ), 
        molecular_weights = (
            methanol = 32.04u"g/mol",
            water = 18.02u"g/mol",
            carbon_monoxide = 28.01u"g/mol",
            hydrogen = 2.02u"g/mol",
            carbon_dioxide = 44.01u"g/mol",
            air = 28.97u"g/mol"
        ), 
        reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
        reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0u"kg/m^3", MD_rxn = 1250.0u"kg/m^3", WGS_rxn = 1250.0u"kg/m^3"),), 
    )

    permeability = (reforming_area_properties.particle_diameter^2 * reforming_area_properties.bed_void_fraction^3) / (150.0 * (1.0 - reforming_area_properties.bed_void_fraction)^2)

    reforming_area_properties = merge_properties(reforming_area_properties, ComponentVector(permeability = permeability))

    return reforming_area_properties
end

