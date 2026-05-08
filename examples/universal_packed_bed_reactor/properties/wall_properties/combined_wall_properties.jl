function return_combined_wall_properties(pipe_length, n_cells, cell_lengths_along_pipe, wall_lengths_along_pipe)
    #we use wall_lenghts_along_pipe here because the thermocouple is in the walls
    TC1_position_along_reactor = 1.5u"inch"
    TC1_closest_cell_id = argmin(abs.(wall_lengths_along_pipe .- TC1_position_along_reactor))

    TC2_position_along_reactor = 3.0u"inch"
    TC2_closest_cell_id = argmin(abs.(wall_lengths_along_pipe .- TC2_position_along_reactor))

    TC3_position_along_reactor = 5.0u"inch"
    TC3_closest_cell_id = argmin(abs.(wall_lengths_along_pipe .- TC3_position_along_reactor))

    TC4_position_along_reactor = 7.5u"inch"
    TC4_closest_cell_id = argmin(abs.(wall_lengths_along_pipe .- TC4_position_along_reactor))

    TC5_position_along_reactor = 10.0u"inch"
    TC5_closest_cell_id = argmin(abs.(wall_lengths_along_pipe .- TC5_position_along_reactor))

    reforming_area_properties = ComponentVector(
        k = 0.026u"W/(m*K)", 
        cp = 4186u"J/(kg*K)",
        rho = 998.2u"kg/m^3",
        R_gas = 8.314u"J/(mol*K)",

        wall_inside_diameter = 0.5u"inch" |> u"m",
        wall_area = pi * (wall_inside_diameter / 2)^2,
        wall_thickness = (0.625u"inch" - 0.5u"inch") |> u"m",
        
        wall_length = pipe_length,
        per_cell_wall_length = pipe_length / n_cells,
        wall_lengths_along_pipe = wall_lengths_along_pipe,

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

