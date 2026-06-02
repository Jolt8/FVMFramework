function dry_run_trial()
    dry_run_config = deepcopy(config)
    dry_run_properties = deepcopy(common_properties)

    general_trial_properties = get_empty_trial_properties() 
    dry_run_properties = merge_properties(dry_run_properties, general_trial_properties)

    air_properties = get_air_properties()
    dry_run_properties = merge_properties(dry_run_properties, air_properties)

    region_names = [dry_run_config.regions[i].name for i in eachindex(dry_run_config.regions)]
    regions = Dict(region_names .=> dry_run_config.regions)

    function dry_run_inlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end

    regions["pipe_inlet"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["pipe_inlet"].initial_conditions.mass_fractions = dry_run_properties.inlet_mass_fractions
    regions["pipe_inlet"].region_function = dry_run_inlet!
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, air_properties)
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, general_trial_properties)
    update_region!(dry_run_config, "pipe_inlet")

    regions["silicon_carbide_preheater"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["silicon_carbide_preheater"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, air_properties)
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, general_trial_properties)
    update_region!(dry_run_config, "silicon_carbide_preheater")

    regions["copper_mesh_reformer"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["copper_mesh_reformer"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, air_properties)
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, general_trial_properties)
    update_region!(dry_run_config, "copper_mesh_reformer")

    regions["pipe_outlet"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["pipe_outlet"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, air_properties)
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, general_trial_properties)
    update_region!(dry_run_config, "pipe_outlet")

    regions["steel_pipe_wall"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["steel_pipe_wall"].properties = merge_properties(regions["steel_pipe_wall"].properties, general_trial_properties)
    update_region!(dry_run_config, "steel_pipe_wall")

    regions["thermocouple_and_jacket"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["thermocouple_and_jacket"].properties = merge_properties(regions["thermocouple_and_jacket"].properties, general_trial_properties)
    update_region!(dry_run_config, "thermocouple_and_jacket")

    regions["heating_wire"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["heating_wire"].properties = merge_properties(regions["heating_wire"].properties, general_trial_properties)
    update_region!(dry_run_config, "heating_wire")

    regions["insulation"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["insulation"].properties = merge_properties(regions["insulation"].properties, general_trial_properties)
    update_region!(dry_run_config, "insulation")

    thermocouple_data = get_dry_run_thermocouple_data()

    function heater_1_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * p.heater_weight_1[1]) / n_heater_1_cells
    end

    function heater_2_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * p.heater_weight_2[1]) / n_heater_2_cells
    end

    function heater_3_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * p.heater_weight_3[1]) / n_heater_3_cells
    end

    function heater_4_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * (1 - p.heater_weight_3[1]) * p.heater_weight_4[1]) / n_heater_4_cells
    end

    function heater_5_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * (1 - p.heater_weight_3[1]) * (1 - p.heater_weight_4[1])) / n_heater_5_cells
    end

    return dry_run_config, dry_run_properties, thermocouple_data, heater_1_wattage_per_cell, heater_2_wattage_per_cell, heater_3_wattage_per_cell, heater_4_wattage_per_cell, heater_5_wattage_per_cell
end
