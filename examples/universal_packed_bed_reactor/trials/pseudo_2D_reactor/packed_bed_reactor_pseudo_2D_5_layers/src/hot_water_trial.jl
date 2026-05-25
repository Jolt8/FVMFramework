function hot_water_trial()
    hot_water_config = deepcopy(config)
    hot_water_trial_properties = deepcopy(common_properties)

    values_of_note = get_packed_bed_reactor_water_flow_trial_values_of_note()
    pump_shutoff_timestamp = ustrip(upreferred(values_of_note.pump_shut_off_time))

    general_trial_properties = get_hot_water_trial_properties(values_of_note)
    hot_water_trial_properties = merge_properties(hot_water_trial_properties, general_trial_properties)

    water_properties = get_water_properties()
    hot_water_trial_properties = merge_properties(hot_water_trial_properties, water_properties)

    inlet_and_outlet_temperatures = get_inlet_and_outlet_temperature_correlations()
    inlet_temp_interp = inlet_and_outlet_temperatures.inlet_temp_interp

    region_names = [hot_water_config.regions[i].name for i in eachindex(hot_water_config.regions)]
    regions = Dict(region_names .=> hot_water_config.regions)

    function hot_water_inlet!(du, u, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]
        
        du.heat[cell_id] *= 0.0

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end

    regions["pipe_inlet"].initial_conditions.temp = inlet_temp_interp(0.0)
    regions["pipe_inlet"].initial_conditions.mass_fractions = hot_water_trial_properties.inlet_mass_fractions
    regions["pipe_inlet"].region_function = hot_water_inlet!
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, water_properties)
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, general_trial_properties)
    update_region!(hot_water_config, "pipe_inlet")

    regions["silicon_carbide_preheater"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["silicon_carbide_preheater"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, water_properties)
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, general_trial_properties)
    update_region!(hot_water_config, "silicon_carbide_preheater")

    regions["copper_mesh_reformer"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["copper_mesh_reformer"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, water_properties)
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, general_trial_properties)
    update_region!(hot_water_config, "copper_mesh_reformer")

    regions["pipe_outlet"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["pipe_outlet"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, water_properties)
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, general_trial_properties)
    update_region!(hot_water_config, "pipe_outlet")

    regions["steel_pipe_wall"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["steel_pipe_wall"].properties = merge_properties(regions["steel_pipe_wall"].properties, general_trial_properties)
    update_region!(hot_water_config, "steel_pipe_wall")

    regions["thermocouple_and_jacket"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["thermocouple_and_jacket"].properties = merge_properties(regions["thermocouple_and_jacket"].properties, general_trial_properties)
    update_region!(hot_water_config, "thermocouple_and_jacket")

    regions["heating_wire"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["heating_wire"].properties = merge_properties(regions["heating_wire"].properties, general_trial_properties)
    update_region!(hot_water_config, "heating_wire")

    regions["insulation"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["insulation"].properties = merge_properties(regions["insulation"].properties, general_trial_properties)
    update_region!(hot_water_config, "insulation")

    inlet_cell_id = collect(grid.cellsets["pipe_inlet"])[1]

    function pump_shut_off(du, u, cell_id, t)
        if t <= pump_shutoff_timestamp # pump on
            du.temp[inlet_cell_id] *= 0.0
            du.temp[inlet_cell_id] += DataInterpolations.derivative(inlet_temp_interp, ForwardDiff.value(t))

            u.pipe_mass_flow[cell_id] = ustrip(upreferred(hot_water_trial_properties.pipe_mass_flow))
        else # pump shut off
            u.pipe_mass_flow[cell_id] = 0.0
        end
    end

    thermocouple_data = get_hot_water_thermocouple_data()

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

    return hot_water_config, hot_water_trial_properties, pump_shut_off, thermocouple_data, heater_1_wattage_per_cell, heater_2_wattage_per_cell, heater_3_wattage_per_cell, heater_4_wattage_per_cell, heater_5_wattage_per_cell
end
