fluid_regions = ["pipe_inlet", "silicon_carbide_preheater", "copper_mesh_reformer", "pipe_outlet"]
advecting_fluid_cells = vcat(collect(grid.cellsets["pipe_inlet"]), collect(grid.cellsets["silicon_carbide_preheater"]), collect(grid.cellsets["copper_mesh_reformer"]), collect(grid.cellsets["pipe_outlet"]))

function trial_independent_solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, system.p_axes)

    for cell_id in eachindex(geo.cell_volumes)
        u.insulation_to_air_overall_heat_transfer_coefficient_to_environment[cell_id] = p.insulation_to_air_overall_heat_transfer_coefficient_to_environment[1]
        u.pipe_endcaps_to_air_thermal_conductance[cell_id] = p.pipe_endcaps_to_air_thermal_conductance[1]
        u.fluid_to_steel_pipe_convective_heat_transfer_coefficient[cell_id] = p.fluid_to_steel_pipe_convective_heat_transfer_coefficient[1]
        u.steel_thermal_mass_multiplier[cell_id] = p.steel_thermal_mass_multiplier[1]
    end

    # VERY IMPORTANT: since most software uses 0-based indexing, you need to adjust the cell id by +1
    # for example, if you mouse over cell_id 5161 in ParaView, you need to use 5162 in the code because julia uses 1-based indexing 

    # we could probably just add a pre_calculations function for this in each region group, but this works for now
    for reg in system.region_groups
        if reg.name == "pipe_inlet" || reg.name == "silicon_carbide_preheater"
            for cell_id in reg.region_cells
                update_silicon_carbide_bed_packing_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        elseif reg.name == "copper_mesh_reformer" || reg.name == "pipe_outlet"
            for cell_id in reg.region_cells
                update_copper_mesh_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        elseif reg.name in fluid_regions
            for cell_id in reg.region_cells
                update_fluid_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        elseif reg.name == "steel_pipe_wall"
            for cell_id in reg.region_cells
                update_steel_pipe_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        else
            for cell_id in reg.region_cells
                update_solid_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
            end
        end
    end

    for i in 1:length(advecting_fluid_cells) - 1 # we don't take mass out of the outlet
        idx_a = advecting_fluid_cells[i]
        idx_b = advecting_fluid_cells[i + 1]
        
        du.mass_face[idx_a, 6] -= u.pipe_mass_flow[idx_a]
        du.mass_face[idx_b, 1] += u.pipe_mass_flow[idx_a]
    end

    for cell_id in TC1_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC1_thermal_resistance[1]
    end
    for cell_id in TC2_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC2_thermal_resistance[1]
    end
    for cell_id in TC3_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC3_thermal_resistance[1]
    end
    for cell_id in TC4_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC4_thermal_resistance[1]
    end
    for cell_id in TC5_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC5_thermal_resistance[1]
    end

    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system)
end

function dry_run_solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, p_axes)

    for cell_id in heater_1_cells
        du.heat[cell_id] += dry_run_heater_1_wattage_per_cell(p, t)
    end

    for cell_id in heater_2_cells
        du.heat[cell_id] += dry_run_heater_2_wattage_per_cell(p, t)
    end

    for cell_id in heater_3_cells
        du.heat[cell_id] += dry_run_heater_3_wattage_per_cell(p, t)
    end

    for cell_id in heater_4_cells
        du.heat[cell_id] += dry_run_heater_4_wattage_per_cell(p, t)
    end

    for cell_id in heater_5_cells
        du.heat[cell_id] += dry_run_heater_5_wattage_per_cell(p, t)
    end

    trial_independent_solve_system!(du, u, p_vec, t, geo, system)
end

function hot_water_solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, p_axes)

    for cell_id in advecting_fluid_cells
        hot_water_pump_shut_off(du, u, cell_id, t)
    end

    for cell_id in heater_1_cells
        du.heat[cell_id] += hot_water_heater_1_wattage_per_cell(p, t)
    end

    for cell_id in heater_2_cells
        du.heat[cell_id] += hot_water_heater_2_wattage_per_cell(p, t)
    end

    for cell_id in heater_3_cells
        du.heat[cell_id] += hot_water_heater_3_wattage_per_cell(p, t)
    end

    for cell_id in heater_4_cells
        du.heat[cell_id] += hot_water_heater_4_wattage_per_cell(p, t)
    end

    for cell_id in heater_5_cells
        du.heat[cell_id] += hot_water_heater_5_wattage_per_cell(p, t)
    end
    
    trial_independent_solve_system!(du, u, p_vec, t, geo, system)
end

function built_trial_implicit_prob(f_closure, du0_vec, u0_vec, thermocouple_data, p_guess)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tMax = ustrip(upreferred(thermocouple_data.timestamps[end]))
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end
