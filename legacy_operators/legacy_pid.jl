function pid_temp_controller(du, u, controller_id, monitored_cells, affected_cells, cell_volumes)
        println(u.controllers.desired_value_comp_vector[controller_id])
        field = propertynames(u.controllers.desired_value_comp_vector[controller_id])[1]
        desired = getproperty(u.controllers.desired_value_comp_vector[controller_id], field) #desired will be a single scalar value
        measured_vec = getproperty(u, field)
        measured_du_vec = getproperty(du, field)

        measured_avg = 0.0
        measured_du_avg = 0.0

        @batch for monitored_cell_id in monitored_cells
            measured_avg += measured_vec[monitored_cell_id]
            measured_du_avg += measured_du_vec[monitored_cell_id]
        end

        measured_avg /= length(monitored_cells)
        measured_du_avg /= length(monitored_cells)

        error = measured_avg - desired

        du.integral_error[controller_id] = error

        corrected_volumetric_addition = (
            u.controllers.initial_volumetric_input[controller_id] +
            (u.controllers.proportional_gain[controller_id] * error) +
            (u.controllers.integral_time[controller_id] * u.integral_error[controller_id]) +
            (u.controllers.derivative_time[controller_id] * measured_du_avg)
        )

        corrected_volumetric_addition = clamp(corrected_volumetric_addition, u.controllers.min_volumetric_input[controller_id], u.controllers.max_volumetric_input[controller_id])

        du_field_vec = getproperty(du, field)

        @batch for affected_cell_id in affected_cells
            du_field_vec[affected_cell_id] += corrected_volumetric_addition * cell_volumes[affected_cell_id]
        end
    end