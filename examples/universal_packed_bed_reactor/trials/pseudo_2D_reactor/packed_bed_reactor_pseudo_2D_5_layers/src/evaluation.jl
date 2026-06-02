function get_trial_temperature_variance(thermocouple_data)
    all_recorded_temperatures = []

    for t in ustrip.(thermocouple_data.timestamps)
        push!(all_recorded_temperatures, thermocouple_data.TC1_temps_interp(t))
        push!(all_recorded_temperatures, thermocouple_data.TC2_temps_interp(t))
        push!(all_recorded_temperatures, thermocouple_data.TC3_temps_interp(t))
        push!(all_recorded_temperatures, thermocouple_data.TC4_temps_interp(t))
        push!(all_recorded_temperatures, thermocouple_data.TC5_temps_interp(t))
    end

    mean_recorded_temperatures = sum(all_recorded_temperatures) / length(all_recorded_temperatures)

    sum_of_squares = 0.0

    for i in eachindex(all_recorded_temperatures)
        sum_of_squares += (all_recorded_temperatures[i] - mean_recorded_temperatures)^2
    end

    return sum_of_squares / length(all_recorded_temperatures)
end

dry_run_temperature_variance = get_trial_temperature_variance(dry_run_thermocouple_data)
hot_water_temperature_variance = get_trial_temperature_variance(hot_water_thermocouple_data)

function dry_run_loss(θ)
    prob = get!(task_local_storage(), :dry_run_implicit_prob) do
        # Deepcopy the mutable state and geometry objects
        system_copy = deepcopy(dry_run_system)
        geo_copy = deepcopy(dry_run_geo)
        
        # Build a new closure bound to the thread-isolated copies
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, dry_run_solve_system!, geo_copy, system_copy)
        
        # Re-generate the Jacobian sparsity and the implicit ODEProblem
        built_trial_implicit_prob(f_closure, dry_run_du0_vec, dry_run_u0_vec, dry_run_thermocouple_data, θ)
    end
    
    dry_run_loss_prob = remake(prob, p = θ)

    dry_run_sol = solve(
        dry_run_loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        saveat = dry_run_save_interval,
        tstops = dry_run_tstops,
    )

    if length(dry_run_sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(dry_run_sol.u[i], dry_run_state_axes) for i in eachindex(dry_run_sol.u)]

    mean_squared_error = 0.0

    for i in eachindex(dry_run_sol.t)
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC1_cell_id]) - ustrip(dry_run_thermocouple_data.TC1_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC2_cell_id]) - ustrip(dry_run_thermocouple_data.TC2_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC3_cell_id]) - ustrip(dry_run_thermocouple_data.TC3_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC4_cell_id]) - ustrip(dry_run_thermocouple_data.TC4_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC5_cell_id]) - ustrip(dry_run_thermocouple_data.TC5_temps_interp(dry_run_sol.t[i])))
    end

    loss = (mean_squared_error / dry_run_temperature_variance) / length(dry_run_sol.t)

    return dry_run_sol.t, u_named, loss
end

function viewable_dry_run_loss(θ)
    sol_t, u_named, loss = dry_run_loss(θ)

    simulated_thermocouple_data = (
        dry_run_times = sol_t,
        dry_run_TC1 = [u_named[i].temp[TC1_cell_id] for i in eachindex(sol_t)],
        dry_run_TC2 = [u_named[i].temp[TC2_cell_id] for i in eachindex(sol_t)],
        dry_run_TC3 = [u_named[i].temp[TC3_cell_id] for i in eachindex(sol_t)],
        dry_run_TC4 = [u_named[i].temp[TC4_cell_id] for i in eachindex(sol_t)],
        dry_run_TC5 = [u_named[i].temp[TC5_cell_id] for i in eachindex(sol_t)],
        loss = loss
    )

    return simulated_thermocouple_data
end

function pure_dry_run_loss(θ)
    sol_t, u_named, loss = dry_run_loss(θ)

    return loss 
end

function hot_water_loss(θ)
    prob = get!(task_local_storage(), :hot_water_implicit_prob) do
        # Deepcopy the mutable state and geometry objects
        system_copy = deepcopy(hot_water_system)
        geo_copy = deepcopy(hot_water_geo)
        
        # Build a new closure bound to the thread-isolated copies
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, hot_water_solve_system!, geo_copy, system_copy)
        
        # Re-generate the Jacobian sparsity and the implicit ODEProblem
        built_trial_implicit_prob(f_closure, hot_water_du0_vec, hot_water_u0_vec, hot_water_thermocouple_data, θ)
    end
    
    hot_water_loss_prob = remake(prob, p = θ)

    hot_water_sol = solve(
        hot_water_loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        saveat = hot_water_save_interval
    )

    # this is to prevent the optimizer from crashing the simulation to get a lower mean_squared_error
    if length(hot_water_sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(hot_water_sol.u[i], hot_water_state_axes) for i in eachindex(hot_water_sol.u)]

    mean_squared_error = 0.0

    for i in eachindex(hot_water_sol.t)
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC1_cell_id]) - ustrip(hot_water_thermocouple_data.TC1_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC2_cell_id]) - ustrip(hot_water_thermocouple_data.TC2_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC3_cell_id]) - ustrip(hot_water_thermocouple_data.TC3_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC4_cell_id]) - ustrip(hot_water_thermocouple_data.TC4_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC5_cell_id]) - ustrip(hot_water_thermocouple_data.TC5_temps_interp(hot_water_sol.t[i])))
    end

    loss = (mean_squared_error / hot_water_temperature_variance) / length(hot_water_sol.t)

    return hot_water_sol.t, u_named, loss
end

function viewable_hot_water_loss(θ)
    sol_t, u_named, loss = hot_water_loss(θ)

    simulated_thermocouple_data = (
        hot_water_times = sol_t,
        hot_water_TC1 = [u_named[i].temp[TC1_cell_id] for i in eachindex(sol_t)],
        hot_water_TC2 = [u_named[i].temp[TC2_cell_id] for i in eachindex(sol_t)],
        hot_water_TC3 = [u_named[i].temp[TC3_cell_id] for i in eachindex(sol_t)],
        hot_water_TC4 = [u_named[i].temp[TC4_cell_id] for i in eachindex(sol_t)],
        hot_water_TC5 = [u_named[i].temp[TC5_cell_id] for i in eachindex(sol_t)],
        loss = loss
    )
    
    return simulated_thermocouple_data
end

function pure_hot_water_loss(θ)
    sol_t, u_named, loss = hot_water_loss(θ)

    return loss
end

function loss(θ)
    return (pure_dry_run_loss(θ) + pure_hot_water_loss(θ)) / 2.0
end

plot_output_dir = joinpath(@__DIR__, "..", "..", "graphs")

function generate_validation_plots(dry_run_p_best, hot_water_p_best)
    # Ensure plots directory exists
    mkpath(plot_output_dir)

    # 1. Dry Run Plots
    dry_run_simulated = viewable_dry_run_loss(dry_run_p_best)
    dry_run_times = dry_run_simulated.dry_run_times

    TC1_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC1, label="Sim TC1", linewidth=2)
    plot!(TC1_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC1_temps_interp.(dry_run_times)), label="Exp TC1", linewidth=2)
    savefig(TC1_plt, joinpath(plot_output_dir, "TC1_dry_run_loss.png"))

    TC2_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC2, label="Sim TC2", linewidth=2)
    plot!(TC2_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC2_temps_interp.(dry_run_times)), label="Exp TC2", linewidth=2)
    savefig(TC2_plt, joinpath(plot_output_dir, "TC2_dry_run_loss.png"))

    TC3_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC3, label="Sim TC3", linewidth=2)
    plot!(TC3_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC3_temps_interp.(dry_run_times)), label="Exp TC3", linewidth=2)
    savefig(TC3_plt, joinpath(plot_output_dir, "TC3_dry_run_loss.png"))

    TC4_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC4, label="Sim TC4", linewidth=2)
    plot!(TC4_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC4_temps_interp.(dry_run_times)), label="Exp TC4", linewidth=2)
    savefig(TC4_plt, joinpath(plot_output_dir, "TC4_dry_run_loss.png"))

    TC5_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC5, label="Sim TC5", linewidth=2)
    plot!(TC5_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC5_temps_interp.(dry_run_times)), label="Exp TC5", linewidth=2)
    savefig(TC5_plt, joinpath(plot_output_dir, "TC5_dry_run_loss.png"))

    overall_plot = plot(TC1_plt, TC2_plt, TC3_plt, TC4_plt, TC5_plt, layout=(5, 1), size=(1000, 250*5), ylims=(250, 600))
    savefig(overall_plot, joinpath(plot_output_dir, "overall_dry_run_loss.png"))

    # 2. Hot Water Plots
    hot_water_simulated = viewable_hot_water_loss(hot_water_p_best)
    hot_water_times = hot_water_simulated.hot_water_times

    TC1_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC1, label="Sim TC1", linewidth=2)
    plot!(TC1_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC1_temps_interp.(hot_water_times)), label="Exp TC1", linewidth=2)
    savefig(TC1_hw_plt, joinpath(plot_output_dir, "TC1_hot_water_loss.png"))

    TC2_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC2, label="Sim TC2", linewidth=2)
    plot!(TC2_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC2_temps_interp.(hot_water_times)), label="Exp TC2", linewidth=2)
    savefig(TC2_hw_plt, joinpath(plot_output_dir, "TC2_hot_water_loss.png"))

    TC3_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC3, label="Sim TC3", linewidth=2)
    plot!(TC3_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC3_temps_interp.(hot_water_times)), label="Exp TC3", linewidth=2)
    savefig(TC3_hw_plt, joinpath(plot_output_dir, "TC3_hot_water_loss.png"))

    TC4_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC4, label="Sim TC4", linewidth=2)
    plot!(TC4_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC4_temps_interp.(hot_water_times)), label="Exp TC4", linewidth=2)
    savefig(TC4_hw_plt, joinpath(plot_output_dir, "TC4_hot_water_loss.png"))

    TC5_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC5, label="Sim TC5", linewidth=2)
    plot!(TC5_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC5_temps_interp.(hot_water_times)), label="Exp TC5", linewidth=2)
    savefig(TC5_hw_plt, joinpath(plot_output_dir, "TC5_hot_water_loss.png"))

    overall_hw_plot = plot(TC1_hw_plt, TC2_hw_plt, TC3_hw_plt, TC4_hw_plt, TC5_hw_plt, layout=(5, 1), size=(1000, 250*5), ylims=(290, 340))
    savefig(overall_hw_plot, joinpath(plot_output_dir, "overall_hot_water_loss.png"))
end
