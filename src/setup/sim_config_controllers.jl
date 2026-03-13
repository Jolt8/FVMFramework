
function add_controller!(
    config, name;
    controller,
    monitored_cellset,
    affected_cellset,
    controller_function
)
    monitored_cells = collect(getcellset(config.grid, monitored_cellset))
    affected_cells = collect(getcellset(config.grid, affected_cellset))

    controller = ControllerSetupInfo(name, controller, monitored_cellset, affected_cellset, controller_function, monitored_cells, affected_cells)

    if controller in config.controllers
        existing_controller_idx = findfirst(x -> x == controller, config.controllers)

        config.controllers[existing_controller_idx] = controller
    else
        push!(config.controllers, controller)
    end
    return
end