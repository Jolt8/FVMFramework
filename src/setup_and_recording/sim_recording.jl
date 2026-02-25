function write_to_vtk_helper!(vtk, u_named_step, n_cells)
    for field in keys(u_named_step)
        data = u_named_step[field]

        if data isa ComponentVector || data isa NamedTuple
            write_to_vtk_helper!(vtk, data, n_cells)
        elseif data isa Vector && length(data) == n_cells 
            #the length check is to handle cases where a field is not present on all cells, which write_to_vtk wouldn't accept
            #this is tricky because we can't exactly write stuff to cells that they're not on
            #one way to handle this is to create an array of length n_cells and then fill it with zeros and then write to the respective cell_ids of the field
            #this doens't work for things that are indexed by controller_id for example, which is what the current logic is handling
            write_cell_data(vtk, data, String(field))
        elseif data isa SubArray
            write_cell_data(vtk, data[1], String(field))
        else
            println("not processed: ", field)
        end
    end
end

function sol_to_vtk(sol, u_named, grid, sim_file)
    date_and_time = Dates.format(now(), "I.MM.SS p yyyy-mm-dd")
    #date_and_time = Dates.format(now(), "I.MM.SS p")

    root_dir = "C://Users//wille//Desktop//Julia_cfd_output_files"

    project_name = replace(basename(sim_file), r".jl" => "")

    sim_folder_name = project_name * " " * date_and_time

    output_dir = joinpath(root_dir, sim_folder_name)

    mkpath(output_dir)

    pvd_filename = joinpath(output_dir, "solution_collection")

    pvd = paraview_collection(pvd_filename)

    step_filename = joinpath(output_dir, "timestep")

    #this step may become a significant bottleneck
    #update: This is a very big bottleneck
    #we get a component array of the simulation's state at this step

    n_cells = length(grid.cells)

    for (step, t) in enumerate(sol.t)
        VTKGridFile(step_filename * " $step" * " at $t.vtu", grid) do vtk
            write_to_vtk_helper!(vtk, u_named[step], n_cells)
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end