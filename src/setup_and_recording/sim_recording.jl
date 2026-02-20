function write_to_vtk_helper!(vtk, u_named_step)
    for field in propertynames(u_named_step)
        data = u_named_step[field]

        if data isa ComponentVector
            write_to_vtk_helper!(vtk, data)
        elseif data isa Vector || data isa SubArray
            write_cell_data(vtk, data, String(field))
        else
            println("not_processed: ", field)
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

    for (step, t) in enumerate(sol.t)
        VTKGridFile(step_filename * " $step" * " at $t.vtu", grid) do vtk
            write_to_vtk_helper!(vtk, u_named[step])
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end
