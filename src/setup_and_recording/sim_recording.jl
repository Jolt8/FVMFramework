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
            for field in propertynames(u_named[step])
                data = u_named[step][field]
                if field != :mass_fractions
                    write_cell_data(vtk, data, String(field))
                else
                    if data isa ComponentVector
                        for species in propertynames(data)
                            write_cell_data(vtk, data[species], String(species))
                        end
                    elseif ndims(data) == 2 # Matrix of (species, cells) 
                        for i in eachindex(data[:, 1])
                            write_cell_data(vtk, data[i, :], "mass_fraction_$i")
                        end
                    else
                        # Fallback for other structures
                        write_cell_data(vtk, data, String(field))
                    end
                end
            end
            pvd[t] = vtk
        end
    end
    println("hi")
    vtk_save(pvd)
end