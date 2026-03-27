function sum_mass_flux_face_to_cell!(du, u, cell_id)
    #for mass_face in du.mass_face[cell_id, :]
        #du.mass[cell_id] += mass_face
    #end
    #du.mass[cell_id] += sum(du.mass_face[cell_id, :])
    #du.mass[cell_id] = du.mass_face[cell_id, 1] + du.mass_face[cell_id, 2] + du.mass_face[cell_id, 3] + du.mass_face[cell_id, 4] + du.mass_face[cell_id, 5] + du.mass_face[cell_id, 6]
    #du.mass[cell_id] += sum(du.mass_face[cell_id, :])
    #for i in 1:6
        #du.mass[cell_id] += du.mass_face[cell_id, i]
    #end
    du.mass[cell_id] = sum(view(du.mass_face, cell_id, :))
end