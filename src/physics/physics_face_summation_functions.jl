function sum_mass_flux_face_to_cell!(du, u, cell_id)
    du.mass[cell_id] += sum(view(du.mass_face, cell_id, :))
end