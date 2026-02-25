function sum_mass_flux_face_to_cell!(du, u, cell_id)
    for face_idx in eachindex(du.mass_face[cell_id])
        du.mass[cell_id] += du.mass_face[cell_id][face_idx]
    end
end