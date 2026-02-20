function sum_mass_flux_face_to_cell!(du, u, cell_id, n_faces)
    for face_idx in 1:n_faces
        du.mass[cell_id] += du.mass_face[cell_id][face_idx]
    end
end