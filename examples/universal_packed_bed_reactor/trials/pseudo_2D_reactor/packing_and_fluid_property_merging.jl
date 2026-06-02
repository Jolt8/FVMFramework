

function get_fluid_and_copper_mesh_packing_rho!(du, u, cell_id, vol)
    u.rho[cell_id] = (u.copper_rho[cell_id] * (1 - u.bed_void_fraction[cell_id])) + (u.fluid_rho[cell_id] * u.bed_void_fraction[cell_id])
end

function get_fluid_and_silicon_carbide_bed_packing_rho!(du, u, cell_id, vol)
    u.rho[cell_id] = (u.silicon_carbide_rho[cell_id] * (1 - u.bed_void_fraction[cell_id])) + (u.fluid_rho[cell_id] * u.bed_void_fraction[cell_id])
end

function get_fluid_and_copper_mesh_packing_k!(du, u, cell_id, vol)
    u.k[cell_id] = (1 - u.bed_void_fraction[cell_id]) * u.copper_k[cell_id] + u.bed_void_fraction[cell_id] * u.fluid_k[cell_id]
end

function get_fluid_and_silicon_carbide_bed_packing_k!(du, u, cell_id, vol)
    u.k[cell_id] = (1 - u.bed_void_fraction[cell_id]) * u.silicon_carbide_k[cell_id] + u.bed_void_fraction[cell_id] * u.fluid_k[cell_id]
end

function get_fluid_and_copper_mesh_packing_cp!(du, u, cell_id, vol)
    u.cp[cell_id] = ((1 - u.bed_void_fraction[cell_id]) * u.copper_rho[cell_id] * u.copper_cp[cell_id] + u.bed_void_fraction[cell_id] * u.fluid_rho[cell_id] * u.fluid_cp[cell_id]) / u.rho[cell_id]
end

function get_fluid_and_silicon_carbide_bed_packing_cp!(du, u, cell_id, vol)
    u.cp[cell_id] = ((1 - u.bed_void_fraction[cell_id]) * u.silicon_carbide_rho[cell_id] * u.silicon_carbide_cp[cell_id] + u.bed_void_fraction[cell_id] * u.fluid_rho[cell_id] * u.fluid_cp[cell_id]) / u.rho[cell_id]
end
