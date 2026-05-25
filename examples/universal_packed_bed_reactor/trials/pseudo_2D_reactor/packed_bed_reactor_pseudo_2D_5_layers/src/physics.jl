struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end
struct ThermocoupleAndJacket <: AbstractPhysics end
struct HeatingWire <: AbstractPhysics end

function update_copper_mesh_properties!(du, u, cell_id, vol, system)
    mw_avg!(u, cell_id)
    get_fluid_and_copper_mesh_packing_rho!(du, u, cell_id, vol)
    get_fluid_and_copper_mesh_packing_k!(du, u, cell_id, vol)
    get_fluid_and_copper_mesh_packing_cp!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
end

function update_silicon_carbide_bed_packing_properties!(du, u, cell_id, vol, system)
    mw_avg!(u, cell_id)
    get_fluid_and_silicon_carbide_bed_packing_rho!(du, u, cell_id, vol)
    get_fluid_and_silicon_carbide_bed_packing_k!(du, u, cell_id, vol)
    get_fluid_and_silicon_carbide_bed_packing_cp!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
end

function update_fluid_properties!(du, u, cell_id, vol, system)
    mw_avg!(u, cell_id)

    properties = ComponentVector(system.properties_vec, system.properties_axes)
    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]

    molar_concentrations!(u, cell_id)
end

function update_solid_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]
end

function update_steel_pipe_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id] * u.steel_thermal_mass_multiplier[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]
end

function fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    sum_mass_flux_face_to_cell!(du, u, cell_id) # this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
end

function solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

# Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    all_species_advection!(
        du, u, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
    
    enthalpy_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

function fluid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    fluid_to_steel_pipe_convective_heat_transfer!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

function solid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

function no_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    # we do nothing here because we're going to apply thermal resistance based model somewhere else
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_solid_flux!

    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return fluid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!

    typeof(phys_a) <: ThermocoupleAndJacket && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: ThermocoupleAndJacket && return solid_solid_flux!
    typeof(phys_a) <: ThermocoupleAndJacket && typeof(phys_b) <: ThermocoupleAndJacket && return no_flux! 
    # this is made almost completely of fiberglass tape, so axial heat transfer along the pipe between thermocouple cells can be ignored
    
    typeof(phys_a) <: HeatingWire && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: HeatingWire && return solid_solid_flux!
    typeof(phys_a) <: HeatingWire && typeof(phys_b) <: HeatingWire && return no_flux!
    # this is made of heating wire covered with ceramic fiber, so axial heat transfer along the pipe between thermocouple cells can be ignored
    
    typeof(phys_a) <: HeatingWire && typeof(phys_b) <: ThermocoupleAndJacket && return no_flux! # we apply custom thermal resistances between these two, so no flux should occur
    typeof(phys_a) <: ThermocoupleAndJacket && typeof(phys_b) <: HeatingWire && return no_flux!
end
