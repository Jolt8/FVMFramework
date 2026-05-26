using LinearAlgebra

function update_drift_flux_cell_velocities!(u, cell_id::Int)
    gas_phase_holdup_clamped = clamp(u.gas_holdup[cell_id], 0.0, 0.9999)
    u.gas_velocity[cell_id] = u.distribution_parameter[cell_id] * u.mixture_velocity[cell_id] + u.local_gas_drift_velocity[cell_id]
    u.liquid_velocity[cell_id] = (u.mixture_velocity[cell_id] - gas_phase_holdup_clamped * u.gas_velocity[cell_id]) / (1.0 - gas_phase_holdup_clamped)
end

function drift_flux_mass_conservation!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    # Average bed porosity and other properties
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])
    
    mixture_viscosity_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_viscosity[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_viscosity[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_viscosity[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_viscosity[idx_b])
    )
    
    mixture_density_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_density[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_density[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_density[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_density[idx_b])
    )
    
    solid_particle_diameter_average = 0.5 * (u.solid_particle_diameter[idx_a] + u.solid_particle_diameter[idx_b])
    vessel_hydraulic_diameter_average = 0.5 * (u.vessel_hydraulic_diameter[idx_a] + u.vessel_hydraulic_diameter[idx_b])
    
    # Calculate viscous and inertial drag coefficients on the face
    ergun_viscous = 150.0 * (1.0 - bed_porosity_average)^2 / (solid_particle_diameter_average^2 * bed_porosity_average^3) * mixture_viscosity_average
    ergun_inertial = 1.75 * (1.0 - bed_porosity_average) / (solid_particle_diameter_average * bed_porosity_average^3) * mixture_density_average
    
    viscous_drag = (bed_porosity_average < 1.0 ? ergun_viscous : 0.0) + 32.0 * mixture_viscosity_average / (vessel_hydraulic_diameter_average^2)
    inertial_drag = (bed_porosity_average < 1.0 ? ergun_inertial : 0.0) + 2.0 * 0.02 * mixture_density_average / vessel_hydraulic_diameter_average
    
    driving_force = (u.pressure[idx_a] - u.pressure[idx_b]) / dist + mixture_density_average * LinearAlgebra.dot(Ferrite.Vec{3}((0.0, 0.0, -9.81)), norm)
    
    if inertial_drag > 1e-12
        superficial_velocity_magnitude = (-viscous_drag + sqrt(viscous_drag^2 + 4.0 * inertial_drag * abs(driving_force))) / (2.0 * inertial_drag)
    else
        superficial_velocity_magnitude = abs(driving_force) / viscous_drag
    end
    mixture_superficial_velocity = superficial_velocity_magnitude * sign(driving_force)
    mixture_pore_velocity = mixture_superficial_velocity / bed_porosity_average
    
    gas_velocity_face = 0.5 * (u.distribution_parameter[idx_a] + u.distribution_parameter[idx_b]) * mixture_pore_velocity + 0.5 * (u.local_gas_drift_velocity[idx_a] + u.local_gas_drift_velocity[idx_b])
    gas_holdup_clamped = clamp(0.5 * (u.gas_holdup[idx_a] + u.gas_holdup[idx_b]), 0.0, 0.9999)
    liquid_velocity_face = (mixture_pore_velocity - gas_holdup_clamped * gas_velocity_face) / (1.0 - gas_holdup_clamped)
    
    # Upwind phase states
    upwinded_gas_density = gas_velocity_face >= 0.0 ? u.gas_density[idx_a] : u.gas_density[idx_b]
    upwinded_gas_holdup = gas_velocity_face >= 0.0 ? u.gas_holdup[idx_a] : u.gas_holdup[idx_b]
    
    upwinded_liquid_density = liquid_velocity_face >= 0.0 ? u.liquid_density[idx_a] : u.liquid_density[idx_b]
    upwinded_liquid_holdup = liquid_velocity_face >= 0.0 ? (1.0 - u.gas_holdup[idx_a]) : (1.0 - u.gas_holdup[idx_b])
    
    gas_flux_density = bed_porosity_average * upwinded_gas_density * upwinded_gas_holdup * gas_velocity_face
    liquid_flux_density = bed_porosity_average * upwinded_liquid_density * upwinded_liquid_holdup * liquid_velocity_face
    
    # Loop over gas species
    for_fields!(du.species_gas_mass_per_volume, du.volumetric_gas_mass_generation, u.gas_mass_fractions) do species, gas_accumulation, gas_generation, gas_mass_fractions
        upwinded_fraction = gas_velocity_face >= 0.0 ? gas_mass_fractions[species[idx_a]] : gas_mass_fractions[species[idx_b]]
        change_in_volume = - (gas_flux_density * upwinded_fraction * area / vol_a)
        if face_idx == 1
            change_in_volume += gas_generation[species[idx_a]]
        end
        gas_accumulation[species[idx_a]] += change_in_volume
    end
    
    # Loop over liquid species
    for_fields!(du.species_liquid_mass_per_volume, du.volumetric_liquid_mass_generation, u.liquid_mass_fractions) do species, liquid_accumulation, liquid_generation, liquid_mass_fractions
        upwinded_fraction = liquid_velocity_face >= 0.0 ? liquid_mass_fractions[species[idx_a]] : liquid_mass_fractions[species[idx_b]]
        change_in_volume = - (liquid_flux_density * upwinded_fraction * area / vol_a)
        if face_idx == 1
            change_in_volume += liquid_generation[species[idx_a]]
        end
        liquid_accumulation[species[idx_a]] += change_in_volume
        if hasproperty(du, :overall_liquid_mass_per_volume)
            du.overall_liquid_mass_per_volume[idx_a] += change_in_volume
        end
    end
end

function drift_flux_mixture_momentum!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    # Average bed porosity and other properties
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])
    
    mixture_viscosity_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_viscosity[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_viscosity[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_viscosity[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_viscosity[idx_b])
    )
    
    mixture_density_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_density[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_density[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_density[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_density[idx_b])
    )
    
    solid_particle_diameter_average = 0.5 * (u.solid_particle_diameter[idx_a] + u.solid_particle_diameter[idx_b])
    vessel_hydraulic_diameter_average = 0.5 * (u.vessel_hydraulic_diameter[idx_a] + u.vessel_hydraulic_diameter[idx_b])
    
    # Calculate viscous and inertial drag coefficients on the face
    ergun_viscous = 150.0 * (1.0 - bed_porosity_average)^2 / (solid_particle_diameter_average^2 * bed_porosity_average^3) * mixture_viscosity_average
    ergun_inertial = 1.75 * (1.0 - bed_porosity_average) / (solid_particle_diameter_average * bed_porosity_average^3) * mixture_density_average
    
    viscous_drag = (bed_porosity_average < 1.0 ? ergun_viscous : 0.0) + 32.0 * mixture_viscosity_average / (vessel_hydraulic_diameter_average^2)
    inertial_drag = (bed_porosity_average < 1.0 ? ergun_inertial : 0.0) + 2.0 * 0.02 * mixture_density_average / vessel_hydraulic_diameter_average
    
    driving_force = (u.pressure[idx_a] - u.pressure[idx_b]) / dist + mixture_density_average * LinearAlgebra.dot(Ferrite.Vec{3}((0.0, 0.0, -9.81)), norm)
    
    if inertial_drag > 1e-12
        superficial_velocity_magnitude = (-viscous_drag + sqrt(viscous_drag^2 + 4.0 * inertial_drag * abs(driving_force))) / (2.0 * inertial_drag)
    else
        superficial_velocity_magnitude = abs(driving_force) / viscous_drag
    end
    mixture_superficial_velocity = superficial_velocity_magnitude * sign(driving_force)
    
    if hasproperty(du, :mass_face)
        du.mass_face[idx_a, face_idx] -= mixture_density_average * mixture_superficial_velocity * area
    end
end

function drift_flux_energy_balance!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    # Average bed porosity and other properties
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])
    
    mixture_density_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_density[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_density[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_density[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_density[idx_b])
    )
    
    mixture_viscosity_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_viscosity[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_viscosity[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_viscosity[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_viscosity[idx_b])
    )
    
    solid_particle_diameter_average = 0.5 * (u.solid_particle_diameter[idx_a] + u.solid_particle_diameter[idx_b])
    vessel_hydraulic_diameter_average = 0.5 * (u.vessel_hydraulic_diameter[idx_a] + u.vessel_hydraulic_diameter[idx_b])
    
    # Calculate viscous and inertial drag coefficients on the face
    ergun_viscous = 150.0 * (1.0 - bed_porosity_average)^2 / (solid_particle_diameter_average^2 * bed_porosity_average^3) * mixture_viscosity_average
    ergun_inertial = 1.75 * (1.0 - bed_porosity_average) / (solid_particle_diameter_average * bed_porosity_average^3) * mixture_density_average
    
    viscous_drag = (bed_porosity_average < 1.0 ? ergun_viscous : 0.0) + 32.0 * mixture_viscosity_average / (vessel_hydraulic_diameter_average^2)
    inertial_drag = (bed_porosity_average < 1.0 ? ergun_inertial : 0.0) + 2.0 * 0.02 * mixture_density_average / vessel_hydraulic_diameter_average
    
    driving_force = (u.pressure[idx_a] - u.pressure[idx_b]) / dist + mixture_density_average * LinearAlgebra.dot(Ferrite.Vec{3}((0.0, 0.0, -9.81)), norm)
    
    if inertial_drag > 1e-12
        superficial_velocity_magnitude = (-viscous_drag + sqrt(viscous_drag^2 + 4.0 * inertial_drag * abs(driving_force))) / (2.0 * inertial_drag)
    else
        superficial_velocity_magnitude = abs(driving_force) / viscous_drag
    end
    mixture_superficial_velocity = superficial_velocity_magnitude * sign(driving_force)
    mixture_pore_velocity = mixture_superficial_velocity / bed_porosity_average
    
    gas_velocity_face = 0.5 * (u.distribution_parameter[idx_a] + u.distribution_parameter[idx_b]) * mixture_pore_velocity + 0.5 * (u.local_gas_drift_velocity[idx_a] + u.local_gas_drift_velocity[idx_b])
    gas_holdup_clamped = clamp(0.5 * (u.gas_holdup[idx_a] + u.gas_holdup[idx_b]), 0.0, 0.9999)
    liquid_velocity_face = (mixture_pore_velocity - gas_holdup_clamped * gas_velocity_face) / (1.0 - gas_holdup_clamped)
    
    # Upwind phase states and enthalpies
    if gas_velocity_face >= 0.0
        upwinded_gas_density = u.gas_density[idx_a]
        upwinded_gas_holdup = u.gas_holdup[idx_a]
        upwinded_gas_enthalpy = u.gas_enthalpy[idx_a]
    else
        upwinded_gas_density = u.gas_density[idx_b]
        upwinded_gas_holdup = u.gas_holdup[idx_b]
        upwinded_gas_enthalpy = u.gas_enthalpy[idx_b]
    end
    
    if liquid_velocity_face >= 0.0
        upwinded_liquid_density = u.liquid_density[idx_a]
        upwinded_liquid_holdup = 1.0 - u.gas_holdup[idx_a]
        upwinded_liquid_enthalpy = u.liquid_enthalpy[idx_a]
    else
        upwinded_liquid_density = u.liquid_density[idx_b]
        upwinded_liquid_holdup = 1.0 - u.gas_holdup[idx_b]
        upwinded_liquid_enthalpy = u.liquid_enthalpy[idx_b]
    end
    
    # Calculate energy fluxes (in Watts/m^2)
    gas_energy_flux_density = bed_porosity_average * upwinded_gas_density * upwinded_gas_holdup * gas_velocity_face * upwinded_gas_enthalpy
    liquid_energy_flux_density = bed_porosity_average * upwinded_liquid_density * upwinded_liquid_holdup * liquid_velocity_face * upwinded_liquid_enthalpy
    total_energy_flux_density = gas_energy_flux_density + liquid_energy_flux_density
    
    # Energy transport volumetric change rate (W/m^3)
    change_in_enthalpy_per_volume = - (total_energy_flux_density * area / vol_a)
    
    # Update cell energy
    du.heat[idx_a] += change_in_enthalpy_per_volume * vol_a
end
