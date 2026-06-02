using LinearAlgebra


#NOTE: there are three different types of densities that are tracked by the solver
    #- gas/liquid_densities is the solved for variable that tracks that density of each species in its gaseous/liquid phase if it were to exist on its own
    #- overall_gas/liquid_density is the overall density of the gas/liquid phase in the cell (summed from gas/liquid_densities)
    #- eos_gas/liquid_density is the density of the gas/liquid phase in the cell as determined by the equation of state for the gaseous/liquid phase
    #- fluid_rho is the overall density of the fluid in the cell (determined with overall_gas_density and overall_liquid_density)
    #- solid_rho is the density of the solid packing (if any)
    #- rho is the combination of fluid_rho and solid_rho based on bed_proosity

function liquid_and_gas_density!(du, u, cell_id, vol)
    u.overall_gas_density[cell_id] *= 0.0
    u.overall_liquid_density[cell_id] *= 0.0
    for_fields!(u.gas_densities, u.liquid_densities) do species, u_gas_densities, u_liquid_densities
        u.overall_gas_density[cell_id] += u_gas_densities[species[cell_id]]
        u.overall_liquid_density[cell_id] += u_liquid_densities[species[cell_id]]
    end
end

function liquid_and_gas_mass_fractions!(du, u, cell_id, vol)
    for_fields!(u.gas_densities, u.liquid_densities) do species, u_gas_densities, u_liquid_densities
        u.gas_mass_fractions[species[cell_id]] = u_gas_densities[species[cell_id]] / u.overall_gas_density[cell_id]
        u.liquid_mass_fractions[species[cell_id]] = u_liquid_densities[species[cell_id]] / u.overall_liquid_density[cell_id]
    end
end


function gas_and_liquid_mw_avg!(du, u, cell_id, vol)
    sum_gas_moles = 0.0
    sum_liquid_moles = 0.0
    for_fields!(u.gas_densities, u.liquid_densities, u.molecular_weights) do species, u_gas_densities, u_liquid_densities, u_molecular_weights
        sum_gas_moles += u_gas_densities[species[cell_id]] / u_molecular_weights[species[cell_id]]
        sum_liquid_moles += u_liquid_densities[species[cell_id]] / u_molecular_weights[species[cell_id]]
    end
    if sum_gas_moles > 0.0
        u.gas_mw_avg[cell_id] = u.overall_gas_density[cell_id] / sum_gas_moles
    else
        u.gas_mw_avg[cell_id] = 0.01802
    end
    if sum_liquid_moles > 0.0
        u.liquid_mw_avg[cell_id] = u.overall_liquid_density[cell_id] / sum_liquid_moles
    else
        u.liquid_mw_avg[cell_id] = 0.01802
    end
end

function two_phase_mw_avg!(du, u, cell_id, vol)
    total_density = u.overall_gas_density[cell_id] + u.overall_liquid_density[cell_id]
    total_molar_density = u.overall_gas_density[cell_id] / u.gas_mw_avg[cell_id] + u.overall_liquid_density[cell_id] / u.liquid_mw_avg[cell_id]
    if total_molar_density > 0.0
        u.mw_avg[cell_id] = total_density / total_molar_density
    else
        u.mw_avg[cell_id] = 0.01802
    end
end


#update_eos_densities! would go here

function liquid_and_gas_holdup!(du, u, cell_id, vol)
    gas_volume = u.overall_gas_density[cell_id] / u.eos_gas_density[cell_id]
    liquid_volume = u.overall_liquid_density[cell_id] / u.eos_liquid_density[cell_id]

    sum_volume = gas_volume + liquid_volume
    if sum_volume > 0.0
        u.gas_holdup[cell_id] = gas_volume / sum_volume
    else
        u.gas_holdup[cell_id] = 0.0
    end
    
    u.liquid_holdup[cell_id] = 1.0 - u.gas_holdup[cell_id]
end

function fluid_rho!(du, u, cell_id, vol)
    u.fluid_rho[cell_id] = u.overall_gas_density[cell_id] + u.overall_liquid_density[cell_id]
end

function gas_and_liquid_mole_fractions!(du, u, cell_id, vol)
    for_fields!(u.gas_mole_fractions, u.liquid_mole_fractions, u.gas_mass_fractions, u.liquid_mass_fractions, u.molecular_weights) do species, u_gas_mole_fractions, u_liquid_mole_fractions, u_gas_mass_fractions, u_liquid_mass_fractions, u_molecular_weights
        u_gas_mole_fractions[species[cell_id]] = u_gas_mass_fractions[species[cell_id]] * u.gas_mw_avg[cell_id] / u_molecular_weights[species[cell_id]]
        u_liquid_mole_fractions[species[cell_id]] = u_liquid_mass_fractions[species[cell_id]] * u.liquid_mw_avg[cell_id] / u_molecular_weights[species[cell_id]]
    end
end

#update_K_vle! would go here

function phase_change_kinetic_mass_transfer!(du, u, cell_id, vol)
    for_fields!(du.gas_densities, du.liquid_densities, u.liquid_densities, u.gas_densities, u.K_vle, u.liquid_mole_fractions, u.gas_mole_fractions, u.specific_heat_of_vaporizations, u.molecular_weights) do species, du_gas_densities, du_liquid_densities, u_liquid_densities, u_gas_densities, u_K_vle, u_liquid_mole_fractions, u_gas_mole_fractions, u_specific_heat_of_vaporization, u_molecular_weights
        driving_force = u_K_vle[species[cell_id]] * u_liquid_mole_fractions[species[cell_id]] - u_gas_mole_fractions[species[cell_id]] 
        effective_bubble_area = u.bubble_area_per_m3[cell_id] * u.gas_holdup[cell_id] * 0.00001 #maybe the 0.00001 could be an experimentally derived value
        kinetic_constant = u.phase_change_mass_transfer_coefficient[cell_id] * effective_bubble_area
        #hmm, whenever I multiply the kinetic_constant by 0.00001, it can get to 3000 seconds before the stiffness sets in
        
        #stiff 8 seconds in
        
        # Use a single rate scaled by available source-phase density
        
        if driving_force >= 0.0
            # Evaporation: limited by available liquid
            evaporation_rate = kinetic_constant * u_liquid_densities[species[cell_id]] * driving_force
        else
            # Condensation: limited by available gas  
            evaporation_rate = kinetic_constant * u_gas_densities[species[cell_id]] * driving_force
        end
        
        

        
        #stiff 0.05 seconds in
        #=
        ε = 1e-6
        evap_weight = 0.5 * (1.0 + tanh(driving_force / ε))
        reference_density = evap_weight * u_liquid_densities[species[cell_id]] + (1.0 - evap_weight) * u_gas_densities[species[cell_id]]
        evaporation_rate = kinetic_constant * reference_density * driving_force
        =#
        

        #stiff 1.66 seconds in
        #=
        reference_density = sqrt(
            max(u_liquid_densities[species[cell_id]], 1e-20) *
            max(u_gas_densities[species[cell_id]], 1e-20)
        )
        evaporation_rate = kinetic_constant * reference_density * driving_force
        =#
        
        #stiff 5 seconds in
        #=
        liq = max(u_liquid_densities[species[cell_id]], 1e-20)
        gas = max(u_gas_densities[species[cell_id]], 1e-20)
        reference_density = 2.0 * liq * gas / (liq + gas)
        evaporation_rate = kinetic_constant * reference_density * driving_force
        =#
        

        #bubble_area_per_m3 represents the surface area of the bubbles in contact with the liquid phase
        #stiff 7 seconds in
        #evaporation_rate = u.phase_change_mass_transfer_coefficient[cell_id] * u.bubble_area_per_m3[cell_id] * u.gas_densities[species[cell_id]] * (u_K_vle[species[cell_id]] * u_liquid_mole_fractions[species[cell_id]] - u_gas_mole_fractions[species[cell_id]])

        #remember, evaporation_rate has units kg/(m^3*s)
        du.heat[cell_id] -= evaporation_rate * vol * u_specific_heat_of_vaporization[species[cell_id]]

        #stiff 0.2 seconds in
        #= 
        c_ref = u.eos_liquid_density[cell_id] / u.liquid_mw_avg[cell_id]
        kinetic_constant = u.phase_change_mass_transfer_coefficient[cell_id] * u.bubble_area_per_m3[cell_id]
        # Molar rate [mol/(m³·s)]
        molar_rate = kinetic_constant * c_ref * driving_force
        # Mass rate for species i [kg/(m³·s)]
        evaporation_rate = molar_rate * u_molecular_weights[species[cell_id]]
        =#

        #just for reference, if phase change is turned off, the system solves the entire 3000 second problem in 3.2 seconds

        #@show u_K_vle[species[cell_id]] u_K_vle comes back as 2 and 10 usually
        #@show driving_force driving force comes back as 0.79 or 3.0
        
        du_gas_densities[species[cell_id]] += evaporation_rate
        du_liquid_densities[species[cell_id]] -= evaporation_rate
    end
end

#NOTE: stuff that must be updated externally:
    # - distribution_parameter
    # - local_gas_drift_velocity
    # - void_fraction/bed_porosity
    # - solid_particle_diameter
    # - vessel_hydraulic_diameter

#NOTE: stuff required for cached variables
    #- viscous_drag (per face)
    #- inertial_drag (per face)
    #- driving_force (per face)
    #- mixture_superficial_velocity (per face)
    #- mixture_pore_velocity (per face)
    #- driving_force (per face)
    #- mass_face (per face)
    
    #- mass (per cell)
    #- gas_holdup (per cell)
    #- liquid_holdup (per cell)
    #- heat (per cell)

#variables that are solved for
    #- gas_densities (per cell)
    #- liquid_densities (per cell)
    #- pressure (per cell)
    #- temp (per cell)



#Cumulative list of everything empirical
    #- distribution_parameter
    #- local_gas_drift_velocity
    #- update_eos_densities!
    #- update_K_vle!

#fluxes
function update_face_ergun_drag!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])

    mixture_viscosity_average = 0.5 * (
        (u.gas_holdup[idx_a] * u.gas_viscosity[idx_a] + (1.0 - u.gas_holdup[idx_a]) * u.liquid_viscosity[idx_a]) + 
        (u.gas_holdup[idx_b] * u.gas_viscosity[idx_b] + (1.0 - u.gas_holdup[idx_b]) * u.liquid_viscosity[idx_b])
    )
    
    mixture_density_average = 0.5 * (u.fluid_rho[idx_a] + u.fluid_rho[idx_b])
    
    solid_particle_diameter_average = 0.5 * (u.solid_particle_diameter[idx_a] + u.solid_particle_diameter[idx_b])
    vessel_hydraulic_diameter_average = 0.5 * (u.vessel_hydraulic_diameter[idx_a] + u.vessel_hydraulic_diameter[idx_b])
    
    # Calculate viscous and inertial drag coefficients on the face
    ergun_viscous = 150.0 * (1.0 - bed_porosity_average)^2 / (solid_particle_diameter_average^2 * bed_porosity_average^3) * mixture_viscosity_average
    ergun_inertial = 1.75 * (1.0 - bed_porosity_average) / (solid_particle_diameter_average * bed_porosity_average^3) * mixture_density_average
    
    if bed_porosity_average < 1.0
        u.viscous_drag[idx_a, face_idx] = ergun_viscous + 32.0 * mixture_viscosity_average / (vessel_hydraulic_diameter_average^2)
        u.inertial_drag[idx_a, face_idx] = ergun_inertial + 2.0 * 0.02 * mixture_density_average / vessel_hydraulic_diameter_average
    else
        u.viscous_drag[idx_a, face_idx] = 32.0 * mixture_viscosity_average / (vessel_hydraulic_diameter_average^2)
        u.inertial_drag[idx_a, face_idx] = 2.0 * 0.02 * mixture_density_average / vessel_hydraulic_diameter_average
    end
end

function update_face_driving_force_and_velocities!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])
    
    mixture_density_average = 0.5 * (u.fluid_rho[idx_a] + u.fluid_rho[idx_b])

    u.driving_force[idx_a, face_idx] = (u.pressure[idx_a] - u.pressure[idx_b]) / dist + mixture_density_average * LinearAlgebra.dot(Ferrite.Vec{3}((0.0, 0.0, -9.81)), norm)
    
    if u.inertial_drag[idx_a, face_idx] > 1e-12
        superficial_velocity_magnitude = (-u.viscous_drag[idx_a, face_idx] + sqrt(u.viscous_drag[idx_a, face_idx]^2 + 4.0 * u.inertial_drag[idx_a, face_idx] * abs(u.driving_force[idx_a, face_idx]))) / (2.0 * u.inertial_drag[idx_a, face_idx])
    else
        superficial_velocity_magnitude = abs(u.driving_force[idx_a, face_idx]) / u.viscous_drag[idx_a, face_idx]
    end
    
    u.mixture_superficial_velocity[idx_a, face_idx] = superficial_velocity_magnitude * sign(u.driving_force[idx_a, face_idx])
    u.mixture_pore_velocity[idx_a, face_idx] = u.mixture_superficial_velocity[idx_a, face_idx] / bed_porosity_average

    gas_holdup_clamped = clamp(0.5 * (u.gas_holdup[idx_a] + u.gas_holdup[idx_b]), 0.0, 0.9999) #we need to come back to this to see if gas_holdup actually does go above 1.0
    
    average_distribution_parameter = 0.5 * (u.distribution_parameter[idx_a] + u.distribution_parameter[idx_b])
    drift_velocity = 0.5 * (u.local_gas_drift_velocity[idx_a] + u.local_gas_drift_velocity[idx_b]) * norm[3]

    u.gas_velocity_face[idx_a, face_idx] = average_distribution_parameter * u.mixture_pore_velocity[idx_a, face_idx] + drift_velocity * (1.0 - gas_holdup_clamped)
    
    # Algebraically simplified v_l to avoid division by near-zero when gas_holdup approaches 1:
    u.liquid_velocity_face[idx_a, face_idx] = ((1.0 - gas_holdup_clamped * average_distribution_parameter) / (1.0 - gas_holdup_clamped)) * u.mixture_pore_velocity[idx_a, face_idx] - gas_holdup_clamped * drift_velocity
end

function drift_flux_mass_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    #remember to call update_face_ergun_drag!() and update_face_driving_force_and_velocities!()
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])

    # Upwind phase states
    if u.gas_velocity_face[idx_a, face_idx] >= 0.0 
        upwinded_gas_density = u.overall_gas_density[idx_a]
        upwinded_gas_holdup = u.gas_holdup[idx_a]
    else
        upwinded_gas_density = u.overall_gas_density[idx_b]
        upwinded_gas_holdup = u.gas_holdup[idx_b]
    end
    
    if u.liquid_velocity_face[idx_a, face_idx] >= 0.0
        upwinded_liquid_density = u.overall_liquid_density[idx_a]
        upwinded_liquid_holdup = 1.0 - u.gas_holdup[idx_a]
    else
        upwinded_liquid_density = u.overall_liquid_density[idx_b]
        upwinded_liquid_holdup = 1.0 - u.gas_holdup[idx_b]
    end
    
    gas_flux_density = bed_porosity_average * upwinded_gas_density * u.gas_velocity_face[idx_a, face_idx]
    liquid_flux_density = bed_porosity_average * upwinded_liquid_density * u.liquid_velocity_face[idx_a, face_idx]
    
    # Loop over gas species
    for_fields!(du.gas_densities, u.gas_mass_fractions) do species, du_gas_densities, u_gas_mass_fractions
        if u.gas_velocity_face[idx_a, face_idx] >= 0.0
            upwinded_fraction = u_gas_mass_fractions[species[idx_a]]
        else
            upwinded_fraction = u_gas_mass_fractions[species[idx_b]]
        end

        gas_accumulation_rate = - (gas_flux_density * upwinded_fraction * area / (vol_a * u.bed_porosity[idx_a]))

        du_gas_densities[species[idx_a]] += gas_accumulation_rate
    end
    
    # Loop over liquid species
    for_fields!(du.liquid_densities, u.liquid_mass_fractions) do species, du_liquid_densities, u_liquid_mass_fractions
        if u.liquid_velocity_face[idx_a, face_idx] >= 0.0
            upwinded_fraction = u_liquid_mass_fractions[species[idx_a]]
        else
            upwinded_fraction = u_liquid_mass_fractions[species[idx_b]]
        end
        
        liquid_accumulation_rate = - (liquid_flux_density * upwinded_fraction * area / (vol_a * u.bed_porosity[idx_a]))

        du_liquid_densities[species[idx_a]] += liquid_accumulation_rate
    end
end

function drift_flux_mixture_momentum!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    #remember to call update_face_ergun_drag!() and update_face_driving_force_and_velocities!()

    mixture_density_average = 0.5 * (u.fluid_rho[idx_a] + u.fluid_rho[idx_b])

    du.mass_face[idx_a, face_idx] -= mixture_density_average * u.mixture_superficial_velocity[idx_a, face_idx] * area
end

function drift_flux_energy_balance!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    bed_porosity_average = 0.5 * (u.bed_porosity[idx_a] + u.bed_porosity[idx_b])

    # Upwind phase states and enthalpies
    if u.gas_velocity_face[idx_a, face_idx] >= 0.0
        upwinded_gas_density = u.overall_gas_density[idx_a]
        upwinded_gas_enthalpy = u.gas_enthalpy[idx_a]
    else
        upwinded_gas_density = u.overall_gas_density[idx_b]
        upwinded_gas_enthalpy = u.gas_enthalpy[idx_b]
    end
    
    if u.liquid_velocity_face[idx_a, face_idx] >= 0.0
        upwinded_liquid_density = u.overall_liquid_density[idx_a]
        upwinded_liquid_enthalpy = u.liquid_enthalpy[idx_a]
    else
        upwinded_liquid_density = u.overall_liquid_density[idx_b]
        upwinded_liquid_enthalpy = u.liquid_enthalpy[idx_b]
    end
    
    # Calculate energy fluxes (in Watts/m^2)
    gas_energy_flux_density = bed_porosity_average * upwinded_gas_density  * u.gas_velocity_face[idx_a, face_idx] * upwinded_gas_enthalpy
    liquid_energy_flux_density = bed_porosity_average * upwinded_liquid_density * u.liquid_velocity_face[idx_a, face_idx] * upwinded_liquid_enthalpy
    total_energy_flux_density = gas_energy_flux_density + liquid_energy_flux_density
    
    # Energy transport volumetric change rate (W/m^3)
    change_in_enthalpy_per_volume = - (total_energy_flux_density * area / vol_a)
    
    # Update cell energy
    du.heat[idx_a] += change_in_enthalpy_per_volume * vol_a #note that heat is not W/m^3, but just W, it's capped by the cell's volume later
end

Revise.includet(joinpath(@__DIR__, "eos_stuff.jl")) #for update_eos_densities! and update_K_vle!

function overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model)
    liquid_and_gas_density!(du, u, cell_id, vol)
    liquid_and_gas_mass_fractions!(du, u, cell_id, vol)
    gas_and_liquid_mw_avg!(du, u, cell_id, vol)
    two_phase_mw_avg!(du, u, cell_id, vol)
    gas_and_liquid_mole_fractions!(du, u, cell_id, vol)
    update_eos_densities!(du, u, cell_id, vol, clapeyron_model) #from eos_stuff.jl
    liquid_and_gas_holdup!(du, u, cell_id, vol)
    fluid_rho!(du, u, cell_id, vol)
    update_K_vle!(du, u, cell_id, vol, clapeyron_model) #from eos_stuff.jl
    phase_change_kinetic_mass_transfer!(du, u, cell_id, vol)
end

function overall_drift_flux_model!(du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    update_face_ergun_drag!(
        du, u,
        idx_a, idx_b, face_idx,
        area, norm, dist,
        vol_a, vol_b
    )
    update_face_driving_force_and_velocities!(
        du, u,
        idx_a, idx_b, face_idx,
        area, norm, dist,
        vol_a, vol_b
    )
    drift_flux_mass_flux!(
        du, u,
        idx_a, idx_b, face_idx,
        area, norm, dist,
        vol_a, vol_b
    )
    drift_flux_mixture_momentum!(
        du, u,
        idx_a, idx_b, face_idx,
        area, norm, dist,
        vol_a, vol_b
    )
    drift_flux_energy_balance!(
        du, u,
        idx_a, idx_b, face_idx,
        area, norm, dist,
        vol_a, vol_b
    )
end
