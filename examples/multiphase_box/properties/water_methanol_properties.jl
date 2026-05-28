function get_water_methanol_properties()
    water_methanol_properties = ComponentVector(
        R_gas = 8.314u"J/(mol*K)",

        molecular_weights = (
            methanol = 0.03204u"kg/mol",
            water = 0.01802u"kg/mol",
            #carbon_monoxide = 28.01u"g/mol",
            #hydrogen = 2.02u"g/mol",
            #carbon_dioxide = 44.01u"g/mol",
            #air = 28.97u"g/mol"
        ),
        distribution_parameter = 1.0,
        local_gas_drift_velocity = 0.25u"m/s",

        #if we had experimental data, these two coefficients would probably be lumped together and a neural network would be used to figure this out
        phase_change_mass_transfer_coefficient = 1e-3u"m/s", 
        bubble_area_per_m3 = 300.0u"m^2/m^3",

        k = 0.6u"W/(m*K)",
        cp = 3000u"J/(kg*K)",
        gas_viscosity = 1e-5u"Pa*s",
        liquid_viscosity = 1e-3u"Pa*s",

        solid_rho = 8960u"kg/m^3",
        bed_porosity = 0.4,
        solid_particle_diameter = 10.0u"mm",
        vessel_hydraulic_diameter = 1.0u"m",

        specific_heat_of_vaporizations = (
            methanol = 2260u"J/g",
            water = 1100u"J/g",
        )
    )

    return water_methanol_properties
end