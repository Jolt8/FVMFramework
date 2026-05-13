function get_water_properties()
    water_properties = ComponentVector(
        fluid_k = 0.6u"W/(m*K)", 
        fluid_cp = 4186u"J/(kg*K)",
        dynamic_viscosity = 1.0e-3u"Pa*s",
        fluid_rho = 998.2u"kg/m^3",

        #saturation_temp = 100.0u"°C",
        #liquid_rho = 997.0u"kg/m^3",
        #gas_rho = 0.4u"kg/m^3",
        #mass_transfer_coeff_vap = 0.001u"kg/(m^2*s*K)",
        #heat_of_vaporization = 2260.0u"kJ/kg",

        diffusion_coefficients = (
            methanol = 1e-5u"m^2/s",
            water = 1e-5u"m^2/s",
            carbon_monoxide = 1e-5u"m^2/s",
            hydrogen = 1e-5u"m^2/s",
            carbon_dioxide = 1e-5u"m^2/s",
            air = 1e-5u"m^2/s"
        ), 
    )

    return water_properties
end