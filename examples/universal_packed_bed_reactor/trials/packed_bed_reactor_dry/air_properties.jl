function get_air_properties()
    air_properties = ComponentVector(
        fluid_k = 0.0256u"W/(m*K)", 
        fluid_cp = 1007u"J/(kg*K)",
        dynamic_viscosity = 1.80e-5u"Pa*s",
        fluid_rho = 1.225u"kg/m^3",

        diffusion_coefficients = (
            methanol = 1e-5u"m^2/s",
            water = 1e-5u"m^2/s",
            carbon_monoxide = 1e-5u"m^2/s",
            hydrogen = 1e-5u"m^2/s",
            carbon_dioxide = 1e-5u"m^2/s",
            air = 1e-5u"m^2/s"
        ), 
    )

    return air_properties
end