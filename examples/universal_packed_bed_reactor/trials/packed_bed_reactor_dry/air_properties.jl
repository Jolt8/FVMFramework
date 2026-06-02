function get_air_properties()
    air_properties = ComponentVector(
        fluid_k = 0.0256u"W/(m*K)", 
        fluid_cp = 1007u"J/(kg*K)",
        dynamic_viscosity = 1.80e-5u"Pa*s",
        fluid_rho = 1.225u"kg/m^3",
    )

    return air_properties
end