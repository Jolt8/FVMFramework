function get_water_properties()
    water_properties = ComponentVector(
        R_gas = 8.314u"J/(mol*K)",

        molecular_weights = (
            methanol = 0.03204u"kg/mol",
            water = 0.01802u"kg/mol",
            #carbon_monoxide = 28.01u"g/mol",
            #hydrogen = 2.02u"g/mol",
            #carbon_dioxide = 44.01u"g/mol",
            #air = 28.97u"g/mol"
        ),
        
        k = 0.6u"W/(m*K)",
        cp = 3000u"J/(kg*K)",
        fluid_rho = 1000u"kg/m^3",
        compressibility = 4.47e-10u"Pa^-1"
    )

    return water_properties
end