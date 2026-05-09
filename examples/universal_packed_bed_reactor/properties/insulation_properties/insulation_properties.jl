function get_insulation_properties()
    insulation_properties = ComponentVector(
        k = 0.15u"W/(m*K)",
        cp = 1000.0u"J/(kg*K)",
        rho = 100.0u"kg/m^3"
    )
    return insulation_properties
end

