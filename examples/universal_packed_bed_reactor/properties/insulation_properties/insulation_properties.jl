function get_insulation_properties()
    insulation_properties = ComponentVector(
        k = 0.10u"W/(m*K)",
        cp = 950.0u"J/(kg*K)",
        rho = 80.0u"kg/m^3"
    )
    return insulation_properties
end

