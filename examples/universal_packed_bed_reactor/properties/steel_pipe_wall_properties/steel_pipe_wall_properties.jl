function get_steel_pipe_wall_properties()
    steel_pipe_wall_properties = ComponentVector(
        k = 16.2u"W/(m*K)",
        cp = 500.0u"J/(kg*K)",
        rho = 7900.0u"kg/m^3",
        
    )
    return steel_pipe_wall_properties
end

