function get_silicon_carbide_sand_properties()
    silicon_carbide_sand_properties = ComponentVector(
        silicon_carbide_k = 150.0u"W/(m*K)", #this could be changed to packing_k later, but this is fine for now
        silicon_carbide_cp = 700.0u"J/(kg*K)",
        silicon_carbide_rho = 3150.0u"kg/m^3",
        
        bed_void_fraction = 0.504,
        particle_diameter = 0.200u"mm"
    )

    return silicon_carbide_sand_properties
end
