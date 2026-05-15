Revise.includet(joinpath(@__DIR__, "..", "..", "insulation_properties", "insulation_properties.jl"))
insulation_properties = get_insulation_properties()

insulation_k = insulation_properties.k
insulation_cp = insulation_properties.cp
insulation_rho = insulation_properties.rho

function get_heating_wire_properties(grid, pipe_length, n_cells_axial, cell_lengths_along_pipe)

    air_k = 0.0257u"W/(m*K)"
    air_cp = 1.007u"J/(g*K)"
    air_rho = 1.225u"kg/m^3"

    #we may acutally want to use the properties of the insulation because it was wrapped on top of the heating wire, 
    #filling the gaps between the wires

    kanthal_heating_wire_k = 18.0u"W/(m*K)"
    kanthal_heating_wire_cp = 500.0u"J/(kg*K)"
    kanthal_heating_wire_rho = 7800.0u"kg/m^3"

    heating_wire_spacing_between_turns = 10.0u"inch" / 30.6 
    heating_wire_thickness = 0.7u"mm" #thickness of the heating wire
    heating_wire_length = 2.7432u"m"

    heating_wire_volume = pi * heating_wire_thickness^2 / 4 * heating_wire_length
    heating_section_volume = pi * (common_properties.thermocouple_region_outside_diameter/2)^2 * pipe_length - 
        pi * (common_properties.insulation_inside_diameter/2)^2 * pipe_length

    heating_wire_section_non_heater_void = (heating_section_volume - heating_wire_volume) / heating_section_volume

    heating_wire_properties = ComponentVector(
        k = heating_wire_section_non_heater_void * air_k + (1 - heating_wire_section_non_heater_void) * kanthal_heating_wire_k,
        cp = heating_wire_section_non_heater_void * air_cp + (1 - heating_wire_section_non_heater_void) * kanthal_heating_wire_cp,
        rho = heating_wire_section_non_heater_void * air_rho + (1 - heating_wire_section_non_heater_void) * kanthal_heating_wire_rho,
    )

    return heating_wire_properties
end