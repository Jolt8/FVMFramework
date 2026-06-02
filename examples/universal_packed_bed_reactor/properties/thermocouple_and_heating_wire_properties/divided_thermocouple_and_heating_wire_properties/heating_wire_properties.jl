Revise.includet(joinpath(@__DIR__, "..", "..", "insulation_properties", "insulation_properties.jl"))
insulation_properties = get_insulation_properties()

insulation_k = insulation_properties.k
insulation_cp = insulation_properties.cp
insulation_rho = insulation_properties.rho

function get_heating_wire_properties(grid, pipe_length, common_properties)
    kanthal_heating_wire_k = 18.0u"W/(m*K)"
    kanthal_heating_wire_cp = 500.0u"J/(kg*K)"
    kanthal_heating_wire_rho = 7800.0u"kg/m^3"

    heating_wire_spacing_between_turns = 10.0u"inch" / 30.6 
    heating_wire_diameter = common_properties.heating_wire_outside_diameter - common_properties.heating_wire_inside_diameter
    heating_wire_length = 2.7432u"m" #length of heating wire actually used
    heating_section_length = 10.0u"inch" #length along the pipe the heating wire is coiled around on

    heating_wire_volume = pi * (heating_wire_diameter/2)^2 * heating_wire_length 

    total_heating_section_volume = pi * (common_properties.heating_wire_outside_diameter/2)^2 * common_properties.pipe_length - 
        pi * (common_properties.heating_wire_inside_diameter/2)^2 * common_properties.pipe_length
    #the heating wire is only coiled for 10 inches along the pipe, so the 1 inch on the bottom and top of the pipe are filled in with the 
    #total heating section volume

    heating_wire_section_non_heater_void = (total_heating_section_volume - heating_wire_volume) / total_heating_section_volume

    heating_wire_properties = ComponentVector(
        k = heating_wire_section_non_heater_void * insulation_k + (1 - heating_wire_section_non_heater_void) * kanthal_heating_wire_k,
        cp = heating_wire_section_non_heater_void * insulation_cp + (1 - heating_wire_section_non_heater_void) * kanthal_heating_wire_cp,
        rho = heating_wire_section_non_heater_void * insulation_rho + (1 - heating_wire_section_non_heater_void) * kanthal_heating_wire_rho,
    )

    return heating_wire_properties
end