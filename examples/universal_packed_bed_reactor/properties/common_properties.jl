#huh, I guess these are the only properties that are really shared between the different property groups
#I guess I just won't use it

function get_common_properties(pipe_length, cell_lengths_along_pipe, n_cells_axial)
    #Reforming Area
    reforming_area_outside_diameter = 16u"mm"
    reforming_area_cell_volumes = pi * (reforming_area_outside_diameter / 2)^2 * (pipe_length / n_cells_axial)

    #Reforming Area to Steel Pipe Wall 
    reforming_area_to_steel_pipe_wall_area = pi * reforming_area_outside_diameter * pipe_length
    reforming_area_to_steel_pipe_wall_cell_areas = reforming_area_to_steel_pipe_wall_area / n_cells_axial

    #Steel Pipe Wall 
    steel_pipe_wall_inside_diameter = 16u"mm"
    steel_pipe_wall_outside_diameter = 18u"mm"
    steel_pipe_wall_cell_volumes = pi * (steel_pipe_wall_outside_diameter / 2)^2 * (pipe_length / n_cells_axial) - pi * (steel_pipe_wall_inside_diameter / 2)^2 * (pipe_length / n_cells_axial)

    #Steel Pipe Wall to Thermocouple Region
    steel_pipe_wall_to_thermocouple_area = pi * steel_pipe_wall_outside_diameter * pipe_length
    steel_pipe_wall_to_thermocouple_cell_areas = steel_pipe_wall_to_thermocouple_area / n_cells_axial

    #Thermocouple Region
    thermocouple_region_inside_diameter = 18u"mm"
    thermocouple_region_outside_diameter = 30u"mm" #just an estimate, not very certain on this one 
    thermocouple_region_cell_volumes = pi * (thermocouple_region_outside_diameter / 2)^2 * (pipe_length / n_cells_axial) - pi * (thermocouple_region_inside_diameter / 2)^2 * (pipe_length / n_cells_axial)

    #Thermocouple Region to Insulation
    thermocouple_region_to_insulation_area = pi * thermocouple_region_outside_diameter * pipe_length
    thermocouple_region_to_insulation_cell_areas = thermocouple_region_to_insulation_area / n_cells_axial

    #Insulation
    insulation_inside_diameter = 30u"mm" #there's heating wire that may cause a slight gap, but it's negligible
    insulation_outside_diameter = insulation_inside_diameter + 1.0u"inch" #this will be accurate assuming insulation_inside_diameter is accurate
    insulation_cell_volumes = pi * (insulation_outside_diameter / 2)^2 * (pipe_length / n_cells_axial) - pi * (insulation_inside_diameter / 2)^2 * (pipe_length / n_cells_axial)

    #Insulation to Environment
    insulation_to_environment_area = pi * insulation_outside_diameter * pipe_length
    insulation_to_environment_cell_areas = insulation_to_environment_area / n_cells_axial
    
    common_properties = ComponentVector(
        R_gas = 8.314u"J/(mol*K)",

        molecular_weights = (
            methanol = 32.04u"g/mol",
            water = 18.02u"g/mol",
            carbon_monoxide = 28.01u"g/mol",
            hydrogen = 2.02u"g/mol",
            carbon_dioxide = 44.01u"g/mol",
            air = 28.97u"g/mol"
        ),

        #General properties of the pipe
        pipe_length = pipe_length,
        cell_lengths_along_pipe = cell_lengths_along_pipe,
        per_cell_pipe_length = cell_lengths_along_pipe / n_cells_axial,

        reforming_area_outside_diameter = reforming_area_outside_diameter,
        reforming_area_cell_volumes = reforming_area_cell_volumes,
        reforming_area_to_steel_pipe_wall_cell_areas = reforming_area_to_steel_pipe_wall_cell_areas,

        steel_pipe_wall_inside_diameter = steel_pipe_wall_inside_diameter,
        steel_pipe_wall_outside_diameter = steel_pipe_wall_outside_diameter,
        steel_pipe_wall_cell_volumes = steel_pipe_wall_cell_volumes,
        steel_pipe_wall_to_thermocouple_cell_areas = steel_pipe_wall_to_thermocouple_cell_areas,

        thermocouple_region_inside_diameter = thermocouple_region_inside_diameter,
        thermocouple_region_outside_diameter = thermocouple_region_outside_diameter,
        thermocouple_region_cell_volumes = thermocouple_region_cell_volumes,
        thermocouple_region_to_heating_wire_cell_areas = thermocouple_region_to_heating_wire_cell_areas,

        heating_wire_inside_diameter = heating_wire_inside_diameter,
        heating_wire_outside_diameter = heating_wire_outside_diameter,
        heating_wire_cell_volumes = heating_wire_cell_volumes,
        heating_wire_to_insulation_cell_areas = heating_wire_to_insulation_cell_areas,

        insulation_inside_diameter = insulation_inside_diameter,
        insulation_outside_diameter = insulation_outside_diameter,
        insulation_cell_volumes = insulation_cell_volumes,
        insulation_to_environment_cell_areas = insulation_to_environment_cell_areas,
    )

    return common_properties
end

