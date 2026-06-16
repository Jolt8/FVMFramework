function get_hot_water_trial_properties(values_of_note)
    water_density = 998u"kg/m^3"
    pipe_mass_flow = values_of_note.flow_rate_throughout_trial * water_density

    inlet_mass_fractions = ComponentVector(
        #methanol = 1e-20u"kg/kg",
        water = 1.0u"kg/kg",
        #carbon_monoxide = 0.0u"kg/kg",
        #hydrogen = 0.0u"kg/kg",
        #carbon_dioxide = 0.0u"kg/kg",
        air = 0.0u"kg/kg"
    )

    total_inlet_mass_fractions = sum(inlet_mass_fractions)
    inlet_mass_fractions = inlet_mass_fractions ./ total_inlet_mass_fractions

    empty_mass_fractions = ComponentVector(
        #methanol = 1e-20u"kg/kg",
        water = 1e-20u"kg/kg",
        #carbon_monoxide = 1e-20u"kg/kg",
        #hydrogen = 1e-6u"kg/kg",
        #carbon_dioxide = 1e-20u"kg/kg",
        air = 1.0u"kg/kg"
    )

    total_empty_mass_fractions = sum(empty_mass_fractions)
    empty_mass_fractions = empty_mass_fractions ./ total_empty_mass_fractions

    hot_water_trial_properties = ComponentVector(
        pipe_mass_flow = pipe_mass_flow,

        room_temperature = 20u"°C" |> u"K",

        inlet_mass_fractions = inlet_mass_fractions,
        empty_mass_fractions = empty_mass_fractions
    )

    return hot_water_trial_properties
end
