#huh, I guess these are the only properties that are really shared between the different property groups
#I guess I just won't use it

function get_empty_trial_properties()
    empty_trial_properties = ComponentVector(
        pipe_mass_flow = 0.0u"g/minute",

        room_temperature = 14.5u"°C" |> u"K",

        inlet_mass_fractions = ComponentVector(
            #methanol = 1e-20u"kg/kg",
            water = 1e-20u"kg/kg",
            #carbon_monoxide = 1e-20u"kg/kg",
            #hydrogen = 1e-6u"kg/kg",
            #carbon_dioxide = 1e-20u"kg/kg",
            air = 1.0u"kg/kg"
        ),

        empty_mass_fractions = ComponentVector(
            #methanol = 1e-20u"kg/kg",
            water = 1e-20u"kg/kg",
            #carbon_monoxide = 1e-20u"kg/kg",
            #hydrogen = 1e-6u"kg/kg",
            #carbon_dioxide = 1e-20u"kg/kg",
            air = 1.0u"kg/kg"
        ),
    )

    return empty_trial_properties
end