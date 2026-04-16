using CoolProp
using Unitful

#mixture = "methanol[$(methanol_molar_fraction)]&water[$(water_molar_fraction)]"

mixture = "water"

T_in = 21u"°C"
pressure = 1u"atm"

mixture_volumetric_flow = 1u"ml/minute"

rho = PropsSI("D", "T", T_in, "P", pressure, mixture)

mass_flow = mixture_volumetric_flow * rho |> u"kg/s"

T_sat = PropsSI("T", "P", pressure, "Q", 1, mixture)

delta_H = PropsSI("H", "T", T_sat, "Q", 1, mixture) - PropsSI("H", "T", T_in, "Q", 0, mixture)

required_heater_section_wattage = delta_H * mass_flow |> u"W"

heater_approximate_efficiency = 0.95

required_heater_section_wattage /= heater_approximate_efficiency

heater_section_percentage_of_total_length = 3u"inch" / 12u"inch"

total_reactor_wattage = required_heater_section_wattage / heater_section_percentage_of_total_length

wire_resistance = 32.0u"Ω"

required_voltage = sqrt(total_reactor_wattage * wire_resistance) |> u"V"