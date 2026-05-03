using XLSX
using Unitful

#I don't indent anything here so the function can be easily removed later for inspection
function get_trial_1_values_of_note()

experimental_data_path = joinpath(@__DIR__, "PBR_hot_water_test.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

sheet = xf["processed data"]


#the values given by the thermocouples will be used to create an offset for each
#this was observed basically at the start
room_temperature_at_start = sheet["G1"] .* u"°C" .|> u"K"

TC1_temp_at_start = sheet["G2"] .* u"°C" .|> u"K"
TC2_temp_at_start = sheet["G3"] .* u"°C" .|> u"K"
TC3_temp_at_start = sheet["G4"] .* u"°C" .|> u"K"
TC4_temp_at_start = sheet["G5"] .* u"°C" .|> u"K"
TC5_temp_at_start = sheet["G6"] .* u"°C" .|> u"K"

#TC1_offset = room_temperature_at_start - TC1_temp_at_start
#TC2_offset = room_temperature_at_start - TC2_temp_at_start
#TC3_offset = room_temperature_at_start - TC3_temp_at_start
#TC4_offset = room_temperature_at_start - TC4_temp_at_start
#TC5_offset = room_temperature_at_start - TC5_temp_at_start
#these are probably not accurate


multimeter_offset = sheet["G7"] .* u"°C" .|> u"K"
multimeter_offset = Float64(ustrip(multimeter_offset)) * u"K"
#I hate that for some reason unitful converts something like 23°C into something like 5923//20 K, so I always add 0.13 K 

flow_rate_throughout_trial = sheet["G8"] * u"ml/minute" .|> u"m^3/s"
flow_rate_throughout_trial = Float64(ustrip(flow_rate_throughout_trial)) * u"m^3/s"

first_drops_at_outlet_time = sheet["G10"] * u"minute" .|> u"s"

pump_shut_off_time = sheet["G11"] * u"minute" .|> u"s"

room_temp_at_53_minutes = sheet["G12"] * u"°C" .|> u"K"
room_temp_observation_time_2 = 3200u"s"

room_temp_at_103_minutes = sheet["G13"] * u"°C" .|> u"K"
room_temp_observation_time_3 = 6200u"s"

return (
    room_temperature_at_start = room_temperature_at_start,
    TC1_temp_at_start = TC1_temp_at_start,
    TC2_temp_at_start = TC2_temp_at_start,
    TC3_temp_at_start = TC3_temp_at_start,
    TC4_temp_at_start = TC4_temp_at_start,
    TC5_temp_at_start = TC5_temp_at_start,
    TC1_offset = TC1_offset,
    TC2_offset = TC2_offset,
    TC3_offset = TC3_offset,
    TC4_offset = TC4_offset,
    TC5_offset = TC5_offset,
    multimeter_offset = multimeter_offset,
    flow_rate_throughout_trial = flow_rate_throughout_trial,
    first_drops_at_outlet_time = first_drops_at_outlet_time,
    pump_shut_off_time = pump_shut_off_time,
    room_temp_at_53_minutes = room_temp_at_53_minutes,
    room_temp_observation_time_2 = room_temp_observation_time_2,
    room_temp_at_103_minutes = room_temp_at_103_minutes,
    room_temp_observation_time_3 = room_temp_observation_time_3
)

end
