using DataFrames, CSV
using Unitful

#options for reading data:

#1. "reactor_data_2026_05_02__15_56_13.csv"
    #this was the trial where I flowed water into the reactor at a known temperature to roughly estimate the thermal mass of the reactor
    #this was also to determinet he overall heat transfer coefficient of the reactor once the pump was shut off at 47 minutes into the trial

function get_thermocouple_data(experimental_data_path::String)

trial = CSV.read(experimental_data_path, DataFrame)

timestamps = (trial.Time_ms .* u"ms") .|> u"s"
timestamps = Float64.(ustrip(timestamps)) * u"s"

temp_tc1 = (trial.TC1_C .* u"°C") .|> u"K"
temp_tc2 = (trial.TC2_C .* u"°C") .|> u"K"
temp_tc3 = (trial.TC3_C .* u"°C") .|> u"K"
temp_tc4 = (trial.TC4_C .* u"°C") .|> u"K"
temp_tc5 = (trial.TC5_C .* u"°C") .|> u"K"

#screw it, we're using these numbers for creating the offset
room_temp = 20u"°C" .|> u"K"
room_temp = Float64(ustrip(room_temp)) * u"K"

temp_tc1_offset = sum(room_temp .- temp_tc1[1:10])/10
temp_tc2_offset = sum(room_temp .- temp_tc2[1:10])/10
temp_tc3_offset = sum(room_temp .- temp_tc3[1:10])/10
temp_tc4_offset = sum(room_temp .- temp_tc4[1:10])/10
temp_tc5_offset = sum(room_temp .- temp_tc5[1:10])/10

temp_tc1 .+= temp_tc1_offset
temp_tc2 .+= temp_tc2_offset
temp_tc3 .+= temp_tc3_offset
temp_tc4 .+= temp_tc4_offset
temp_tc5 .+= temp_tc5_offset

return (
    timestamps = timestamps,
    temp_tc1 = temp_tc1,
    temp_tc2 = temp_tc2,
    temp_tc3 = temp_tc3,
    temp_tc4 = temp_tc4,
    temp_tc5 = temp_tc5,
    temp_tc1_offset = temp_tc1_offset,
    temp_tc2_offset = temp_tc2_offset,
    temp_tc3_offset = temp_tc3_offset,
    temp_tc4_offset = temp_tc4_offset,
    temp_tc5_offset = temp_tc5_offset,
)

end