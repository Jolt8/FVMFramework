using DataFrames, CSV
using Unitful

#options for reading data:

#1. "hot_water_flow_tc_temps.csv"
    #this was the trial where I flowed water into the reactor at a known temperature to roughly estimate the thermal mass of the reactor
    #this was also to determinet he overall heat transfer coefficient of the reactor once the pump was shut off at 47 minutes into the trial
    
#2. "heated_trial_tc_temps.csv"
    #this is the trial whre I heated up the reactor periodically to 300*C periodically and then waited for the temperatures to fall back down to a lower level
    #this was to determine the heat transfer coefficient of the heater to the thermocouples as well as the thermal mass of the reactor
    #the temperature of the room during the experiment was 14.5°C which stayed pretty constant throughout the experiment 

function heater_power_interp(t)
    return 0.0
end
    
function get_hot_water_thermocouple_data()

experimental_data_path = joinpath(@__DIR__, "hot_water_flow_tc_temps.csv")
trial = CSV.read(experimental_data_path, DataFrame)

timestamps = (trial.Time_ms .* u"ms") .|> u"s"
timestamps = Float64.(ustrip(timestamps)) * u"s"

TC1_temps_raw = (trial.TC1_C .* u"°C") .|> u"K"
TC2_temps_raw = (trial.TC2_C .* u"°C") .|> u"K"
TC3_temps_raw = (trial.TC3_C .* u"°C") .|> u"K"
TC4_temps_raw = (trial.TC4_C .* u"°C") .|> u"K"
TC5_temps_raw = (trial.TC5_C .* u"°C") .|> u"K"

#screw it, we're using these numbers for creating the offset
room_temp = 20u"°C" .|> u"K"
room_temp = Float64(ustrip(room_temp)) * u"K"

TC1_temp_offset = sum(room_temp .- TC1_temps_raw[1:10])/10
TC2_temp_offset = sum(room_temp .- TC2_temps_raw[1:10])/10
TC3_temp_offset = sum(room_temp .- TC3_temps_raw[1:10])/10
TC4_temp_offset = sum(room_temp .- TC4_temps_raw[1:10])/10
TC5_temp_offset = sum(room_temp .- TC5_temps_raw[1:10])/10

TC1_temps = TC1_temps_raw .+ TC1_temp_offset
TC2_temps = TC2_temps_raw .+ TC2_temp_offset
TC3_temps = TC3_temps_raw .+ TC3_temp_offset
TC4_temps = TC4_temps_raw .+ TC4_temp_offset
TC5_temps = TC5_temps_raw .+ TC5_temp_offset

TC1_temps_interp = LinearInterpolation(ustrip.(vcat(TC1_temps[1], TC1_temps)), ustrip.(vcat(0.0, timestamps)))
TC2_temps_interp = LinearInterpolation(ustrip.(vcat(TC2_temps[1], TC2_temps)), ustrip.(vcat(0.0, timestamps)))
TC3_temps_interp = LinearInterpolation(ustrip.(vcat(TC3_temps[1], TC3_temps)), ustrip.(vcat(0.0, timestamps)))
TC4_temps_interp = LinearInterpolation(ustrip.(vcat(TC4_temps[1], TC4_temps)), ustrip.(vcat(0.0, timestamps)))
TC5_temps_interp = LinearInterpolation(ustrip.(vcat(TC5_temps[1], TC5_temps)), ustrip.(vcat(0.0, timestamps)))

return (
    timestamps = timestamps,
    TC1_temps_interp = TC1_temps_interp,
    TC2_temps_interp = TC2_temps_interp,
    TC3_temps_interp = TC3_temps_interp,
    TC4_temps_interp = TC4_temps_interp,
    TC5_temps_interp = TC5_temps_interp,
    TC1_temp_offset = TC1_temp_offset,
    TC2_temp_offset = TC2_temp_offset,
    TC3_temp_offset = TC3_temp_offset,
    TC4_temp_offset = TC4_temp_offset,
    TC5_temp_offset = TC5_temp_offset,
    heater_power_interp = heater_power_interp
)

end