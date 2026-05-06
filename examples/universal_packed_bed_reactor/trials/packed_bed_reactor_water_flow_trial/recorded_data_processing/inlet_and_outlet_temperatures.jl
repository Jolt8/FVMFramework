using XLSX
using Unitful
using DataFrames, CSV, GLM
#using LsqFit
using DataInterpolations
using Plots

function get_inlet_and_outlet_temperature_correlations()

#function get_inlet_and_outlet_correlations
experimental_data_path = joinpath(@__DIR__, "PBR_hot_water_test.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

timestamps = (vec(xf["processed data"]["B2:B96"]) .* u"minute") .|> u"s"

#Inlet Temperatures
#note, for some reason my multimeter was reading obscenely low values from the thermocouple, so I have to add a value that I found to get it right
#also turn it negative because I didn't record that it was negative at all
inlet_temperatures_unprocessed = vec(xf["processed data"]["C2:C96"])
inlet_temperatures_with_missing = ((inlet_temperatures_unprocessed .* -1) .+ 160.13) .* u"°C" .|> u"K"

inlet_timestamps = []
inlet_temperatures = []

for i in eachindex(inlet_temperatures_with_missing)
    if ismissing(inlet_temperatures_with_missing[i]) == false && timestamps[i] > 540u"s"
        #before 540 seconds I was taking temperatures at the wrong spots, I'm going to have to interpolate the temperatures now
        #hopefully a linear correlation is enough since it looks pretty linear from the plot although that's not how heat transfer works
        push!(inlet_timestamps, timestamps[i])
        push!(inlet_temperatures, inlet_temperatures_with_missing[i])
    end
end

inlet_timestamps_and_temperatures = DataFrame(
    inlet_timestamps = ustrip.(inlet_timestamps),
    inlet_temperatures = ustrip.(inlet_temperatures)
)

model_inlet_temperatures = lm(@formula(inlet_temperatures ~ inlet_timestamps), inlet_timestamps_and_temperatures)

coeffs_inlet_temperatures = coef(model_inlet_temperatures)

if coeffs_inlet_temperatures[2] < 0
    equation_inlet_temperatures = "inlet_temperatures = $(coeffs_inlet_temperatures[1]) - $(abs(coeffs_inlet_temperatures[2])) * inlet_timestamps"
else
    equation_inlet_temperatures = "inlet_temperatures = $(coeffs_inlet_temperatures[1]) + $(coeffs_inlet_temperatures[2]) * inlet_timestamps"
end

interpolated_inlet_temperatures = [coeffs_inlet_temperatures[1] + coeffs_inlet_temperatures[2] * t for t in ustrip.(timestamps)] .* u"K"

plot(ustrip.(inlet_timestamps), ustrip.(inlet_temperatures))
plot!(ustrip.(timestamps), ustrip.(interpolated_inlet_temperatures))

inlet_temperatures_interpolation = CubicSpline(ustrip.(interpolated_inlet_temperatures), ustrip.(timestamps))


#Outlet Temperatures
outlet_temperatures_unprocessed = vec(xf["processed data"]["D2:D96"])
outlet_temperatures_with_missing = ((outlet_temperatures_unprocessed .* -1) .+ 160.13) .* u"°C" .|> u"K"

outlet_timestamps = []
outlet_temperatures = []

for i in eachindex(outlet_temperatures_with_missing)
    if ismissing(outlet_temperatures_with_missing[i]) == false
        push!(outlet_timestamps, timestamps[i])
        push!(outlet_temperatures, outlet_temperatures_with_missing[i])
    end
end

plot(ustrip.(outlet_timestamps), ustrip.(outlet_temperatures))

#this one doesn't need to be interpolated because the curve is much more complex
#note that the first droplets of this liquid came out at around 5 minutes in
#thus, we're going to have to account for that in the fit

outlet_temperatures_interpolation = CubicSpline(ustrip.(vcat(outlet_temperatures[1], outlet_temperatures)), ustrip.(vcat(timestamps[1], outlet_timestamps)))

plot(ustrip.(outlet_timestamps), ustrip.(outlet_temperatures))

small_interpolation_timestamps = range(timestamps[1], outlet_timestamps[end], 500)
plot!(ustrip.(small_interpolation_timestamps), outlet_temperatures_interpolation.(ustrip.(small_interpolation_timestamps))) 
#yep, this looks good!\

return (
    inlet_temp_interp = inlet_temperatures_interpolation,
    outlet_temp_interp = outlet_temperatures_interpolation
)

end