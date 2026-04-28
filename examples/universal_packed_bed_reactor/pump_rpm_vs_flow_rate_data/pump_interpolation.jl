using XLSX
using Polynomials
using DataFrames, GLM
using Unitful
using Plots

experimental_data_path = joinpath(@__DIR__, "Pump Flow Rate Data.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

trial_sheet_names = ["10 RPM", "5 RPM", "1 RPM", "0.3 RPM", "50 RPM", "25 RPM"]

rpms = Float64[]
average_flow_rates_mL_per_sec = Float64[]
average_flow_rates_mL_per_min = Float64[]

for trial_sheet_name in trial_sheet_names
    ts = xf[trial_sheet_name]
    
    rpm = ts["B1"] 
    push!(rpms, rpm)

    times = vec(ts["A3:A12"]) * u"s"
    grams_accumulated = vec(ts["B3:B12"]) * u"g"

    missing_indices = Int64[]

    for i in eachindex(grams_accumulated)
        if ismissing(grams_accumulated[i])
            push!(missing_indices, i)
        end
    end

    deleteat!(times, missing_indices)
    deleteat!(grams_accumulated, missing_indices)

    seconds_times = Float64.(ustrip.(times))
    minutes_times = seconds_times ./ 60 
    
    grams_accumulated = Float64.(ustrip.(grams_accumulated)) * u"g" #some values were interpreted as Int64

    water_density = 1.0u"g/mL"

    volumes = grams_accumulated ./ water_density

    seconds_data = DataFrame(X = ustrip(seconds_times), Y = ustrip(volumes))
    seconds_model = lm(@formula(Y ~ X), seconds_data)

    minutes_data = DataFrame(X = ustrip(minutes_times), Y = ustrip(volumes))
    minutes_model = lm(@formula(Y ~ X), minutes_data)

    push!(average_flow_rates_mL_per_sec, coef(seconds_model)[2])
    push!(average_flow_rates_mL_per_min, coef(minutes_model)[2])
end

#seconds
rpm_vs_flow_rate_data_seconds = DataFrame(
    RPM = rpms,
    Flow_Rate = average_flow_rates_mL_per_sec
)

plot(rpms, average_flow_rates_mL_per_sec, seriestype = :scatter)

model_seconds = lm(@formula(RPM ~ Flow_Rate), rpm_vs_flow_rate_data_seconds)

coeffs_seconds = coef(model_seconds)

if coeffs_seconds[2] < 0
    equation_seconds = "RPM = $(coeffs_seconds[1]) - $(abs(coeffs_seconds[2])) * Flow_Rate"
else
    equation_seconds = "RPM = $(coeffs_seconds[1]) + $(coeffs_seconds[2]) * Flow_Rate"
end

#minutes
rpm_vs_flow_rate_data_minutes = DataFrame(
    RPM = rpms,
    Flow_Rate = average_flow_rates_mL_per_min
)

plot(rpms, average_flow_rates_mL_per_min, seriestype = :scatter)

model_minutes = lm(@formula(RPM ~ Flow_Rate), rpm_vs_flow_rate_data_minutes)

coeffs_minutes = coef(model_minutes)

if coeffs_minutes[2] < 0
    equation_minutes = "RPM = $(coeffs_minutes[1]) - $(abs(coeffs_minutes[2])) * Flow_Rate"
else
    equation_minutes = "RPM = $(coeffs_minutes[1]) + $(coeffs_minutes[2]) * Flow_Rate"
end

#printed together:

println("linear equation for ml per second: ", equation_seconds)
# RPM = -0.0723492318793547 + 196.6253272135435 * Flow_Rate

println("linear equation for ml per minute: ", equation_minutes)
# RPM = -0.07234923187935079 + 3.2770887868923912 * Flow_Rate