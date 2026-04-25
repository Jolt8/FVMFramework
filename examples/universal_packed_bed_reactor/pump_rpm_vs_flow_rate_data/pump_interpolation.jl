using XLSX
using Polynomials
using DataFrames, GLM
using Unitful
using Plots

experimental_data_path = joinpath(@__DIR__, "Pump Flow Rate Data.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

trial_sheet_names = ["10 RPM", "5 RPM", "1 RPM", "0.3 RPM", "50 RPM", "25 RPM"]

rpms = Float64[]
average_flow_rates = Float64[]

for trial_sheet_name in trial_sheet_names
    ts = xf[trial_sheet_name]
    
    rpm = ts["B1"] 
    push!(rpms, rpm)

    times = vec(ts["A3:A12"]) * u"s"
    grams_accumulated = vec(ts["B3:B12"]) * u"g"

    missing_indicies = []

    for i in eachindex(grams_accumulated)
        if ismissing(grams_accumulated[i])
            push!(missing_indicies, i)
        end
    end

    deleteat!(times, missing_indicies)
    deleteat!(grams_accumulated, missing_indicies)

    times = Float64.(ustrip.(times)) * u"s"
    grams_accumulated = Float64.(ustrip.(grams_accumulated)) * u"g" #some values were interpreted as Int64

    water_density = 1.0u"g/mL"

    volumes = grams_accumulated ./ water_density

    data = DataFrame(X = ustrip(times), Y = ustrip(volumes))

    model = lm(@formula(Y ~ X), data)

    push!(average_flow_rates, coef(model)[2])
end

rpm_vs_flow_rate_data = DataFrame(
    RPM = rpms,
    Flow_Rate = average_flow_rates
)

plot(rpms, flow_rates, seriestype = :scatter)

model = lm(@formula(Flow_Rate ~ RPM), rpm_vs_flow_rate_data)

coeffs = coef(model)

if coeffs[2] < 0
    equation = "Flow_Rate = $(coeffs[1]) - $(abs(coeffs[2])) * RPM"
else
    equation = "Flow_Rate = $(coeffs[1]) + $(coeffs[2]) * RPM"
end

println(equation)