
using XLSX
using Unitful
using Plots
using DataFrames
using GLM

experimental_data_path = joinpath(@__DIR__, "esterification_data_processing_for_julia.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

typeof(float.(vec(xf["40C"]["A2:A10"])))

struct TrialPreprocessedData
    temperature_timestamps::Vector{Float64}
    temperatures::Vector{Float64}
    moles_timestamps::Vector{Float64}
    moles_respective_temp_idxs::Vector{Int64}
    acetic_acid_moles::Vector{Float64}
    ethanol_moles::Vector{Float64}
    ethyl_acetate_moles::Vector{Float64}
    water_moles::Vector{Float64}
end

T1_temp_end_idx = xf["30C Temps"]["I3"]
T1_temp_timestamps = float.(vec(xf["30C Temps"]["K2:K$T1_temp_end_idx"]))
T1_moles_timestamps = float.(vec(xf["30C"]["A2:A11"]))
T1_moles_respective_temp_idxs = findall(x -> x in T1_moles_timestamps, T1_temp_timestamps)
trial_1_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["30C Temps"]["K2:K$T1_temp_end_idx"])),
    #temperatures
    float.(vec(xf["30C Temps"]["L2:L$T1_temp_end_idx"])),

    #moles_timestamps 
    float.(vec(xf["30C"]["A2:A11"])),
    #moles_respective_temp_idxs
    T1_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["30C"]["G2:G11"])),
    #ethanol_moles
    float.(vec(xf["30C"]["H2:H11"])),
    #ethyl_acetate_moles
    float.(vec(xf["30C"]["J2:J11"])),
    #water_moles
    float.(vec(xf["30C"]["I2:I11"])),
)


T2_temp_end_idx = xf["40C Temps"]["I3"]
T2_temp_timestamps = float.(vec(xf["40C Temps"]["K2:K$T2_temp_end_idx"]))
T2_moles_timestamps = float.(vec(xf["40C"]["A2:A11"]))
T2_moles_respective_temp_idxs = findall(x -> x in T2_moles_timestamps, T2_temp_timestamps)
trial_2_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["40C Temps"]["K2:K$T2_temp_end_idx"])),
    #temperatures
    float.(vec(xf["40C Temps"]["L2:L$T2_temp_end_idx"])),

    #moles_timestamps 
    float.(vec(xf["40C"]["A2:A11"])),
    #moles_respective_temp_idxs
    T2_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["40C"]["E2:E11"])),
    #ethanol_moles
    float.(vec(xf["40C"]["F2:F11"])),
    #ethyl_acetate_moles
    float.(vec(xf["40C"]["H2:H11"])),
    #water_moles
    float.(vec(xf["40C"]["G2:G11"])),
)

T3_temp_end_idx = xf["50C Temps"]["I3"]
T3_temp_timestamps = float.(vec(xf["50C Temps"]["K2:K$T3_temp_end_idx"]))
T3_moles_timestamps = float.(vec(xf["50C"]["A2:A14"]))
T3_moles_respective_temp_idxs = findall(x -> x in T3_moles_timestamps, T3_temp_timestamps)
trial_3_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["50C Temps"]["K2:K$T3_temp_end_idx"])),
    #temperatures
    float.(vec(xf["50C Temps"]["L2:L$T3_temp_end_idx"])),

    #moles_timestamps
    float.(vec(xf["50C"]["A2:A14"])),
    #moles_respective_temp_idxs
    T3_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["50C"]["G2:G14"])),
    #ethanol_moles
    float.(vec(xf["50C"]["H2:H14"])),
    #ethyl_acetate_moles
    float.(vec(xf["50C"]["J2:J14"])),
    #water_moles
    float.(vec(xf["50C"]["I2:I14"])),
)

all_trials_data = [trial_1_data, trial_2_data, trial_3_data]

struct Trial
    temperature_timestamps::Vector{typeof(1.0u"s")}
    temperatures::Vector{typeof(1.0u"K")}

    moles_timestamps::Vector{typeof(1.0u"s")}
    moles_respective_temp_idxs::Vector{Int64}
    molar_concentrations_matrix::Array{typeof(1.0u"mol/L"), 2}
end

function get_mass_fraction(species_moles, species_molecular_weight, rho, mixture_volume)
    (((species_moles / mixture_volume) * species_molecular_weight) / rho)
end

trials = Trial[]

all_trials_data[1].temperature_timestamps
all_trials_data[1].temperatures
all_trials_data[1].moles_timestamps
all_trials_data[1].moles_respective_temp_idxs
all_trials_data[1].acetic_acid_moles
all_trials_data[1].ethanol_moles
all_trials_data[1].ethyl_acetate_moles
all_trials_data[1].water_moles

for trial in all_trials_data
    temp_timestamps = trial.temperature_timestamps .* u"s"

    temperatures_deg_C = trial.temperatures .* u"°C"

    trial_temperatures = temperatures_deg_C .|> u"K"

    moles_timestamps = trial.moles_timestamps .* u"s"
    moles_respective_temp_idxs = trial.moles_respective_temp_idxs

    acetic_acid_moles = trial.acetic_acid_moles .* u"mol"
    ethanol_moles = trial.ethanol_moles .* u"mol"
    ethyl_acetate_moles = trial.ethyl_acetate_moles .* u"mol"
    water_moles = trial.water_moles .* u"mol"

    mixture_volume = 50u"ml"

    acetic_acid_molar_concentrations = acetic_acid_moles ./ mixture_volume
    ethanol_molar_concentrations = ethanol_moles ./ mixture_volume
    ethyl_acetate_molar_concentrations = ethyl_acetate_moles ./ mixture_volume
    water_molar_concentrations = water_moles ./ mixture_volume

    molar_concentrations_matrix = zeros(typeof(1.0u"mol/L"), length(moles_timestamps), 4)

    for i in eachindex(moles_timestamps)
        molar_concentrations_matrix[i, :] .= [acetic_acid_molar_concentrations[i], ethanol_molar_concentrations[i], ethyl_acetate_molar_concentrations[i], water_molar_concentrations[i]]
    end
    #molar_concentrations_matrix is formatted as [time_idx, species_idx]

    temperatures = zeros(typeof(1.0u"K"), length(temp_timestamps))

    for i in eachindex(temp_timestamps)
        temperatures[i] += trial_temperatures[i]
    end
    #temperatures is formatted as [time_idx]

    trial = Trial(temp_timestamps, temperatures, moles_timestamps, moles_respective_temp_idxs, molar_concentrations_matrix)

    push!(trials, trial)
end

trials[1].moles_timestamps
trials[1].moles_respective_temp_idxs
trials[1].molar_concentrations_matrix
trials[1].temperatures
trials[1].temperature_timestamps

trials[1].temperatures[trials[1].moles_respective_temp_idxs[1]]

#arrhenius fitting
struct TrialArrheniusData
    average_temperatures::Vector{typeof(1.0u"K")}
    rate_constants::Vector{typeof(1.0u"L/(mol*s)")}
end

all_trials_arrhenius_data = TrialArrheniusData[]

for trial in trials
    n_timestamps = length(trial.moles_timestamps)
    trial_arrhenius_data = TrialArrheniusData(zeros(typeof(1.0u"K"), n_timestamps-1), zeros(typeof(1.0u"L/(mol*s)"), n_timestamps-1))
    
    for i in 1:length(trial.moles_respective_temp_idxs)-1
        temp_start_idx = trial.moles_respective_temp_idxs[i]
        temp_end_idx = trial.moles_respective_temp_idxs[i+1]
        
        avg_temp = sum(trial.temperatures[temp_start_idx:temp_end_idx]) / (length(trial.temperatures[temp_start_idx:temp_end_idx]))
    
        trial_arrhenius_data.average_temperatures[i] = avg_temp
        trial_arrhenius_data.rate_constants[i] = (1/(trial.molar_concentrations_matrix[i+1, 1]) - 1/(trial.molar_concentrations_matrix[i, 1])) / (trial.moles_timestamps[i+1] - trial.moles_timestamps[i])
    end
    push!(all_trials_arrhenius_data, trial_arrhenius_data)
end

all_trials_arrhenius_data[1].average_temperatures
all_trials_arrhenius_data[1].rate_constants

all_trials_arrhenius_data[2].average_temperatures
all_trials_arrhenius_data[2].rate_constants

all_trials_arrhenius_data[3].average_temperatures
all_trials_arrhenius_data[3].rate_constants

one_over_temperatures = Float64[]
ln_rate_constants = Float64[]

for trial in all_trials_arrhenius_data
    for (temperature, rate_constant) in zip(trial.average_temperatures, trial.rate_constants)
        push!(one_over_temperatures, 1.0/ustrip(uconvert(u"K", temperature)))
        push!(ln_rate_constants, log(ustrip(rate_constant)))
    end
end

one_over_temperatures 
ln_rate_constants

plot(one_over_temperatures, ln_rate_constants,
    seriestype=:scatter,
    xlabel="1/T (K⁻¹)",
    ylabel="ln(k)",
    title="Arrhenius Plot for Esterification",
    legend=false,
    framestyle=:box,
    #xlims = (0, one_over_temperatures[end]) #if you want to see how you would get that y intercept
)

data = DataFrame(
    one_over_T = one_over_temperatures,
    ln_k = ln_rate_constants,
)

model = lm(@formula(ln_k ~ one_over_T), data)

coefficients = coef(model)

y_intercept = coefficients[1]  #this is ln(pre exponential factor)
pre_exponential_factor = exp(y_intercept) #units of L/(mol*s)

slope = coefficients[2] #this is our Ea/R_GAS
Ea = -slope * 8.314 #comes out as 38209.4 J/mol

