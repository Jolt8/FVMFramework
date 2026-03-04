using Unitful
using CoolProp
using DataInterpolations
using Plots
using Roots

temp = 270.13u"°C" |> u"K"
pressure = 1u"atm" |> u"Pa"

methanol_molar_flow = 0.0005u"mol/s"
water_molar_flow = methanol_molar_flow * 1.3

total_molar_flow = methanol_molar_flow + water_molar_flow

methanol_molar_fraction = methanol_molar_flow / total_molar_flow
water_molar_fraction = water_molar_flow / total_molar_flow

mixture = "methanol[$(methanol_molar_fraction)]&water[$(water_molar_fraction)]"

molar_mass = PropsSI("M", mixture) * u"kg/mol"
density = PropsSI("D", "T", temp, "Q", 1, mixture) 

total_mass_flow = total_molar_flow * molar_mass |> u"g/s"
total_volumetric_flow = total_mass_flow / density |> u"cm^3/s"
pipe_width = 0.5u"cm"
pipe_height = 0.5u"cm"

pipe_area = pipe_width * pipe_height

dynamic_viscosity = PropsSI("V", "T", temp, "Q", 1, mixture) 

pipe_length = 1.0u"cm"

function calculate_Re_from_pipe_area(pipe_area)
    velocity = total_volumetric_flow / pipe_area |> u"m/s"

    Re = (density * velocity * pipe_length) / dynamic_viscosity

    return Re |> u"s/s"
end

calculate_Re_from_pipe_area(pipe_area) 

pipe_area_scan = collect(range(0.0001u"cm^2", 10.0u"cm^2", 1000))
Re_scan = calculate_Re_from_pipe_area.(pipe_area_scan)

p = sortperm(Re_scan)
pipe_area_scan = pipe_area_scan[p]
Re_scan = Re_scan[p]

#plot(pipe_area_scan, calculate_Re_from_pipe_area.(pipe_area_scan))

Re_to_pipe_area = LinearInterpolation(ustrip(pipe_area_scan), Re_scan)

function darcy_weisbach_pressure_drop(pipe_area, pipe_length)
    superficial_mass_velocity = total_mass_flow / pipe_area

    density = PropsSI("D", "T", temp, "Q", 1, mixture)

    gas_viscosity = PropsSI("V", "T", temp, "Q", 1, mixture)

    catalyst_particle_diameter = 1.0u"mm"
    bed_void_fraction = 0.80

    term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (density * catalyst_particle_diameter^2)
    term2 = 1.75 * (1 - bed_void_fraction) / bed_void_fraction^3 * superficial_mass_velocity^2 / (catalyst_particle_diameter * density)
    pressure_drop = -(term1 + term2) * pipe_length

    return pressure_drop |> u"Pa"
end

desired_residence_time = 10u"s"
#5.40u"s"

pipe_width = 1.0u"cm"
pipe_area = pipe_width^2

velocity = total_volumetric_flow / pipe_area |> u"m/s"

gas_density = PropsSI("D", "T", temp, "Q", 1, mixture)
gas_dynamic_viscosity = PropsSI("V", "T", temp, "Q", 1, mixture)

liquid_density = 700u"kg/m^3"
liquid_dynamic_viscosity = 0.00075u"Pa*s"

function desired_pipe_length_pressure_drop_gas_Re_liquid_Re_from_pipe_width(pipe_width)
    pipe_area = pipe_width^2

    velocity = total_volumetric_flow / pipe_area |> u"m/s"

    desired_pipe_length = desired_residence_time * velocity |> u"cm"

    pressure_drop = darcy_weisbach_pressure_drop(pipe_area, desired_pipe_length) 

    gas_Re = (gas_density * velocity * pipe_length) / gas_dynamic_viscosity |> u"kg/kg"

    liquid_Re = (liquid_density * velocity * pipe_length) / liquid_dynamic_viscosity |> u"kg/kg"

    return desired_pipe_length, pressure_drop, gas_Re, liquid_Re
end

n_samples = 100

pipe_width_scan = range(0.1u"cm", 1.5u"cm", n_samples)

desired_pipe_lengths = zeros(n_samples)u"cm"
pressure_drops = zeros(n_samples)u"Pa"
gas_Res = zeros(n_samples)u"kg/kg"
liquid_Res = zeros(n_samples)u"kg/kg"

for i in eachindex(pipe_width_scan)
    desired_pipe_length, pressure_drop, gas_Re, liquid_Re = desired_pipe_length_pressure_drop_gas_Re_liquid_Re_from_pipe_width(pipe_width_scan[i])
    desired_pipe_lengths[i] = desired_pipe_length
    pressure_drops[i] = pressure_drop
    gas_Res[i] = gas_Re
    liquid_Res[i] = liquid_Re
end

pipe_width_scan
desired_pipe_lengths

plot(pipe_width_scan, desired_pipe_lengths, xlabel = "pipe width", ylabel="pipe length")
plot(pipe_width_scan, pressure_drops, xlabel = "pipe width", ylabel="pressure drop")
plot(pipe_width_scan, gas_Res, xlabel = "pipe width", ylabel="gas Re")
plot(pipe_width_scan, liquid_Res, xlabel = "pipe width", ylabel="liquid Re")

desired_pipe_length_pressure_drop_gas_Re_liquid_Re_from_pipe_width(0.35u"cm")


#max_pipe_length = find_zero(pipe_length -> darcy_weisbach_pressure_drop(pipe_area, pipe_length) + max_head_loss, pipe_length)