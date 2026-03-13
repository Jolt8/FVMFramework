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
gas_density = 0.517u"kg/m^3"
liquid_density = 700u"kg/m^3"
#PropsSI("D", "T", temp, "Q", 1, mixture) 

total_mass_flow = total_molar_flow * molar_mass |> u"g/s"
gas_total_volumetric_flow = total_mass_flow / gas_density |> u"cm^3/s"
liquid_total_volumetric_flow = total_mass_flow / liquid_density |> u"cm^3/s"
pipe_width = 0.5u"cm"
pipe_height = 0.5u"cm"

pipe_area = pipe_width * pipe_height

function darcy_weisbach_pressure_drop(pipe_area, pipe_length)
    superficial_mass_velocity = gas_total_volumetric_flow / pipe_area |> u"m/s"

    density = gas_density

    gas_viscosity = PropsSI("V", "T", temp, "Q", 1, mixture)

    catalyst_particle_diameter = 1.0u"mm"
    bed_void_fraction = 0.80

    term1 = 150 * (1 - bed_void_fraction)^2 * gas_viscosity * superficial_mass_velocity / (catalyst_particle_diameter^2 * bed_void_fraction^3)
    term2 = 1.75 * (1 - bed_void_fraction) * density * superficial_mass_velocity^2 / (catalyst_particle_diameter * bed_void_fraction^3)
    pressure_drop = -(term1 + term2) * pipe_length

    return pressure_drop |> u"Pa"
end

darcy_weisbach_pressure_drop(pipe_area, 1.0u"cm")

desired_residence_time = 1.5u"s"

pipe_width = 1.0u"cm"
pipe_area = pipe_width^2

liquid_velocity = liquid_total_volumetric_flow / pipe_area |> u"m/s"
gas_velocity = gas_total_volumetric_flow / pipe_area |> u"m/s"

gas_density = 0.517u"kg/m^3"
#PropsSI("D", "T", temp, "Q", 1, mixture)
gas_dynamic_viscosity = PropsSI("V", "T", temp, "Q", 1, mixture)

liquid_density = 700u"kg/m^3"
liquid_dynamic_viscosity = 1.0e-3u"Pa*s"

function desired_pipe_length_pressure_drop_gas_Re_liquid_Re_from_pipe_width(pipe_width)
    pipe_area = pipe_width^2

    local_gas_velocity = gas_total_volumetric_flow / pipe_area |> u"m/s"
    local_liquid_velocity = liquid_total_volumetric_flow / pipe_area |> u"m/s"

    desired_pipe_length = desired_residence_time * local_gas_velocity |> u"cm"

    pressure_drop = darcy_weisbach_pressure_drop(pipe_area, desired_pipe_length) 

    gas_Re = (gas_density * local_gas_velocity * desired_pipe_length) / gas_dynamic_viscosity |> u"kg/kg"

    liquid_Re = (liquid_density * local_liquid_velocity * desired_pipe_length) / liquid_dynamic_viscosity |> u"kg/kg"

    return desired_pipe_length, pressure_drop, gas_Re, liquid_Re
end

n_samples = 100

pipe_width_scan = range(0.5u"cm", 5.0u"cm", n_samples)

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

desired_pipe_length_pressure_drop_gas_Re_liquid_Re_from_pipe_width(2.0u"cm")


#max_pipe_length = find_zero(pipe_length -> darcy_weisbach_pressure_drop(pipe_area, pipe_length) + max_head_loss, pipe_length)