using CoolProp
using Unitful
using ComponentArrays
using OrdinaryDiffEq
using Plots

function Nu_packed_bed_Gnielinski(Re, Pr, bed_void_fraction)
    #borrowed from ht.py (love that package)
    #Calculates Nusselt number of a fluid passing over a bed of particles

    Nu_lam = 0.664 * Re^0.5 * Pr^(1/3.)
    
    Nu_turb = 0.037 * Re^0.8 * Pr/(1 + 2.443 * Re^-0.1 * (Pr^(2/3.) - 1))
    
    Nu_sphere = 2 + (Nu_lam^2 + Nu_turb^2)^0.5
    
    fa = 1.0 + 1.5 * (1.0 - bed_void_fraction)

    Nu = Nu_sphere * fa
    
    return Nu
end

function calc_Re(rho, velocity, characteristic_dimension, viscosity)
    return (rho * velocity * characteristic_dimension) / viscosity
end

function calc_Pr(viscosity, cp, thermal_conductivity)
    return (viscosity * cp) / thermal_conductivity
end

methanol_mw = 32.04u"g/mol"
water_mw = 18.02u"g/mol"

methanol_density = 0.791u"g/cm^3"

external_temperature = 275u"°C"
pressure = 1u"bar"

methanol_volumetric_flow = 1u"ml/minute"

methanol_molar_flow = (methanol_volumetric_flow * methanol_density) / methanol_mw |> u"mol/s"
water_molar_flow = methanol_molar_flow * 1.3 |> u"mol/s"

total_molar_flow = methanol_molar_flow + water_molar_flow 

methanol_molar_fraction = methanol_molar_flow / total_molar_flow
water_molar_fraction = water_molar_flow / total_molar_flow

avg_mw = methanol_molar_fraction * methanol_mw + water_molar_fraction * water_mw

mixture = "methanol[$(methanol_molar_fraction)]&water[$(water_molar_fraction)]"

mixture_mw = PropsSI("M", mixture) * u"kg/mol"
rho = PropsSI("D", "T", external_temperature, "Q", 0, mixture)
total_volumetric_flow = (total_molar_flow * mixture_mw / rho) |> u"ml/s"
total_mass_flow = total_molar_flow * mixture_mw |> u"kg/s"

pipe_outer_diameter = 0.5u"inch"
pipe_thickness = 2.0u"mm"
pipe_inner_diameter = pipe_outer_diameter - 2 * pipe_thickness

pipe_area = pi * (pipe_inner_diameter / 2)^2

superficial_velocity = total_volumetric_flow / pipe_area |> u"m/s"

bed_void_fraction = 0.4
interstitial_velocity = superficial_velocity / bed_void_fraction

packed_media_diameter = 0.5u"mm"

viscosity = PropsSI("V", "T", external_temperature, "Q", 0, mixture)
cp = PropsSI("C", "T", external_temperature, "Q", 0, mixture)
k = PropsSI("conductivity", "T", external_temperature, "Q", 0, mixture)

Re = calc_Re(rho, superficial_velocity, packed_media_diameter, viscosity) |> u"kg/kg"
Pr = calc_Pr(viscosity, cp, k)

Nu_conv = Nu_packed_bed_Gnielinski(Re, Pr, bed_void_fraction)

heat_transfer_coeff = (Nu_conv * k) / packed_media_diameter |> u"W/(m^2*K)"

pipe_length = 3.0u"inch"
pipe_inner_surface_area = pi * pipe_inner_diameter * pipe_length |> u"cm^2"
pipe_outer_surface_area = pi * pipe_outer_diameter * pipe_length |> u"cm^2"

convective_resistance = 1 / (heat_transfer_coeff * pipe_inner_surface_area) |> u"K/W"

pipe_k = 20u"W/(m*K)"

total_cond_resistance = log(pipe_outer_diameter / pipe_inner_diameter) / (2 * pi * pipe_k * pipe_length)

UA_prime = 1 / (convective_resistance + total_cond_resistance) |> u"W/K"

heat_transfer_per_m_per_K = UA_prime / pipe_length |> u"W/(m*K)"

inlet_temp = 21.0u"°C"
desired_temp = 275.0u"°C"

delta_H_to_desired = PropsSI("H", "T", desired_temp, "Q", 1, mixture) - PropsSI("H", "T", inlet_temp, "Q", 0, mixture)
T_sat = PropsSI("T", "P", pressure, "Q", 0, mixture)
latent_heat = PropsSI("H", "T", T_sat, "Q", 1, mixture) - PropsSI("H", "T", T_sat, "Q", 0, mixture)

watts_for_desired = delta_H_to_desired * total_mass_flow |> u"W"

heat_transfer_per_m_per_K * pipe_length |> u"W/K"

p = ComponentVector(
    pipe_length = pipe_length,
    heat_transfer_per_m_per_K = heat_transfer_per_m_per_K,
    outside_temp = external_temperature,
    T_sat = T_sat,
    total_mass_flow = total_mass_flow,
    liquid_cp = PropsSI("C", "T", external_temperature, "Q", 0, mixture),
    vapor_cp = PropsSI("C", "T", external_temperature, "Q", 1, mixture),
    latent_heat = latent_heat
)

u = ComponentVector(
    fluid_temperature = 21.0u"°C" |> u"K",
    fluid_quality = 0.0,
)
u_vec = Vector(u)
u_axes = getaxes(u)

function ode_system!(du, u, p, z, u_axes)
    u = ComponentVector(u, u_axes)
    du .*= 0.0
    du = ComponentVector(du, u_axes)
    #println(u)
    
    du_heat = p.heat_transfer_per_m_per_K * ((p.outside_temp |> u"K") - (u.fluid_temperature |> u"K"))
    
    if u.fluid_temperature > p.T_sat && u.fluid_quality < 0.99
        du.fluid_quality = du_heat / (p.total_mass_flow * p.latent_heat)
    elseif u.fluid_quality <= 0.01
        #p.liquid_cp = PropsSI("C", "T", u.fluid_temperature, "Q", 0, mixture)
        du.fluid_temperature = du_heat / (p.total_mass_flow * p.liquid_cp)
    else
        #p.vapor_cp = PropsSI("C", "T", u.fluid_temperature, "Q", 1, mixture)
        du.fluid_temperature = du_heat / (p.total_mass_flow * p.vapor_cp)
    end

    #=if u.fluid_quality > 1.0 #conndense if it goes above
        du_heat = (u.fluid_quality - 1.0) * (p.total_mass_flow * p.latent_heat)
        du.fluid_temperature += du_heat / (p.total_mass_flow * p.cp)
    end=#
end

ode_closure = (du, u, p, z) -> ode_system!(du, u, p, z, u_axes)

prob = ODEProblem(ode_closure, u_vec, (0.0u"m", 12u"inch"), p)
sol = solve(prob, Tsit5(), dt = 12u"inch" / 1000, adaptive = false)
#(12.0u"inch" / 100))

fluid_temperatures = [ComponentVector(sol.u[i], u_axes).fluid_temperature |> u"°C" for i in eachindex(sol.t)]
fluid_qualities = [ComponentVector(sol.u[i], u_axes).fluid_quality for i in eachindex(sol.t)]

good_enough_temp = 260u"°C"

function find_good_enough_length(fluid_temperatures, good_enough_temp)
    min_difference = 10000u"K"
    min_idx = 1
    for i in eachindex(fluid_temperatures)
        if abs(fluid_temperatures[i] - good_enough_temp) < min_difference
            min_difference = abs(fluid_temperatures[i] - good_enough_temp)
            min_idx = i
        end
    end

    return min_idx
end

min_idx = find_good_enough_length(fluid_temperatures, good_enough_temp)
sol.t[min_idx] |> u"inch"

plot(sol.t .|> u"inch", fluid_temperatures, xlabel = "Length", ylabel = "Temperature")
plot(sol.t .|> u"inch", fluid_qualities, xlabel = "Length", ylabel = "Quality")

#Hmm, 6 inches is quite a lot I only have 6 inches of reactor to work with. Is there an experiment that I could perform beforehand to create a correlation to not need a pre-heater?