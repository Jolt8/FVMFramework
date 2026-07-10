using ComponentArrays
using Unitful, UnitfulAssets
using FVMFramework
using NonlinearSolve

stream_names = (:feed, :mix_out, :c1_out, :recycle, :product)
energy_names = (:electricity_1, )

component_names = (:compressor_1, :mixer_1, :splitter_1)

species_symbols = (:carbon_dioxide, )

material_stream_template = ComponentVector(
    mass_flow = 0.0u"kg/s",
    mass_flows = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"kg/s" for _ in eachindex(species_symbols)])...)
    ),
    molar_flow = 0.0u"mol/s",
    molar_flows = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"mol/s" for _ in eachindex(species_symbols)])...)
    ),
    volumetric_flow = 0.0u"m^3/s",
    volumetric_flows = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"m^3/s" for _ in eachindex(species_symbols)])...)
    ),
    temp = 0.0u"K" ,
    pressure = 0.0u"bar" ,
    mass_fractions = (
        ComponentVector(; NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...)
    ),
    molar_fractions = (
        ComponentVector(; NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...)
    ),
    molar_concentrations = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"mol/m^3" for _ in eachindex(species_symbols)])...)
    ),
    specific_heat_ratio = 0.0 ,
    molecular_weights = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"g/mol" for _ in eachindex(species_symbols)])...)
    ),
    mw_avg = 0.0u"g/mol",
    rho = 0.0u"kg/m^3",
    species_cps = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"J/(kg*K)" for _ in eachindex(species_symbols)])...)
    ),
    cp_avg = 0.0u"J/(kg*K)",
    species_cvs = (
        ComponentVector(; NamedTuple{species_symbols}([0.0u"J/(kg*K)" for _ in eachindex(species_symbols)]) ...)
    ),
    cv_avg = 0.0u"J/(kg*K)",

    #constant properties
    R_gas = 8.314u"J/(mol*K)"
)

material_stream_template.molecular_weights.carbon_dioxide = 40.02u"g/mol"
material_stream_template.species_cps.carbon_dioxide = 837.0u"J/(kg*K)"
material_stream_template.species_cvs.carbon_dioxide = 637.0u"J/(kg*K)"

energy_stream_template = ComponentVector(
    #energy stream properties
    electrical_wattage = 0.0u"W",
    heat_wattage = 0.0u"W",
)

streams = ComponentVector(; 
    NamedTuple{stream_names}([deepcopy(material_stream_template) for _ in eachindex(stream_names)])..., 
    NamedTuple{energy_names}([deepcopy(energy_stream_template) for _ in eachindex(energy_names)])...
)

streams.feed.mass_flow = 1.0u"kg/s"
streams.feed.mass_fractions.carbon_dioxide = 1.0
streams.feed.temp = 25.0u"°C"
streams.feed.pressure = 100u"kPa"

components = ComponentVector(
    compressor_1 = (
        pressure_increase = 300u"kPa",
        efficiency = 0.75
    ),
    mixer_1 = 0.0,
    splitter_1 = (
        split_ratios = (
            stream_1 = 0.5,
            stream_2 = 0.5
        ),
    ),
)

#indexed like components.heat_exchanger_1.area

#system will likely become a struct later
system = ComponentVector(
    total_plant_cost = 0.0u"USD",
    upfront_cost = 0.0u"USD",
    interest_rate = 0.0u"USD"
)

abstract type AbstractComponent{CalcMode} end

struct Compressor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

struct TwoMixer{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in_1::Symbol
    stream_in_2::Symbol
    stream_out::Symbol
end

struct TwoSplitter{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out_1::Symbol
    stream_out_2::Symbol
end

# Define your plant layout ONCE, outside the loop
topo = (
    mix_1 = TwoMixer{:AveragedPressure}(:mixer_1, :feed, :recycle, :mix_out),
    compressor_1 = Compressor{:PressureRatio}(:compressor_1, :mix_out, :c1_out, :electricity_1),
    splitter_1 = TwoSplitter{:Ratio}(:splitter_1, :c1_out, :recycle, :product)
)

function copy_stream!(from, to)
    to .= from
end

function run_component(
    type::TwoMixer{:AveragedPressure},
    streams, components, system
)
    stream_in_1 = getproperty(streams, type.stream_in_1)
    stream_in_2 = getproperty(streams, type.stream_in_2)
    stream_out = getproperty(streams, type.stream_out)
    component = getproperty(components, type.id)

    if stream_in_1.mass_flow == 0.0u"kg/s" || stream_in_2.mass_flow == 0.0u"kg/s"
        copy_stream!(stream_in_1, stream_in_2) #this is just temporary until we get a non linear solver in here
    end

    resulting_stream_total_heat_capacity = stream_in_1.cp_avg * stream_in_1.mass_flow + stream_in_2.cp_avg * stream_in_2.mass_flow
    stream_1_capacity_ratio = stream_in_1.cp_avg * stream_in_1.mass_flow / resulting_stream_total_heat_capacity
    stream_2_capacity_ratio = stream_in_2.cp_avg * stream_in_2.mass_flow / resulting_stream_total_heat_capacity

    stream_out.temp = ((stream_in_1.temp |> u"K") * stream_1_capacity_ratio + (stream_in_2.temp |> u"K") * stream_2_capacity_ratio) |> u"°C"
    stream_out.pressure = 0.5*(stream_in_1.pressure + stream_in_2.pressure)
    stream_out.mass_flow = stream_in_1.mass_flow + stream_in_2.mass_flow
    stream_out.mass_fractions = 0.5*(stream_in_1.mass_fractions + stream_in_2.mass_fractions)
    update_stream_properties!(stream_out)
end


function run_component(
    meta::Compressor{:PressureRatio},
    streams, components, system
)
    stream_in = getproperty(streams, meta.stream_in)
    stream_out = getproperty(streams, meta.stream_out)
    energy_in = getproperty(streams, meta.energy_in)
    component = getproperty(components, meta.id)

    copy_stream!(stream_in, stream_out)
    
    stream_out.pressure = stream_in.pressure + component.pressure_increase

    stream_out.temp = (stream_in.temp |> u"K") * (stream_out.pressure / stream_in.pressure) ^ ((stream_in.specific_heat_ratio - 1) / stream_in.specific_heat_ratio)

    ((stream_out.temp |> u"K") - (stream_in.temp |> u"K"))

    ideal_wattage_required = stream_in.molar_flow * stream_in.R_gas * ((stream_out.temp |> u"K") - (stream_in.temp |> u"K")) / (stream_in.specific_heat_ratio - 1.0)

    energy_in.electrical_wattage = ideal_wattage_required / component.efficiency
    update_stream_properties!(stream_out)
end

function run_component(
    type::TwoSplitter{:Ratio},
    streams, components, system
)
    stream_in = getproperty(streams, type.stream_in)
    stream_out_1 = getproperty(streams, type.stream_out_1)
    stream_out_2 = getproperty(streams, type.stream_out_2)
    component = getproperty(components, type.id)

    copy_stream!(stream_in, stream_out_1)
    copy_stream!(stream_in, stream_out_2)

    stream_out_1.mass_flow = stream_in.mass_flow * component.split_ratios.stream_1
    stream_out_2.mass_flow = stream_in.mass_flow * component.split_ratios.stream_2

    update_stream_properties!(stream_out_1)
    update_stream_properties!(stream_out_2)
end

Revise.includet(joinpath(@__DIR__, "new_physics_helper_functions.jl"))

function evaluate_plant_cost(components, streams, system, topo)
    #=
    for stream_name in stream_names
        update_stream_properties!(getproperty(streams, stream_name))
    end

    for stream_name in energy_names
    end
    =#
    update_stream_properties!(streams.feed)

    #TODO: do some graph theory shenanagins to make topo list the components in the right execution order
    foreach(values(topo)) do component
        run_component(component, streams, components, system)
    end

    #this doesn't have to be ordered
    #=
    foreach(values(topo)) do component
        calculate_economics!(component, streams, components, system)
    end
    =#
    
    #return system.total_plant_cost

    return
end

evaluate_plant_cost(components, streams, system, topo)

initial_recycle_guess = [0.5, 300.0, 1.0] 

# 2. Define the Residual Function: F(guess) = Calculated - Guess
function flowsheet_residual!(resid, u_guess, p)
    # p is our (components, streams, system, topo) tuple
    comps, streams, sys, top = p
    
    # A. Inject the solver's guess into the Recycle Stream
    # We apply units here so your internal flowsheet logic stays pure Unitful
    streams.recycle.mass_flow = u_guess[1] * u"kg/s"
    streams.recycle.temp = u_guess[2] * u"K"
    streams.recycle.pressure = u_guess[3] * u"bar"
    
    # (If you had multiple species, you would also inject mass fractions here)
    update_stream_properties!(streams.recycle)

    # B. Run the flowsheet
    # This cascades the guessed recycle through the mixer, compressor, and splitter
    foreach(values(top)) do component
        run_component(component, streams, comps, sys)
    end
    
    # C. Extract the newly calculated recycle state
    calc_mass_flow = ustrip(u"kg/s", streams.recycle.mass_flow)
    calc_temp = ustrip(u"K", streams.recycle.temp)
    calc_pressure = ustrip(u"bar", streams.recycle.pressure)
    
    # D. Calculate the Residuals (We want these to be exactly 0)
    resid[1] = calc_mass_flow - u_guess[1]
    resid[2] = calc_temp - u_guess[2]
    resid[3] = calc_pressure - u_guess[3]
    
    return resid
end

# 3. Setup and Run the Nonlinear Solver
function solve_recycle!(initial_guess, components, streams, system, topo)
    # Package our plant data to pass to the solver
    p = (components, streams, system, topo)
    
    # Create the problem
    prob = NonlinearProblem(flowsheet_residual!, initial_guess, p)
    
    # Solve it! (NewtonRaphson or TrustRegion are incredibly fast for this)
    sol = solve(prob, NewtonRaphson(autodiff = AutoFiniteDiff()), abstol=1e-6)
    
    #println("Recycle Stream Converged!")
    #println("Converged Mass Flow: ", sol.u[1], " kg/s")
    #println("Converged Temp:      ", sol.u[2], " K")
    #println("Converged Pressure:  ", sol.u[3], " bar")
    
    # The flowsheet is now in its perfect steady-state. 
    return sol.u
end

using BenchmarkTools

@btime converged_state = solve_recycle!($initial_recycle_guess, $components, $streams, $system, $topo) 
#=
function test_2(initial_recycle_guess, components, streams, system, topo)
    for _ in 1:1000
        converged_state = solve_recycle!(initial_recycle_guess, components, streams, system, topo) 
    end
    return
end

VSCodeServer.@profview test_2(initial_recycle_guess, components, streams, system, topo)
=#