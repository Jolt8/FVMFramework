using ComponentArrays
using Unitful, UnitfulAssets
using FVMFramework
using NonlinearSolve
using PreallocationTools
using ForwardDiff
using ImplicitDifferentiation

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

streams_vec = ustrip.(upreferred.(Vector(streams)))
components_vec = ustrip.(upreferred.(Vector(components)))

streams_axes = virtual_merge_axes((streams,))
components_axes = virtual_merge_axes((components,))

initial_recycle_guess = ustrip.(upreferred.([0.5u"kg/s", 300.0u"K", 1.0u"bar"]))
N = ForwardDiff.pickchunksize(length(initial_recycle_guess))
streams_diff_cache = DiffCache(streams_vec, N)

system = (
    streams_diff_cache = streams_diff_cache,
    streams_vec = streams_vec,
    streams_virtual_axes = streams_axes,
    components_vec = components_vec,
    components_virtual_axes = components_axes,
    topo = topo,
    total_plant_cost = 0.0,
    upfront_cost = 0.0,
    interest_rate = 0.0
)

function unpack_process_state(u_guess, sys)
    streams_tmp = get_tmp(sys.streams_diff_cache, first(u_guess))
    streams_tmp .= sys.streams_vec
    
    streams_virt = VirtualFVMArray((streams_tmp,), sys.streams_virtual_axes)
    comps_virt = VirtualFVMArray((sys.components_vec,), sys.components_virtual_axes)
    
    return streams_virt, comps_virt
end

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

    resulting_stream_total_heat_capacity = stream_in_1.cp_avg * stream_in_1.mass_flow + stream_in_2.cp_avg * stream_in_2.mass_flow
    stream_1_capacity_ratio = stream_in_1.cp_avg * stream_in_1.mass_flow / resulting_stream_total_heat_capacity
    stream_2_capacity_ratio = stream_in_2.cp_avg * stream_in_2.mass_flow / resulting_stream_total_heat_capacity

    stream_out.temp = stream_in_1.temp * stream_1_capacity_ratio + stream_in_2.temp * stream_2_capacity_ratio
    stream_out.pressure = 0.5*(stream_in_1.pressure + stream_in_2.pressure)
    stream_out.mass_flow = stream_in_1.mass_flow + stream_in_2.mass_flow
    for_fields!(stream_out.mass_fractions, stream_in_1.mass_fractions, stream_in_2.mass_fractions) do species, stream_out_mass_fractions, stream_in_1_mass_fractions, stream_in_2_mass_fractions
        stream_out_mass_fractions[species] = 0.5 * (stream_in_1_mass_fractions[species] + stream_in_2_mass_fractions[species])
    end
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

    stream_out.temp = stream_in.temp * (stream_out.pressure / stream_in.pressure) ^ ((stream_in.specific_heat_ratio - 1) / stream_in.specific_heat_ratio)

    (stream_out.temp - stream_in.temp)

    ideal_wattage_required = stream_in.molar_flow * stream_in.R_gas * (stream_out.temp - stream_in.temp) / (stream_in.specific_heat_ratio - 1.0)

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

#evaluate_plant_cost(components, streams, system, topo)

# 2. Define the Residual Function: F(guess) = Calculated - Guess
function flowsheet_residual!(resid, u_guess, sys)
    streams_virt, comps_virt = unpack_process_state(u_guess, sys)
    
    # A. Inject the solver's guess into the Recycle Stream
    # Operating purely in unitless base SI
    streams_virt.recycle.mass_flow = u_guess[1]
    streams_virt.recycle.temp = u_guess[2]
    streams_virt.recycle.pressure = u_guess[3]
    
    # (If you had multiple species, you would also inject mass fractions here)
    update_stream_properties!(streams_virt.recycle)

    # B. Run the flowsheet
    # This cascades the guessed recycle through the mixer, compressor, and splitter
    foreach(values(sys.topo)) do component
        run_component(component, streams_virt, comps_virt, sys)
    end
    
    # C. Extract the newly calculated recycle state
    calc_mass_flow = streams_virt.recycle.mass_flow
    calc_temp = streams_virt.recycle.temp
    calc_pressure = streams_virt.recycle.pressure
    
    # D. Calculate the Residuals (We want these to be exactly 0)
    resid[1] = calc_mass_flow - u_guess[1]
    resid[2] = calc_temp - u_guess[2]
    resid[3] = calc_pressure - u_guess[3]
    
    return resid
end

# 3. Setup and Run the Nonlinear Solver
function solve_recycle!(initial_guess, sys)
    # Package our plant data to pass to the solver
    p = sys
    
    # Create the problem
    prob = NonlinearProblem(flowsheet_residual!, initial_guess, p)
    
    # Solve it! (NewtonRaphson or TrustRegion are incredibly fast for this)
    sol = solve(prob, NewtonRaphson(), abstol=1e-6)
    
    #println("Recycle Stream Converged!")
    #println("Converged Mass Flow: ", sol.u[1], " kg/s")
    #println("Converged Temp:      ", sol.u[2], " K")
    #println("Converged Pressure:  ", sol.u[3], " bar")
    
    # The flowsheet is now in its perfect steady-state. 
    return sol.u
end

using BenchmarkTools

@btime converged_state = solve_recycle!($initial_recycle_guess, $system) #84.600 μs (1422 allocations: 73.11 KiB)!

implicit_recycle_solver = ImplicitFunction(solve_recycle!, flowsheet_residual!)

@btime converged_u = implicit_recycle_solver($initial_recycle_guess, $system) #84.900 μs (1423 allocations: 73.41 KiB)!



#=
function test_2(initial_recycle_guess, system)
    for _ in 1:1000
        converged_state = solve_recycle!(initial_recycle_guess, system) 
    end
    return
end

VSCodeServer.@profview test_2(initial_recycle_guess, system)
=#
