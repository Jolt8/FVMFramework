using ComponentArrays
using Unitful, UnitfulAssets

stream_names = (:feed, :mix_out, :c1_out, :recycle, :product)
energy_names = (:electricity_1)

component_names = (:compressor_1, :mixer_1, :splitter_1)

species_symbols = (:carbon_dioxide)

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
stream.feed.temp = 25.0u"°C"
stream.feed.pressure = 100u"kPa"

components = ComponentVector(
    compressor_1 = (
        pressure_ratio = 3.0,
    ),
    mixer_1 = 0.0,
    splitter_1 = (
        split_ratios = (
            stream_1 = 0.5,
            stream_2 = 0.5
        )
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

stream_names = (:feed, :mix_out, :c1_out, :recycle, :product)
energy_names = (:electricity_1)

component_names = (:compressor_1, mixer_1, splitter_1)

# Define your plant layout ONCE, outside the loop
topo = (
    mix_1 = Mixer{:AveragedPressure}(:mixer_1, :feed, :recycle, :mix_out),
    compressor_1 = Compressor{:PressureRatio}(:compressor_1, :mix_out, :c1_out, :electricity_1),
    splitter_1 = TwoSplitter{:Ratio}(:splitter_1, :c1_out, :recycle, :product)
)

function copy_stream!(from, to)
    to .= from
end


function _run_component(
    type::TwoMixer{:AveragedPressure},
    streams, components, system
)
    stream_in_1 = getproperty(streams, type.stream_in_1)
    stream_in_2 = getproperty(streams, type.stream_in_2)
    stream_out = getproperty(streams, type.stream_out)
    component = getproperty(components, type.id)

    copy_stream!(stream_in, stream_out)

    stream_out.pressure = 0.5*(stream_in_1.pressure + stream_in_2.pressure)
    stream_out.mass_flow = stream_in_1.mass_flow + stream_in_2.mass_flow
    stream_out.mass_fractions = 0.5*(stream_in_1.mass_fractions + stream_in_2.mass_fractions)
    update_stream_properties!(stream_out)
end


function _run_component(
    type::Compressor{:PressureIncrease},
    streams, components, system
)
    stream_in = getproperty(streams, type.stream_in)
    stream_out = getproperty(streams, type.stream_out)
    energy_in = getproperty(streams, type.energy_in)
    component = getproperty(components, type.id)

    copy_stream!(stream_in, stream_out)
    
    # Notice how we use the names directly!
    stream_out.pressure = stream_in.pressure * component.pressure_ratio
    
    # ... temperature and energy math ...
end


function _run_component(
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

    stream_out_1.mass_fractions = stream_in.mass_fractions
    stream_out_2.mass_fractions = stream_in.mass_fractions

    update_stream_properties!(stream_out_1)
    update_stream_properties!(stream_out_2)
end

function _run_component(
    type::Compressortype{:PressureRatio},
    streams, components, system
)
    stream_in = getproperty(streams, type.stream_in)
    stream_out = getproperty(streams, type.stream_out)
    energy_in = getproperty(streams, type.energy_in)
    component = getproperty(components, type.id)

    copy_stream!(stream_in, stream_out)
    
    # Notice how we use the names directly!
    stream_out.pressure = stream_in.pressure * component.pressure_ratio
    
    # ... temperature and energy math ...
end

decision_variables = ComponentVector(
    c1_pressure_ratio = 3.0,
    c2_pressure_ratio = 4.0,
    hex1_area = 20.0,
    reactor_volume = 5.0
)

function evaluate_plant_cost(decision_vars, components, streams, system, topo)
    
    # 1. Mutate the components with optimizer guesses (Your VirtualArray magic)
    components.compressor_1.pressure_ratio = decision_vars.c1_pressure_ratio
    components.compressor_2.pressure_ratio = decision_vars.c2_pressure_ratio

    for stream in streams
        update_stream_properties!(stream)
    end

    #TODO: do some graph theory shenanagins to make topo list the components in the right execution order
    foreach(values(topo)) do unit
        _run_component(unit, streams, components, system)
    end

    #this doesn't have to be ordered
    foreach(values(topo)) do unit
        calculate_economics!(unit, streams, components, system)
    end
    
    return system.total_plant_cost
end 