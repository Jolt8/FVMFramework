using ComponentArrays
using Unitful, UnitfulAssets
using FVMFramework

total_streams_required = 20

species_symbols = (:carbon_dioxide, :water)

streams = ComponentVector(
    mass_flow = [0.0u"kg/s" for _ in 1:total_streams_required],
    mass_flows = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"kg/s" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)])...)
    ),
    molar_flow = [0.0u"mol/s" for _ in 1:total_streams_required],
    molar_flows = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"mol/s" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)])...)
    ),
    volumetric_flow = [0.0u"m^3/s" for _ in 1:total_streams_required],
    volumetric_flows = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"m^3/s" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)])...)
    ),
    temp = [0.0u"K" for _ in 1:total_streams_required],
    pressure = [0.0u"bar" for _ in 1:total_streams_required],
    mass_fractions = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0 for _ in 1:total_streams_required] for _ in eachindex(species_symbols)]) ...)
    ),
    molar_fractions = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0 for _ in 1:total_streams_required] for _ in eachindex(species_symbols)]) ...)
    ),
    molar_concentrations = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"mol/m^3" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)]) ...)
    ),
    #=species_specific_heat_ratios = (
        carbon_dioxide = 1.29
    ),=# #this one is likely not necessary because we're averaging it out anyways
    specific_heat_ratio = [0.0 for _ in 1:total_streams_required],
    molecular_weights = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"g/mol" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)]) ...)
    ),
    mw_avg = [0.0u"g/mol" for _ in 1:total_streams_required],
    rho = [0.0u"kg/m^3" for _ in 1:total_streams_required],
    species_cps = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"J/(kg*K)" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)]) ...)
    ),
    cp_avg = [0.0u"J/(kg*K)" for _ in 1:total_streams_required],
    species_cvs = (
        ComponentVector(; NamedTuple{species_symbols}([[0.0u"J/(kg*K)" for _ in 1:total_streams_required] for _ in eachindex(species_symbols)]) ...)
    ),
    cv_avg = [0.0u"J/(kg*K)" for _ in 1:total_streams_required],
    electrical_wattage = [0.0u"W" for _ in 1:total_streams_required],
    heat_wattage = [0.0u"W" for _ in 1:total_streams_required],

    #constant properties
    R_gas = [8.314u"J/(mol*K)" for _ in 1:total_streams_required],
)

total_components_required = 20

components = ComponentVector(
    #commmon stuff
    component_cost = [0.0u"USD" for _ in 1:total_components_required],
    

    #compressor stuff
    rotation_speed = [], #pattern would continue here
    pressure_increase,
    outlet_pressure,
    pressure_ratio,
    adiabatic_efficiency,
    polytropic_efficiency,
    power_required,
    outlet_temperature,
    temperature_change,
    adiabatic_coefficient,
    polytropic_coefficient,
    adiabatic_head,
    polytropic_head,
    
    #heater/cooler stuff
    pressure_drop,
    efficiency,
    outlet_temperature,
    temperature_change,
    outlet_vapor_fraction, 
    heat_wattage, #heaters are negative, coolers are positive

    #conversion reactor
    reactions = (),
    outlet_temperature, 
    pressure_drop,

    #compound_separator
    stream_1_species_separation_fractions = (
        carbon_dioxide = 0.5,
        water = 0.25
    ),
    stream_2_species_separation_fractions = (
        carbon_dioxide = 0.5, 
        water = 0.75
    ),

    #mixers
    #hmm, guess there's nothing heater

    #splitters
    split_ratios = (
        stream_1 = 0.25,
        stream_2 = 0.5, 
        stream_3 = 0.25
    )

    #heat exchangers
    flow_1_pressure_drop,
    flow_2_pressure_drop,
    flow_1_outlet_temp,
    flow_2_outlet_temp,
    heat_transfer_coefficient = 1000u"W/(m^2*K)",
    heat_exchange_area = 20u"m^2",
    heat_exchanged = 1.0u"W",
    min_temp_difference = 1.0u"K",
    heat_loss = 1.0u"W",
    heat_transfer_efficiency = 99.9u"%",
    flow_1_outlet_vapor_fraction = 0.5,
    flow_2_outlet_vapor_fraction = 0.5,

    #gas liquid separator
    #nothing here either

    #pipe segment
    outlet_temperature,
    outlet_pressure,
    inside_diameter,
    outside_diameter,
    relative_roughness,
    length, 
    elevation,
    heat_loss,
    heat_transfer_coefficient,
    pipe_outside_area, #this one will likely always be calculated externally
    
    #valves
    pressure_drop,
    outlet_pressure,
)

#system will likely become a struct later
system = ComponentVector(
    total_plant_cost = 0.0u"USD",
    upfront_cost = 0.0u"USD",
    interest_rate = 0.0u"USD"
)

uncompressed_stream_id = 1
uncompressed_stream = ComponentVector(
    mass_flow = 1.0u"kg/s",
    temp = 25.0u"°C",
    pressure = 1.0u"bar",
    mass_fractions = (
        carbon_dioxide = 1.0,
        water = 0.0
    ),
    molecular_weights = (
        carbon_dioxide = 44.01u"g/mol",
        water = 18.02u"g/mol",
    ),
    species_cps = (
        carbon_dioxide = 846.0u"J/(kg*K)",
        water = 4186u"J/(kg*K)",
    ),
    species_cvs = (
        carbon_dioxide = 657.0u"J/(kg*K)",
        water = 3377u"J/(kg*K)",
    ),
)

merge_stream!(streams, uncompressed_stream, uncompressed_stream_id)

Revise.includet(joinpath(@__DIR__, "process_sim_physics_helper_functions.jl"))
Revise.includet(joinpath(@__DIR__, "stream_copier.jl"))

function update_stream_properties!(streams, stream_id)
    mw_avg!(streams, stream_id)
    rho_ideal!(streams, stream_id)
    molar_concentrations!(streams, stream_id)
    molar_fractions!(streams, stream_id)
    mass_flows!(streams, stream_id)
    molar_flow!(streams, stream_id)
    molar_flows!(streams, stream_id)
    volumetric_flow!(streams, system, stream_id)
    volumetric_flows!(streams, system, stream_id)
    cp_avg!(streams, stream_id)
    cv_avg!(streams, stream_id)
    specific_heat_ratio!(streams, stream_id)
end

function update_electrical_energy_stream_properties!(streams, stream_id)
end
function update_heat_energy_stream_properties!(streams, stream_id)
end

function compressor!(streams, stream_in_id, stream_out_id, energy_stream_in_id, output_pressure, efficiency)
    copy_stream!(streams, stream_in_id, stream_out_id)

    streams.pressure[stream_out_id] = output_pressure

    streams.temp[stream_out_id] = (streams.temp[stream_in_id] |> u"K") * (streams.pressure[stream_out_id] / streams.pressure[stream_in_id]) ^ ((streams.specific_heat_ratio[stream_in_id] - 1) / streams.specific_heat_ratio[stream_in_id])

    @show streams.molar_flows[stream_in_id]
    @show streams.R_gas[stream_in_id]
    @show streams.temp[stream_out_id]
    @show streams.temp[stream_in_id]
    @show streams.specific_heat_ratio[stream_in_id]

    ((streams.temp[stream_out_id] |> u"K") - (streams.temp[stream_in_id] |> u"K"))

    ideal_wattage_required = streams.molar_flow[stream_in_id] * streams.R_gas[stream_in_id] * ((streams.temp[stream_out_id] |> u"K") - (streams.temp[stream_in_id] |> u"K")) / (streams.specific_heat_ratio[stream_in_id] - 1.0)

    streams.electrical_wattage[energy_stream_in_id] -= ideal_wattage_required / efficiency
end

function cooler!(streams, stream_in_id, stream_out_id, energy_stream_out_id, outlet_temperature)
    copy_stream!(streams, stream_in_id, stream_out_id)

    streams.temp[stream_out_id] = outlet_temperature

    delta_temp = (streams.temp[stream_in_id] |> u"K") - (streams.temp[stream_out_id] |> u"K")

    streams.heat_wattage[energy_stream_out_id] += streams.cp_avg[stream_in_id] * streams.mass_flow[stream_in_id] * delta_temp
end

function multistage_compressor_with_intercoolers!(
    streams,
    ordered_compressor_input_stream_ids, #first is input stream
    ordered_compressor_output_stream_ids, #last it output stream
    input_electrical_energy_stream_ids,
    n_compressors,
    output_heat_energy_stream_ids,
    output_pressure,
    efficiencies,
    desired_temperatures_after_cooling
)
    for stream_id in ordered_compressor_input_stream_ids
        update_stream_properties!(streams, stream_id)
    end
    for stream_id in ordered_compressor_output_stream_ids
        update_stream_properties!(streams, stream_id)
    end
    for stream_id in input_electrical_energy_stream_ids
        update_electrical_energy_stream_properties!(streams, stream_id)
    end
    for stream_id in output_heat_energy_stream_ids
        update_heat_energy_stream_properties!(streams, stream_id)
    end

    for i in 1:n_compressors
        compressor!(streams, ordered_compressor_input_stream_ids[i], ordered_compressor_output_stream_ids[i], input_electrical_energy_stream_ids[i], output_pressure[i], efficiencies[i])
        if i != n_compressors
            cooler!(streams, ordered_compressor_output_stream_ids[i], ordered_compressor_input_stream_ids[i+1], output_heat_energy_stream_ids[i], desired_temperatures_after_cooling[i])
        end
    end
end

multistage_compressor_with_intercoolers!(
    streams,
    [1, 3], [2, 4],
    [5, 6], [7, 8],
    2,
    [3.0u"bar", 10.0u"bar"],
    [1.0, 1.0],
    [25.0u"°C", 25.0u"°C"]
)
