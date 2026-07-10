
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
    ),

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

stream_names = (:feed, :c1_out, :c2_out, :reactor_out, :recycle)
unit_names   = (:comp1, :comp2, :hex1, :r1)
energy_names = (:elec1, :elec2, :heat1, :heat2)

stream_idx = NamedTuple{stream_names}(1:length(stream_names))

unit_idx = NamedTuple{unit_names}(1:length(unit_names))

energy_idx = NamedTuple{energy_names}(1:length(energy_names))

function run_unit!(::Val{:Compressor}, ::Val{:PressureRatio}, streams, components, component_id, stream_in, stream_out, energy_in)
    # The math uses the exact parameters from the specific compressor block
    streams.pressure[stream_out] = streams.pressure[stream_in] * components.pressure_ratio
    # ...
end

decision_variables = ComponentVector(
    c1_pressure_ratio = 3.0,
    c2_pressure_ratio = 4.0,
    hex1_area = 20.0,
    reactor_volume = 5.0
)

function evaluate_plant_cost(streams, components, system, stream_idxs, unit_idxs, energy_idxs)
    # No variable assignment needed! The Duals are already in plant_params.
    
    # Run the units using Value Dispatch
    run_unit!(Val(:Compressor), Val(:PressureRatio), streams, components, unit_idxs.compressor_1, stream_idxs.feed, stream_idxs.c1_out, energy_idxs.elec1)
    run_unit!(Val(:Compressor), Val(:PressureRatio), streams, components, unit_idxs.compressor_2, stream_idxs.c1_out, stream_idxs.c2_out, energy_idxs.elec2)
    run_unit!(Val(:HeatExchanger), Val(:Temporary), streams, components, unit_idxs.hex_1, stream_idxs.c2_out, stream_idxs.recycle, energy_idxs.heat1)
    
    return system.total_plant_cost
end