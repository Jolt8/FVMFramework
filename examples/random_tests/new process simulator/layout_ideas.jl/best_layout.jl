stream_names = (:feed, :c1_out, :c2_out, :reactor_out, :recycle)
component_names   = (:comp1, :comp2, :hex1, :r1)
energy_names = (:elec1, :elec2, :heat1, :heat2)

total_streams_required = length(stream_names)

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

total_components_required = length(component_names)

components = ComponentVector(
    pressure_ratio = ComponentVector(; NamedTuple{component_names}([0.0 for _ in 1:total_streams_required])...),
    area = ComponentVector(; NamedTuple{component_names}([0.0u"m^3" for _ in 1:total_streams_required])...)
    # ...
)

#system will likely become a struct later
system = ComponentVector(
    total_plant_cost = 0.0u"USD",
    upfront_cost = 0.0u"USD",
    interest_rate = 0.0u"USD"
)

struct CompressorMeta{CalcMode}
    id::Int
    stream_in::Int
    stream_out::Int
    energy_in::Int
end

stream_idxs = NamedTuple{stream_names}(1:length(stream_names))

component_idxs = NamedTuple{component_names}(1:length(component_names))

energy_idxs = NamedTuple{energy_names}(1:length(energy_names))

# Define your plant layout ONCE, outside the loop
topo = (
    comp1 = CompressorMeta{:PressureRatio}(component_idxs.comp1, stream_idxs.feed, stream_idxs.c1_out, energy_idxs.elec1),
    comp2 = CompressorMeta{:PressureRatio}(component_idxs.comp2, stream_idxs.c1_out, stream_idxs.c2_out, energy_idxs.elec2)
)

components.pressure_ratio.compressor_1 = 3.0
components.pressure_ratio.compressor_2 = 3.0
components.area.heat_exchanger_1 = 20.0
components.area.reactor_1 = 5.0

function _run_component(
    meta::CompressorMeta{:PressureRatio}, 
    streams, components, system, 
    component_id, stream_in, stream_out
)
    copy_stream!(streams, stream_in, stream_out)
    
    # Notice how we use the names directly!
    streams.pressure[stream_out] = streams.pressure[stream_in] * components.pressure_ratio[component_id]
    
    # ... temperature and energy math ...
end


function run_component!(
    meta::CompressorMeta{:PressureRatio}, 
    streams, components, system
)
    # Extract names for clean reading
    stream_in = meta.stream_in
    stream_out = meta.stream_out
    component_id = meta.id
    
    _run_component(meta, streams, components, system, component_id, stream_in, stream_out)
end

decision_variables = ComponentVector(
    c1_pressure_ratio = 3.0,
    c2_pressure_ratio = 4.0,
    hex1_area = 20.0,
    reactor_volume = 5.0
)

function evaluate_plant_cost(decision_vars, components, streams, system, topo)
    
    # 1. Mutate the components with optimizer guesses (Your VirtualArray magic)
    components.pressure_ratio.comp1 = decision_vars.c1_pressure_ratio
    components.pressure_ratio.comp2 = decision_vars.c2_pressure_ratio

    # 2. Run the flowsheet 
    # (No magic numbers! The topology dictates the connections)
    run_component!(topo.compressor_1, streams, components, system)
    run_component!(topo.compressor_2, streams, components, system)
    
    # 3. Cost calculation
    calculate_economics!(topo.compressor_1, streams, components, system)
    calculate_economics!(topo.compressor_2, streams, components, system)
    
    return system.total_plant_cost
end