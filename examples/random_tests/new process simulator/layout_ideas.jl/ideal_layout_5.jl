stream_names = (:feed, :c1_out, :c2_out, :reactor_out, :recycle)
unit_names   = (:comp1, :comp2, :hex1, :r1)
energy_names = (:elec1, :elec2, :heat1, :heat2)

streams = ComponentVector(
    temp = ComponentVector(; NamedTuple{stream_names}(zeros(5))...),
    pressure = ComponentVector(; NamedTuple{stream_names}(zeros(5))...)
    # ...
)

components = ComponentVector(
    pressure_ratio = ComponentVector(; NamedTuple{unit_names}(zeros(4))...),
    area = ComponentVector(; NamedTuple{unit_names}(zeros(4))...)
    # ...
)

struct CompressorMeta{CalcMode}
    id::Int
    stream_in::Int
    stream_out::Int
    energy_in::Int
end

stream_idx = NamedTuple{stream_names}(1:length(stream_names))

unit_idx = NamedTuple{unit_names}(1:length(unit_names))

energy_idx = NamedTuple{energy_names}(1:length(energy_names))

# Define your plant layout ONCE, outside the loop
plant_topology = (
    comp1 = CompressorMeta{:PressureRatio}(unit_idx.comp1, stream_idx.feed, stream_idx.c1_out, energy_idx.elec1),
    comp2 = CompressorMeta{:PressureRatio}(unit_idx.comp2, stream_idx.c1_out, stream_idx.c2_out, energy_idx.elec2)
)

function _run_unit(streams, components, system, meta::CompressorMeta{:PressureRatio}, component_id, stream_in, stream_out)
    copy_stream!(streams, stream_in, stream_out)
    
    # Notice how we use the names directly!
    streams.pressure[stream_out] = streams.pressure[stream_in] * components.pressure_ratio[component_id]
    
    # ... temperature and energy math ...
end


function run_unit!(meta::CompressorMeta{:PressureRatio}, streams, components, system)
    # Extract names for clean reading
    stream_in = meta.stream_in
    stream_out = meta.stream_out
    component_id = meta.id
    
    _run_unit(streams, components, system, meta::CompressorMeta{:PressureRatio}, component_id, stream_in, stream_out)
end

decision_variables = ComponentVector(
    c1_pressure_ratio = 3.0,
    c2_pressure_ratio = 4.0,
    hex1_area = 20.0,
    reactor_volume = 5.0
)

function evaluate_plant_cost(decision_vars, components, streams, system, plant_topology)
    
    # 1. Mutate the components with optimizer guesses (Your VirtualArray magic)
    components.pressure_ratio.comp1 = decision_vars.c1_pressure_ratio
    components.pressure_ratio.comp2 = decision_vars.c2_pressure_ratio

    # 2. Run the flowsheet 
    # (No magic numbers! The topology dictates the connections)
    run_unit!(plant_topology.c1, streams, components, system)
    run_unit!(plant_topology.c2, streams, components, system)
    
    # 3. Cost calculation
    calculate_economics!(plant_topology.c1, streams, components, system)
    calculate_economics!(plant_topology.c2, streams, components, system)
    
    return system.total_plant_cost
end