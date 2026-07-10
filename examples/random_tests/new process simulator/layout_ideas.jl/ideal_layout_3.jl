struct Compressor
end

function run_unit!(comp::Compressor{:OutletPressure}, streams, stream_in, stream_out, energy_stream, target_pressure)
    copy_stream!(streams, stream_in, stream_out)
    streams.pressure[stream_out] = target_pressure
    # ... do the math, update temp and energy ...
end

# If the compressor is operating in :PressureRatio mode
function run_unit!(comp::Compressor{:PressureRatio}, streams, stream_in, stream_out, energy_stream, pressure_ratio)
    copy_stream!(streams, stream_in, stream_out)
    streams.pressure[stream_out] = streams.pressure[stream_in] * pressure_ratio
    # ... do the math, update temp and energy ...
end

p = ComponentVector(
    compressor_1 = (pressure_ratio = 3.0, efficiency = 0.8),
    compressor_2 = (pressure_ratio = 4.0, efficiency = 0.85),
    hex_1 = (area = 20.0, U = 1000.0),
    reactor_1 = (volume = 5.0, temperature = 300.0)
)

function run_unit!(::Val{:Compressor}, ::Val{:PressureRatio}, streams, params, stream_in, stream_out, energy_in)
    # The math uses the exact parameters from the specific compressor block
    streams.pressure[stream_out] = streams.pressure[stream_in] * params.pressure_ratio
    # ...
end

decision_variables = ComponentVector(
    c1_pressure_ratio = 3.0,
    c2_pressure_ratio = 4.0,
    hex1_area = 20.0,
    reactor_volume = 5.0
)

function evaluate_plant_cost(p, streams, system)
    # No variable assignment needed! The Duals are already in plant_params.
    
    # Run the units using Value Dispatch
    run_unit!(Val(:Compressor), Val(:PressureRatio), streams, p.compressor_1, 1, 2, 1)
    run_unit!(Val(:Compressor), Val(:PressureRatio), streams, p.compressor_2, 2, 3, 2)
    run_unit!(Val(:HeatExchanger), Val(:Temporary), streams, p.hex_1, 3, 4, 3)
    
    return system.total_plant_cost
end