#this is bad, using pairs on a component vectors doesn't actually return symbols and instead returns Ints
#however, we should probably find a way around this because using pairs on a component vector is faster than using pairs on a NamedTuple
temp_controller = (
    proportional_gain = 1.0,
    integral_time = 1.0,
    derivative_time = 1.0,
    desired_value_comp_vector = ComponentVector(temp = ustrip(270.0u"°C" |> u"K"), pressure = [2.0, 2.0]),
    initial_volumetric_input = 1e7,
    min_volumetric_input = 0.0,
    max_volumetric_input = 1e8
)

test = (pressure = [2.0, 2.0], temp = 1.0)

pairs(temp_controller.desired_value_comp_vector)

for (field, _) in pairs(temp_controller.desired_value_comp_vector)
    println(field)
    println(temp_controller.desired_value_comp_vector[field])
    println(test[field])
end