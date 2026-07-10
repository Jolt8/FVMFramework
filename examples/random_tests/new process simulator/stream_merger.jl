
function merge_stream!(streams, stream, stream_id)
    for property_name in propertynames(stream)
        var = getproperty(streams, property_name)
        initial_condition = getproperty(stream, property_name)
        if var isa ComponentVector
            for sub_name in propertynames(var)
                getproperty(var, sub_name)[stream_id] = getproperty(initial_condition, sub_name)
            end
        else
            var[stream_id] = initial_condition
        end
    end
end