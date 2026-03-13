function append_axes!(temp_axes, axes_length, data, property_name, n_cells)
    if !(property_name in keys(temp_axes))
        if data[property_name] isa NamedTuple
            push!(temp_axes, (property_name => Dict{Symbol,Any}()))

            for sub_name in keys(data[property_name])
                axes_length = append_axes!(temp_axes[property_name], axes_length, data[property_name], sub_name, n_cells)
            end
        elseif data[property_name] isa Vector
            if data[property_name][1] isa Vector #stuff like u.mass_face[cell_id][face_idx] that's tracked per face 
                n_faces = length(data[property_name][1])
                start_idx = axes_length + 1

                push!(temp_axes, (property_name => (start_idx, n_faces, n_cells)))
                axes_length += n_cells * n_faces
            elseif data[property_name][1] isa Number #this is for basically everything else including n_cell properties
                n = length(data[property_name])

                push!(temp_axes, (property_name => ((axes_length+1):(axes_length+n))))
                axes_length += n
            end
        elseif data[property_name] isa Number
            push!(temp_axes, property_name => (axes_length+1):(axes_length+1))
            axes_length += 1
            #=else
                @warn "$(property_name) was not handled"
                push!(temp_axes, property_name => (axes_length+1):(axes_length+1))
                axes_length += 1
            =#
        end
    end
    return axes_length
end

function nested_dict_to_named_tuple(d::Dict)
    pairs = map(collect(d)) do (key, value)
        if value isa Dict
            Symbol(key) => nested_dict_to_named_tuple(value)
        else
            Symbol(key) => value
        end
    end
    return NamedTuple(pairs)
end

function create_axes(data, n_cells)
    temp_axes = Dict{Any, Any}() #Pair{Symbol, UnitRange{Int64}}
    axes_length = 0

    for property_name in keys(data)
        axes_length = append_axes!(temp_axes, axes_length, data, property_name, n_cells)
    end

    return nested_dict_to_named_tuple(temp_axes)
end