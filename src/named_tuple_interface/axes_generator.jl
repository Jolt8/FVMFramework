
function append_axes!(temp_axes, axes_length, data, property_name, n_cells)
    if !(property_name in keys(temp_axes))
        if data[property_name] isa NamedTuple
            push!(temp_axes, (property_name => Dict{Symbol, Any}()))

            for sub_name in keys(data[property_name])
                append_axes!(temp_axes[property_name], axes_length, data[property_name], sub_name, n_cells)
                axes_length += n_cells
            end
            #=
            elseif data[property_name] isa Vector #stuff like u.mass_face[cell_id][face_idx] that's tracked per face 
                push!(temp_axes, (property_name => Dict{Int, Any}()))
                for i in eachindex(data[property_name])
                    push!(temp_axes[property_name], i => (axes_length + 1):(axes_length + n_cells))
                    axes_length += n_cells
                end
            =#
        elseif data[property_name] isa Vector #stuff like u.mass_face[cell_id][face_idx] that's tracked per face 
            n_faces = length(data[property_name])

            #Vector version
            #=
            push!(temp_axes, (property_name => UnitRange[]))

            for cell_id in 1:n_cells
                push!(temp_axes[property_name], (axes_length + 1):(axes_length + n_faces))
                axes_length += n_faces
            end
            =#


            #Tuple Version
            #
            temporary_unit_range_vec = []

            for cell_id in 1:n_cells
                push!(temporary_unit_range_vec, (axes_length+1):(axes_length+n_faces))
                axes_length += n_faces
            end

            push!(temp_axes, (property_name => Tuple(temporary_unit_range_vec)))
            #

        elseif data[property_name] isa Number
            push!(temp_axes, property_name => (axes_length+1):(axes_length+n_cells))
            axes_length += n_cells
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