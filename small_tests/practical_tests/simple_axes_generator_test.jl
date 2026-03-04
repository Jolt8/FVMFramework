using Unitful
using BenchmarkTools
using Polyester
using ComponentArrays #this is just used for flatting NamedTuples at this point 

cells = collect(1:1000)
n_cells = length(cells)

inlet_properties = (temp=ustrip(21u"°C" |> u"K") .* ones(n_cells), mass_fractions=(methanol=0.5 .* ones(n_cells), water=0.5 .* ones(n_cells)), mass_face=[[1.0, 1.0, 1.0, 1.0] for _ in 1:n_cells])

typeof(inlet_properties.mass_face)
flat_properties = Vector(ComponentVector(inlet_properties))

function append_axes!(temp_axes, axes_length, data, property_name, n_cells)
    if !(property_name in keys(temp_axes))
        if data[property_name] isa NamedTuple
            push!(temp_axes, (property_name => Dict{Symbol,Any}()))

            for sub_name in keys(data[property_name])
                append_axes!(temp_axes[property_name], axes_length, data[property_name], sub_name, n_cells)
                axes_length += n_cells
            end
        elseif data[property_name] isa Vector
            if data[property_name][1] isa Vector #stuff like u.mass_face[cell_id][face_idx] that's tracked per face 
                n_faces = length(data[property_name][1])
                #=
                temporary_unit_range_vec = []

                for cell_id in 1:n_cells
                    push!(temporary_unit_range_vec, (axes_length+1):(axes_length+n_faces))
                    axes_length += n_faces
                end

                push!(temp_axes, (property_name => Tuple(temporary_unit_range_vec)))
                =#
                matrix_size = n_cells * n_faces
                push!(temp_axes, (property_name => (axes_length+1):(axes_length+matrix_size)))
                axes_length += matrix_size
                #


                #

                #=
                push!(temp_axes, (property_name => UnitRange[]))

                for cell_id in 1:n_cells
                    push!(temp_axes[property_name], (axes_length + 1):(axes_length + n_faces))
                    axes_length += n_faces
                end
                =#
            elseif data[property_name][1] isa Number
                push!(temp_axes, property_name => (axes_length+1):(axes_length+n_cells))
                axes_length += n_cells
            end
        else
            data[property_name] isa Number
            push!(temp_axes, property_name => (axes_length+1):(axes_length+1))
            axes_length += 1
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
    temp_axes = Dict{Any,Any}() #Pair{Symbol, UnitRange{Int64}}
    axes_length = 0

    for property_name in keys(data)
        axes_length = append_axes!(temp_axes, axes_length, data, property_name, n_cells)
    end

    return nested_dict_to_named_tuple(temp_axes), axes_length
end

test_axes, axes_length = create_axes(inlet_properties, n_cells)
test_vector = rand(axes_length)

#i'm suprised the above works, I can't believe that doing Symbol(Int64) actually gets the right index of the Dict

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
#@inline create_views_inline(v, ax::NTuple{N, UnitRange{Int64}}) where N = ntuple(i -> view(v, ax[i]), Val(N)) # 14.5 μs, 1 allocation, 39.25 KiB
#however, the tuple version makes the first time create_views_inline is called take absolutely forever (> 10 seconds)
#actually, I think that was just because it was printing everything to the terminal when it was not called inside a function 
#@inline create_views_inline(v, ax::Vector{UnitRange}) = [view(v, a) for a in ax] # 43.2 μs, 2011 allocations, 102.24 KiB
#@inline create_views_inline(v, ax::Vector{UnitRange}) = map(a -> view(v, a), ax) # 47.1 μs, 2020 allocations, 102.60 KiB
#@inline create_views_inline(v, ax::Vector{UnitRange{Int64}}) = (view(v, a) for a in ax) #this seemed like it worked at one point, but now it doesn't anymore
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

function update_mass_face(properties, cell_id, face_idx)
    n_cells = length(properties.temp)
    idx = cell_id + ((face_idx - 1) * n_cells)
    properties.mass_face[idx] += 1.0
end

function update_cells(vector, axes, cell_ids, face_idxs)
    properties = create_views_inline(vector, axes)

    for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            update_mass_face(properties, cell_id, face_idx)
        end
    end

end

function update_cells_batch(vector, axes, cell_ids, face_idxs)
    properties = create_views_inline(vector, axes)

    @batch for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            update_mass_face(properties, cell_id, face_idx)
        end
    end
end

cell_ids = collect(1:1000)
face_idxs = [[1, 2, 3, 4] for _ in eachindex(cell_ids)]

length(create_views_inline(test_vector, test_axes).mass_face)

@benchmark update_cells(test_vector, test_axes, cell_ids, face_idxs)
@benchmark update_cells_batch(test_vector, test_axes, cell_ids, face_idxs)

function raw_matrix_test(matrix, cell_ids, face_idxs)
    @batch for cell_id in cell_ids
        for face_idx in face_idxs
            matrix[cell_id, face_idx] += 1.0
        end
    end
end

matrix = zeros(1000, 4)
cell_ids = collect(1:1000)
face_idxs = [1, 2, 3, 4]

@benchmark raw_matrix_test(matrix, cell_ids, face_idxs)