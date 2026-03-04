using Unitful
using BenchmarkTools
using Polyester
using OrdinaryDiffEq
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
            #=
            elseif data[property_name] isa Vector #stuff like u.mass_face[cell_id][face_idx] that's tracked per face 
                push!(temp_axes, (property_name => Dict{Int, Any}()))
                for i in eachindex(data[property_name])
                    push!(temp_axes[property_name], i => (axes_length + 1):(axes_length + n_cells))
                    axes_length += n_cells
                end
            =#
        elseif data[property_name] isa Vector
            if data[property_name][1] isa Vector #stuff like u.mass_face[cell_id][face_idx] that's tracked per face 
                n_faces = length(data[property_name][1])
                start_idx = axes_length + 1

                push!(temp_axes, (property_name => (start_idx, n_faces)))
                axes_length += n_cells * n_faces
            elseif data[property_name][1] isa Number
                push!(temp_axes, property_name => (axes_length+1):(axes_length+n_cells))
                axes_length += n_cells
            end
        elseif data[property_name] isa Number
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
    temp_axes = Dict{Any, Any}() #Pair{Symbol, UnitRange{Int64}}
    axes_length = 0

    for property_name in keys(data)
        axes_length = append_axes!(temp_axes, axes_length, data, property_name, n_cells)
    end

    return nested_dict_to_named_tuple(temp_axes), axes_length
end

struct FaceVectorView{V} <: AbstractVector{AbstractVector{Float64}}
    v::V
    n_faces::Int
    start_idx::Int
end

@inline Base.size(A::FaceVectorView) = (length(A.v) ÷ A.n_faces,)
@inline function Base.getindex(A::FaceVectorView, i::Int)
    idx = A.start_idx + (i - 1) * A.n_faces
    return view(A.v, idx:idx+A.n_faces-1)
end

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::Tuple{Int, Int}) = FaceVectorView(v, ax[2], ax[1])
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

test_axes, axes_length = create_axes(inlet_properties, n_cells)
test_vector = rand(axes_length)

println("==== FaceVectorView Creation ====")
@btime FaceVectorView($test_vector, 4, 1)

println("typeof(test_axes) = ", typeof(test_axes))

println("==== _create_views_inline Initial Call ====")
@time test = create_views_inline(test_vector, test_axes)

function update_cells(properties, cell_ids, face_idxs)
    for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            properties.mass_face[cell_id][face_idx] += 1.0
        end
    end
end

function update_cells_batch(properties, cell_ids, face_idxs)
    @batch for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            properties.mass_face[cell_id][face_idx] += 1.0
        end
    end
end

function run_benchmarks(test_vector, test_axes, cell_ids, face_idxs)
    properties = create_views_inline(test_vector, test_axes)

    println("==== benchmark for update_cells! ====")
    @btime update_cells($properties, $cell_ids, $face_idxs)

    println("==== benchmark for update_cells_batch! ====")
    @btime update_cells_batch($properties, $cell_ids, $face_idxs)
end

cell_ids = collect(1:1000)
face_idxs = [[1, 2, 3, 4] for _ in eachindex(cell_ids)]

run_benchmarks(test_vector, test_axes, cell_ids, face_idxs)

@btime create_views_inline($test_vector, $test_axes)
