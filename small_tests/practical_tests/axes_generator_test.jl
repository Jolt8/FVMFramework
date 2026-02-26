using Unitful
using BenchmarkTools
using Polyester
using OrdinaryDiffEq

cells = collect(1:1000)
n_cells = length(cells)

inlet_properties = (temp=ustrip(21u"°C" |> u"K"), mass_fractions=(methanol=0.5, water=0.5), mass_face=[1.0, 1.0, 1.0, 1.0])
inlet_cellset = collect(300:400)

(length(cells)+1):(length(cells)+n_cells)

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

    return nested_dict_to_named_tuple(temp_axes), axes_length
end

test_axes, axes_length = create_axes(inlet_properties, n_cells)
test_vector = rand(axes_length)

#i'm suprised the above works, I can't believe that doing Symbol(Int64) actually gets the right index of the Dict

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::NTuple{N, UnitRange{Int64}}) where N = ntuple(i -> view(v, ax[i]), Val(N)) # 14.5 μs, 1 allocation, 39.25 KiB
#however, the tuple version makes the first time create_views_inline is called take absolutely forever (> 10 seconds)
#actually, I think that was just because it was printing everything to the terminal when it was not called inside a function 
#@inline create_views_inline(v, ax::Vector{UnitRange}) = [view(v, a) for a in ax] # 43.2 μs, 2011 allocations, 102.24 KiB
#@inline create_views_inline(v, ax::Vector{UnitRange}) = map(a -> view(v, a), ax) # 47.1 μs, 2020 allocations, 102.60 KiB
#@inline create_views_inline(v, ax::Vector{UnitRange{Int64}}) = (view(v, a) for a in ax) #this seemed like it worked at one point, but now it doesn't anymore
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

function test!(test_vector, test_axes)
    zoop = create_views_inline(test_vector, test_axes)
    return nothing
end

@benchmark test!(test_vector, test_axes)

@time test_properties = create_views_inline(test_vector, test_axes)

function update_mass_face(properties, cell_id, face_idx)
    properties.mass_face[cell_id][face_idx] += 1.0
end

function update_cells(vector, axes, cell_ids, face_idxs)
    properties = create_views_inline(vector, axes)

    @batch for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            update_mass_face(properties, cell_id, face_idx)
        end
    end
end

cell_ids = collect(1:1000)
face_idxs = [[1, 2, 3, 4] for _ in eachindex(cell_ids)]

@benchmark update_cells(test_vector, test_axes, cell_ids, face_idxs) 
#compared to the actual the solver, this one has 0 allocations and 0 bytes of data used
#ok, I've figured something out, it seems like simply using @batch causes 1 allocation and 39.25 KiB

function test_f!(du_vec, u_vec, p, t, axes, cell_ids, face_idxs)
    du = create_views_inline(du_vec, axes)
    u = create_views_inline(u_vec, axes)

    @batch for cell_id in cell_ids
        du.temp[cell_id] += 1.0
        du.temp[cell_id] += 1.0 #oh no, just adding this increases the number of allocations from 3388 to 3548.
        #du.mass_fractions.methanol[cell_id] += 1.0 #just adding this (not even in addition to the one before) increases allocations from 3388 to 3842
    end

    @batch for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            #update_mass_face(du, cell_id, face_idx)
        end
    end
end

cell_ids = collect(1:100)
face_idxs = [[1, 2, 3, 4] for _ in eachindex(cell_ids)]

inlet_properties = (temp = ustrip(21.13u"°C" |> u"K"), mass_fractions = (methanol = 0.5, water = 0.5), mass_face = [0.0, 0.0, 0.0, 0.0])

cells = collect(1:100)
axes, vector_size = create_axes(inlet_properties, length(cells))

u0_vec = ones(vector_size)

test_f_closure! = (du, u, p, t) -> test_f!(du, u, p, t, axes, cell_ids, face_idxs)

prob = ODEProblem(test_f_closure!, u0_vec, (0.0, 10000.0), [0.0])
@time sol = solve(prob, Tsit5())

#VSCodeServer.@profview sol = solve(prob, Tsit5())

using SparseConnectivityTracer
import ADTypes

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> test_f_closure!(du, u, 0.0, 0.0), copy(u0_vec), u0_vec, detector
)

ode_func = ODEFunction(test_f_closure!, jac_prototype=float.(jac_sparsity))
implicit_prob = ODEProblem(ode_func, u0_vec, (0.0, 1000.0), 0.0)

@time sol = solve(implicit_prob, FBDF(), saveat=1000/1000)
#VSCodeServer.@profview sol = solve(implicit_prob, FBDF())