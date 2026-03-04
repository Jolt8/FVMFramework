using ComponentArrays
using BenchmarkTools
using PreallocationTools
using Polyester

struct FaceVectorView{T, V <: AbstractVector{T}} <: AbstractVector{AbstractVector{T}}
    v::V
    n_faces::Int
    start_idx::Int
end
#=
function FaceVectorView(v::V, n_faces::Int, start_idx::Int) where {V<:AbstractVector}
    T = eltype(V)
    return FaceVectorView{T, V}(v, n_faces, start_idx)
end
=#

@inline Base.size(A::FaceVectorView) = (length(A.v) ÷ A.n_faces,)
@inline Base.IndexStyle(::Type{<:FaceVectorView}) = Base.IndexLinear()

@inline function Base.getindex(A::FaceVectorView, i::Int)
    idx = A.start_idx + (i - 1) * A.n_faces
    return @inbounds view(A.v, idx:idx+A.n_faces-1)
end

@inline function Base.setindex!(A::Vector, v, ax::FaceVectorView)
    idx = ax.start_idx
    A[idx:idx+ax.n_faces-1] = v
end

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::Tuple{Int, Int}) = FaceVectorView(v, ax[2], ax[1])
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

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

                push!(temp_axes, (property_name => (start_idx, n_faces)))
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

n_cells = 1000

du_proto = (
    mass_fractions = (methanol = zeros(n_cells), water = zeros(n_cells), carbon_monoxide = zeros(n_cells), hydrogen = zeros(n_cells), carbon_dioxide = zeros(n_cells)),
    temp = zeros(n_cells)
)

du_proto_nt = (; du_proto...)

u_proto = (
    mass_fractions = (methanol = zeros(n_cells), water = zeros(n_cells), carbon_monoxide = zeros(n_cells), hydrogen = zeros(n_cells), carbon_dioxide = zeros(n_cells)),
    temp = zeros(n_cells)
)

u_proto_nt = (; u_proto...)

du_cache_nt = (
    rho = zeros(n_cells),
)

u_cache_nt = (
    rho = zeros(n_cells),
)

properties = (rho = ones(n_cells) * 1100.0, cp = ones(n_cells) * 2500.0, k = ones(n_cells) * 0.15)

du_vec = Vector(ComponentArray(; du_proto_nt...))
u_vec = Vector(ComponentArray(; u_proto_nt...))
du_cache_vec = Vector(ComponentArray(; du_cache_nt...))
u_cache_vec = Vector(ComponentArray(; u_cache_nt...))

du_proto_axes = create_axes(du_proto, n_cells)
u_proto_axes = create_axes(u_proto, n_cells)
du_cache_axes = create_axes(du_cache_nt, n_cells)
u_cache_axes = create_axes(u_cache_nt, n_cells)

du_diff_cache_vec = DiffCache(du_cache_vec)
u_diff_cache_vec = DiffCache(u_cache_vec)

using ForwardDiff

function test_ode(du_vec, u_vec, p, t, du_diff_cache_vec, u_diff_cache_vec, properties, du_proto_axes, u_proto_axes, du_cache_axes, u_cache_axes, n_cells)
    du_cache_vec = get_tmp(du_diff_cache_vec, u_vec)
    u_cache_vec = get_tmp(u_diff_cache_vec, u_vec)

    #println(typeof(du_cache_vec))

    #do you know what is fucking crazy?
    #the only reason this worked in the past was because I zeroed out the caches
    #otherwise, it just returns #undef for everything
    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = create_views_inline(du_cache_vec, du_cache_axes)
    u_cache_nt = create_views_inline(u_cache_vec, u_cache_axes)

    du_vec .= 0.0
    du_nt = create_views_inline(du_vec, du_proto_axes)
    u_nt = create_views_inline(u_vec, u_proto_axes)

    du = (; du_nt..., du_cache_nt...)
    #u = (; u_nt..., u_cache_nt..., properties...)
    u = (; properties..., u_nt..., u_cache_nt...)

    u.rho .= (u_cache_nt.rho .+ properties.rho)

    @batch for cell_id in 1:n_cells
        du.mass_fractions.methanol[cell_id] += 1.0
        du.temp[cell_id] += 270.0
        du.rho[cell_id] += 1.0
        u.mass_fractions.methanol[cell_id] += 1.0
        u.temp[cell_id] += 270.0
        u.rho[cell_id] += 1.0
        u.cp[cell_id] += 1.0
        map(keys(du.mass_fractions)) do species_name
            #ok, if we run this without anything inside it, we still get zero allocations
            #ahh, so using keys returned by map(keys) when doing @batch causes a lot of allocations
            #otherwise, it works fine
            #that really sucks
        end
        map(keys(du.mass_fractions)) do species_name
            du.mass_fractions[species_name][cell_id] += 1.0
            #when doing batch, 11000 al
        end
        #map(keys(u.mass_fractions)) do species_name #fuck, not using this map makes this way faster
            du.mass_fractions[:methanol][cell_id] += 1.0
            du.mass_fractions[:water][cell_id] += 1.0
            du.mass_fractions[:carbon_monoxide][cell_id] += 1.0
            du.mass_fractions[:hydrogen][cell_id] += 1.0
            du.mass_fractions[:carbon_dioxide][cell_id] += 1.0
        #end
    end
end

test_ode_closure = (du, u, p, t) -> test_ode(
    du, u, p, t,
    du_diff_cache_vec, u_diff_cache_vec, properties,

    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes,
    1000
)

p_guess = [0.0]
t = 0.0

#VSCodeServer.@profview test_ode_closure(du_vec, u_vec, p_guess, t)
test_ode_closure(du_vec, u_vec, p_guess, t)
@benchmark test_ode_closure($du_vec, $u_vec, $p_guess, $t)
VSCodeServer.@profview [test_ode_closure(du_vec, u_vec, p_guess, t) for _ in 1:100]
