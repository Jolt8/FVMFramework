
struct FaceVectorView{V <: AbstractVector} <: AbstractVector{SubArray}
    data::V
    n_faces::Int
    start_idx::Int
end

@inline Base.size(A::FaceVectorView) = (div(length(A.data) - A.start_idx + 1, A.n_faces),)

@inline function Base.getindex(A::FaceVectorView, i::Int)
    @boundscheck checkbounds(A, i)
    idx = A.start_idx + (i - 1) * A.n_faces
    return @inbounds view(A.data, idx:idx + A.n_faces - 1)
end

@inline function Base.setindex!(A::FaceVectorView, v, i::Int)
    @boundscheck checkbounds(A, i)
    idx = A.start_idx + (i - 1) * A.n_faces
    @inbounds A.data[idx:idx + A.n_faces - 1] .= v
    return v
end

Base.IndexStyle(::Type{<:FaceVectorView}) = IndexLinear()

struct VirtualAxis{Src, Ax}
    ax::Ax
end

struct VirtualFVMArray{D <: Tuple, A <: NamedTuple}
    data::D
    axes::A
end

@inline _resolve(data, idx::Int) = view(data, idx:idx)
@inline _resolve(data, idx::UnitRange{Int}) = view(data, idx)
@inline _resolve(data, idx::ComponentArrays.ComponentIndex) = _apply_axis(view(data, idx.idx), idx.ax)
@inline _resolve(data, ax::ComponentArrays.AbstractAxis) = ComponentArray(data, (ax,))

@inline _apply_axis(data, ax::ComponentArrays.AbstractAxis) = ComponentArray(data, (ax,))
@inline _apply_axis(data, ax::ComponentArrays.ShapedAxis) = reshape(data, size(ax))
@inline _apply_axis(data, ax::ComponentArrays.ViewAxis) = _apply_axis(data, ax.ax)

@inline _resolve(data, ax::Tuple{Int, Int}) = FaceVectorView(data, ax[2], ax[1])
@inline _resolve(data, ax::Tuple{Int}) = view(data, ax[1]:ax[1])

@inline function _virtual_resolve(data::Tuple, vax::VirtualAxis{Src}) where {Src}
    return _resolve(data[Src], vax.ax) 
end

@inline function _virtual_resolve(data::Tuple, ax::NamedTuple)
    return VirtualFVMArray(data, ax)
end

@inline function Base.getproperty(A::VirtualFVMArray, s::Symbol)
    ax = getfield(getfield(A, :axes), s)
    return _virtual_resolve(getfield(A, :data), ax)
end

@inline Base.getindex(A::VirtualFVMArray, s::Symbol) = getproperty(A, s)
@inline Base.getindex(A::VirtualFVMArray, ax::VirtualAxis) = _virtual_resolve(getfield(A, :data), ax)
@inline Base.getindex(A::VirtualFVMArray, ax::NamedTuple) = _virtual_resolve(getfield(A, :data), ax)

@inline Base.setindex!(A::VirtualFVMArray, v, s::Symbol) = (getproperty(A, s) .= v)
@inline Base.setindex!(A::VirtualFVMArray, v, ax::VirtualAxis) = (Base.getindex(A, ax) .= v)
@inline Base.setindex!(A::VirtualFVMArray, v, ax::NamedTuple) = (Base.getindex(A, ax) .= v)

@inline Base.keys(A::VirtualFVMArray) = keys(getfield(A, :axes))

@inline Base.view(A::VirtualFVMArray, s::Symbol) = getproperty(A, s)
@inline Base.dotview(A::VirtualFVMArray, s::Symbol) = getproperty(A, s)

function virtual_merge_axes(component_array_list::Tuple)
    merged_entries = Symbol[]
    merged_values = Any[]
    for (i, component_array) in enumerate(component_array_list)
        axis = getaxes(component_array)[1]
        for name in keys(axis)
            if name in merged_entries
                continue 
            end
            push!(merged_entries, name)
            push!(merged_values, _wrap_virtual(axis[name], i))
        end
    end
    return NamedTuple{Tuple(merged_entries)}(Tuple(merged_values))
end

_wrap_virtual(ax::NamedTuple, src::Int) = NamedTuple{keys(ax)}(map(v -> _wrap_virtual(v, src), values(ax)))
_wrap_virtual(ax, src::Int) = VirtualAxis{src, typeof(ax)}(ax)

@inline function _resolve_unified_val(A::VirtualFVMArray, ::Val{s}) where {s}
    return getproperty(A, s)
end


function _show_axes(io::IO, ax::NamedTuple)
    for name in keys(ax)
        val = ax[name]
        if val isa NamedTuple
            println(io, name, " (group):")
            _show_axes(io, val)
        elseif val isa UnitRange
            println(io, name, " → ", length(val), " cells (indices ", val, ")")
        elseif val isa Tuple{Int,Int}
            n_cells = 0  # can't know without data length, just show the tuple
            println(io, name, " → face-indexed (start=", val[1], ", n_faces=", val[2], ")")
        elseif val isa Tuple{Int}
            println(io, name, " → scalar (index ", val[1], ")")
        else
            println(io, name, " → ", val)
        end
    end
end

function Base.show(io::IO, A::VirtualFVMArray)
    print(io, "VirtualFVMArray(")
    ax = getfield(A, :axes)
    for name in keys(ax)
        println(name)
        val = ax[name]
        println(name)
        println(val)
        println(io, _virtual_resolve(getfield(A, :data), val))
    end
end
