struct FaceVectorView{T, V <: AbstractVector{T}} <: AbstractVector{AbstractVector{T}}
    v::V
    n_faces::Int
    start_idx::Int
    len::Int
end
#=
function FaceVectorView(v::V, n_faces::Int, start_idx::Int) where {V<:AbstractVector}
    T = eltype(V)
    return FaceVectorView{T, V}(v, n_faces, start_idx)
end
=#

@inline Base.size(A::FaceVectorView) = (A.len,)
@inline Base.IndexStyle(::Type{<:FaceVectorView}) = Base.IndexLinear()

@inline function Base.getindex(A::FaceVectorView, i::Int)
    #@boundscheck checkbounds(A, i)
    idx = A.start_idx + (i - 1) * A.n_faces
    return @inbounds view(A.v, idx:idx+A.n_faces-1)
end

@inline function Base.setindex!(A::Vector, v, ax::FaceVectorView)
    idx = ax.start_idx
    A[idx:idx+ax.n_faces-1] = v
end

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::Tuple{Int, Int, Int}) = FaceVectorView(v, ax[2], ax[1], ax[3])
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)