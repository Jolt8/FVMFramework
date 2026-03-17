using ComponentArrays

@inline _resolve(data::AbstractVector, ax::UnitRange{Int}) = view(data, ax)

@inline _resolve(data::AbstractVector, ax::NamedTuple) = FVMArray(data, ax)

@inline _resolve(data::AbstractVector, ax::Tuple{Int,Int}) = FaceVectorView(data, ax[2], ax[1])

@inline _resolve(data::AbstractVector, ax::Tuple{Int}) = view(data, ax[1]:ax[1])

@inline function _resolve_val(A::ComponentVector, ::Val{s}) where {s}
    return _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
end

@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{ComponentVector, N}) where {N}
    #println(fieldnames(groups[1].parameters[1]))
    #println(fieldnames(groups[1].parameters[3]))
    println(fieldnames(groups[1].parameters[3]))
    first_ax = fieldnames(groups[1].parameters[3])
    exprs = []
    for n in first_ax
        view_exprs = [:(_resolve_val(groups[$i], Val{$(QuoteNode(n))}())) for i in 1:N]
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

test = ComponentArray(
    mass_fractions = (
        methanol = zeros(100),
        water = zeros(100),
    )
)

getfield(test, :data)
getfield(test, :axes)

foreach_field_at!(1, test.mass_fractions) do species_name, mass_fractions

end

getaxes(test.parent)