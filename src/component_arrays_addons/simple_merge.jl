function merge_properties(a::ComponentArray, b::ComponentArray)
    return ComponentVector(; NamedTuple(a)..., NamedTuple(b)...)
end