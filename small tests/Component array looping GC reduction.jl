using ComponentArrays
using ProfileView

ca = ComponentArray(a = 1.0:10.0, b = 11.0:20.0)

#Base.getindex(x::ComponentVector{Float64, Vector{Float64}, Tuple{Axis}}, f::Symbol) = view(x::ComponentVector{Float64, Vector{Float64}, Tuple{Axis}}, f::Symbol)
#Base.setindex!(x::ComponentVector{Float64, Vector{Float64}, Tuple{Axis}}, val::Float64, f::Symbol) = view(x::ComponentVector{Float64, Vector{Float64}, Tuple{Axis}}, f::Symbol) = val

function test_loop(ca)
    for property_name in propertynames(ca)
        #println(typeof(ca))
        #println(typeof(ca[property_name]))
        ca[property_name][1] = 2.0 #this causes the most GC
        ca[property_name] .= 2.0 #this causes the most GC
    end
end

VSCodeServer.@profview test_loop(ca)