using Unitful
using PreallocationTools
using ComponentArrays
using benchmarkTools

function test_speed(cached_flat_vec, vec_axes)
    cache = get_tmp(cached_flat_vec, 1.0)

    cache_ca = ComponentArray(cache, vec_axes)
end

vec_ca = ComponentArray(kg = 1.0u"kg", mol = 2.0u"mol")
vec_axes = getaxes(vec)

function unitful_DiffCache(v, N)
    dual_vec = reduce(vcat, [v .* 0.0 for i in 1:chunk_size])

    return DiffCache{Vector{Quantity{Float64}}, Vector{Quantity{Float64}}}(Vector(vec_ca), dual_vec, Any[])
end

unitful_DiffCache(vec_ca, 12)

@btime test_speed($cached_flat_vec, $vec_axes)

cache = get_tmp(cached_flat_vec, 1.0)





using ForwardDiff

ForwardDiff.Dual(1.0u"kg", 1.0)
ForwardDiff.Dual(1.0, 1.0)

function PreallocationTools.get_tmp(dc::PreallocationTools.DiffCache, u::T) where {T <: ForwardDiff.Dual}
    if isbitstype(T) 
        nelem = div(sizeof(T), sizeof(eltype(dc.dual_du))) * length(dc.du)
        if nelem > length(dc.dual_du)
            PreallocationTools.enlargediffcache!(dc, nelem)
        end
        PreallocationTools._restructure(dc.du, reinterpret(T, view(dc.dual_du, 1:nelem)))
    else
        #PreallocationTools._restructure(dc.du, zeros(Float64, size(dc.du)))
        #PreallocationTools._restructure(dc.du, zeros(T, size(dc.du)))
        return [ForwardDiff.Dual(dc.du[i], (0.0 * (unit(dc.du[i]) / u"s"))) for i in eachindex(dc.du)]
    end
end

test = [1.0u"kg"]

reinterpret(ForwardDiff.Dual{Nothing, Float64, 1}, view([1.0, 1.0], 1:2))

test = PreallocationTools._restructure(test, reinterpret(ForwardDiff.Dual{Nothing, Float64, 1}, view([1.0, 1.0], 1:2)))

view(cached_flat_vec.du, 1:2)

cache = get_tmp(cached_flat_vec, 1.0)
ForwardDiff.can_dual(::Type{Quantity{ForwardDiff.Dual{Nothing, Float64, 1}}}) = false
ForwardDiff.can_dual(::Type{Quantity{Float64}}) = true
cache = get_tmp(cached_flat_vec, ForwardDiff.Dual(1.0u"kg", 1.0))

cache[1]

test = [1.0u"kg"]

test = [1.0, 1.0]

test_diff_cache = DiffCache(test, 2)

test_cache = get_tmp(test_diff_cache, ForwardDiff.Dual(1.0, 1.0))

test_cache[1] + 1.0

ForwardDiff.Dual{nothing}(1.0u"kg", 1.0u"kg/s")



