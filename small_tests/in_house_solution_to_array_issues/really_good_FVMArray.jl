using ComponentArrays 
using PreallocationTools
using Polyester
using LazyArrays
using BenchmarkTools
using SparseConnectivityTracer
import ADTypes
using OrdinaryDiffEq
include("FVMArray.jl")

struct VirtualAxis{Src, Ax}
    ax::Ax
end
struct VirtualFVMArray{D <: Tuple, A <: NamedTuple}
    data::D
    axes::A
end
# Resolve logic: Redirects to the correct buffer in the Tuple
@inline function _virtual_resolve(data::Tuple, vax::VirtualAxis{Src}) where {Src}
    return _resolve(data[Src], vax.ax) # Uses resolution from FVMArray.jl
end

@inline function _virtual_resolve(data::Tuple, ax::NamedTuple)
    return VirtualFVMArray(data, ax)
end

@inline function Base.getproperty(A::VirtualFVMArray, s::Symbol)
    ax = getfield(getfield(A, :axes), s)
    return _virtual_resolve(getfield(A, :data), ax)
end

@inline Base.getindex(A::VirtualFVMArray, s::Symbol) = getproperty(A, s)
@inline Base.keys(A::VirtualFVMArray) = keys(getfield(A, :axes))

# Helper to transform standard axes into Virtual ones
function virtual_merge_axes(axes_list::Tuple)
    merged_entries = Symbol[]
    merged_values = Any[]
    for (i, ax_nt) in enumerate(axes_list)
        for name in keys(ax_nt)
            push!(merged_entries, name)
            push!(merged_values, _wrap_virtual(ax_nt[name], i))
        end
    end
    return NamedTuple{Tuple(merged_entries)}(Tuple(merged_values))
end

_wrap_virtual(ax::NamedTuple, src::Int) = NamedTuple{keys(ax)}(map(v -> _wrap_virtual(v, src), values(ax)))
_wrap_virtual(ax, src::Int) = VirtualAxis{src, typeof(ax)}(ax)

# Required for @batch compatibility
@inline function _resolve_unified_val(A::VirtualFVMArray, ::Val{s}) where {s}
    return getproperty(A, s)
end

@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{VirtualFVMArray, N}) where {N}
    field_names = fieldnames(groups[1].parameters[2])
    exprs = []
    for n in field_names
        view_exprs = [:(_resolve_unified_val(groups[$i], Val{$(QuoteNode(n))}())) for i in 1:N]
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

function ode_for_testing_f!(
    du_vec, u_vec, p_vec, t,

    virtual_du_axes, virtual_u_axes,

    du_diff_cache, u_diff_cache,

    properties_vec,
    
    cell_volumes,
)
    if (first(u_vec) + first(p_vec)) isa SparseConnectivityTracer.Dual{Float64, SparseConnectivityTracer.GradientTracer{Int64, BitSet}}
        u = VirtualFVMArray((u_vec, (get_tmp(u_diff_cache, first(u_vec) + first(p_vec)) .= 0.0), properties_vec, p_vec), virtual_u_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache, first(u_vec) + first(p_vec)) .= 0.0)), virtual_du_axes)
    else
        u = VirtualFVMArray((u_vec, get_tmp(u_diff_cache, first(u_vec) + first(p_vec)), properties_vec, p_vec), virtual_u_axes)
        du = VirtualFVMArray((du_vec, get_tmp(du_diff_cache, first(u_vec) + first(p_vec))), virtual_du_axes)
    end

    @batch for cell_id in 1:length(cell_volumes)
        du.mass_fractions[:methylene_blue][cell_id] += 1.0 
        du.mass_fractions.methylene_blue[cell_id] += 1.0 

        u.net_rates.reforming_reactions.WGS_rxn[1] += 1.0 

        foreach_field_at!(cell_id, du.mass_fractions) do species_name, mass_fraction
            mass_fraction[species_name] += 1.0 
            #foreach_field_at!(1, du.net_rates.reforming_reactions) do reaction_name, net_rate #this allocates a ton when doing @batch
                #net_rate[reaction_name] += 1.0 
            #end
        end

        foreach_field_at!(1, du.net_rates.reforming_reactions) do reaction_name, net_rate #oh, this is fine when used outside another foreach_field_at!
            net_rate[reaction_name] += 1.0 
        end

        for face_idx in 1:6
            du.mass_face[cell_id][face_idx] += 1.0 
        end

        du.mass[cell_id] += sum(du.mass_face[cell_id]) 
    end
end

#=
other things to keep in mind
    - offset indexing like du.mass_face[cell_id + 1][5]
    - LazyMerging of du_in, cache, and properties 
    - ability to store interpolations like where u.viscosity would call u.viscosity(u.temp, u.pressure) under the hood
    - it's probably going to be called FVMArray or maybe something even more generic because I could see other people using this
=#

n_cells = 100000
n_faces = 6
reaction_names = (:WGS_rxn, :MD_rxn)
N = 2

du_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

du_axes, du_offset = create_axes(du_proto)
du_vec = zeros(du_offset)

u_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_axes, u_offset = create_axes(u_proto)
u_vec = zeros(u_offset)

du_caches = (
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

du_caches_axes, du_caches_offset = create_axes(du_caches)
du_caches_vec = zeros(du_caches_offset)
du_diff_cache_vec = DiffCache(du_caches_vec, N)

u_caches = (
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

u_caches_axes, u_caches_offset = create_axes(u_caches)
u_caches_vec = zeros(u_caches_offset)
u_diff_cache_vec = DiffCache(u_caches_vec, N)

p = (
    diffusion_pre_exponential_factor = [1e-10],
    diffusion_activation_energy = [10000.0]
)

p_axes, p_offset = create_axes(p)
p_vec = zeros(p_offset)

properties = (
    viscosity = zeros(n_cells),
    #diffusion_pre_exponential_factor = zeros(n_cells),
    #diffusion_activation_energy = zeros(n_cells)
)

properties_axes, properties_offset = create_axes(properties)
properties_vec = zeros(properties_offset)

du_merged_buffer = create_merged_buffer(length(du_vec) + length(du_caches_vec))
du_merged_axes = merge_axes((du_axes, du_caches_axes), (du_offset, du_caches_offset))
u_merged_buffer = create_merged_buffer(length(u_vec) + length(u_caches_vec) + length(properties_vec) + length(p_vec))
u_merged_axes = merge_axes((u_axes, u_caches_axes, properties_axes, p_axes), (u_offset, u_caches_offset, properties_offset, p_offset))

cell_volumes = ones(n_cells)

du_view = view(du_merged_buffer, 1:length(du_vec))
u_view = view(u_merged_buffer, 1:length(u_vec))

du_merged_diff_cache = DiffCache(du_merged_buffer, N)
u_merged_diff_cache = DiffCache(u_merged_buffer, N)

caches_unitrange = length(du_vec) + 1 : length(du_vec) + length(du_caches_vec)
property_unitrange = length(u_vec) + length(u_caches_vec) + 1 : length(u_vec) + length(u_caches_vec) + length(properties_vec)
p_unitrange = length(u_vec) + length(u_caches_vec) + length(properties_vec) + 1 : length(u_vec) + length(u_caches_vec) + length(properties_vec) + length(p_vec)

virtual_du_axes = virtual_merge_axes((du_axes, du_caches_axes))
virtual_u_axes = virtual_merge_axes((u_axes, u_caches_axes, properties_axes, p_axes))

test = vcat(du_vec, u_vec)

du_view = view(du_vec, 1:length(du_vec)) #for some reason this is faster than du_view = view(test, 1:length(du_vec))
u_view = view(u_vec, 1:length(u_vec))

@btime ode_for_testing_f!(
    $du_view, $u_view, $p_vec, 0.0,
    
    $virtual_du_axes, $virtual_u_axes,

    $du_diff_cache_vec, $u_diff_cache_vec, 

    $properties_vec,
    
    $cell_volumes,
)
#this takes 22.8 μs with 37 allocations and 2.62 KiB with 1000 cells
#this takes 201.1 μs with 37 allocations and 2.62 KiB with 10000 cells (this is with @batch)
#this takes 267.3 μs with 37 allocations and 2.62 KiB with 10000 cells (this is without @batch)
    #I don't know why @batch isn't speeding this up that much
    #Julia.numthreads is set to auto
#this takes 267.3 μs with 37 allocations and 2.62 KiB with 1000 cells (this is without @batch)
#4.929 ms (37 allocations: 2.62 KiB) with 100000 cells (without @batch)
#4.223 ms (37 allocations: 2.62 KiB) with 100000 cells (with @batch)
#5.554 (22 allocations, 4.58 MiB) with 100000 cells (with DiffCached du_merged and u_merged)
#2.175 ms (1 allocations, 288 bytes) 

ode_for_testing_f!(
    du_view, u_view, p_vec, 0.0,
    
    virtual_du_axes, virtual_u_axes,

    du_diff_cache_vec, u_diff_cache_vec, 

    properties_vec,
    
    cell_volumes
)

using ForwardDiff

zoop(x) = ode_for_testing_f!(
    du_view, u_view, [x[1], x[2]], 0.0,
    
    virtual_du_axes, virtual_u_axes,

    du_diff_cache_vec, u_diff_cache_vec, 

    properties_vec,
    
    cell_volumes
)

ForwardDiff.gradient(zoop, [1.0, 2.0])

function test_raw_loop_performance(du_merged_buffer, u_merged_buffer)
    for i in eachindex(du_merged_buffer)
        du_merged_buffer[i] += 1.0
        u_merged_buffer[i] += 1.0
    end
end

#@btime test_raw_loop_performance(du_merged_buffer, u_merged_buffer)
#1.308 ms (0 allocations: 0 bytes) for 100000 cells (with @batch)
#1.086 ms (0 allocations: 0 bytes) for 100000 cells (without @batch)

f_closure = (du, u, p, t) -> ode_for_testing_f!(
    du, u, p, t,

    virtual_du_axes, virtual_u_axes,

    du_diff_cache_vec, u_diff_cache_vec, 

    properties_vec,
    
    cell_volumes
)

jac_sparsity = SparseConnectivityTracer.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_vec, 0.0), du_view, u_view, SparseConnectivityTracer.TracerLocalSparsityDetector()
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 100000.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u_view, tspan, p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF())

u_view .= 0.0
explicit_prob = ODEProblem(f_closure, u_view, tspan, p_vec)
@time sol = solve(explicit_prob, Tsit5())

u_named = FVMArray(sol.u[1], u_axes)

findall(x -> x > 0.0, u_named)

u_named.mass_fractions.methylene_blue
