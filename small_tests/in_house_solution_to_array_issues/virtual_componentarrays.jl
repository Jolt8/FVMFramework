using ComponentArrays
using PreallocationTools
using Polyester
using BenchmarkTools
using SparseConnectivityTracer
import ADTypes

# ============================================================================
# FaceVectorView — view into a flat vector as [cell_id][face_idx]
# ============================================================================

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

# ============================================================================
# VirtualFVMArray — merges multiple buffers into a single logical array
# ============================================================================

struct VirtualAxis{Src, Ax}
    ax::Ax
end

struct VirtualFVMArray{D <: Tuple, A <: NamedTuple}
    data::D
    axes::A
end

# ---------- Resolution logic ----------

@inline _resolve(data, idx::Int) = view(data, idx:idx)
@inline _resolve(data, idx::UnitRange{Int}) = view(data, idx)
@inline _resolve(data, idx::ComponentArrays.ComponentIndex) = ComponentVector(view(data, idx.idx), (idx.ax,))
@inline _resolve(data, ax::ComponentArrays.AbstractAxis) = ComponentVector(data, (ax,))

# Fallback for old Tuple-based face indexing
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

# Support mutation via Symbol indexing (needed for the 'best' looping style)
@inline Base.setindex!(A::VirtualFVMArray, v, s::Symbol) = (getproperty(A, s) .= v)

@inline Base.keys(A::VirtualFVMArray) = keys(getfield(A, :axes))

# Enable .+= and views for VirtualFVMArray
@inline Base.view(A::VirtualFVMArray, s::Symbol) = getproperty(A, s)
@inline Base.dotview(A::VirtualFVMArray, s::Symbol) = getproperty(A, s)

# ---------- Axes Merging ----------

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

# ---------- Field Iteration (matching best_componentarrays_looping_style.jl) ----------

@inline function _resolve_unified_val(A::VirtualFVMArray, ::Val{s}) where {s}
    return getproperty(A, s)
end

"""
    foreach_field_at!(f, cell_id::Int, groups...)
    
Matches the signature in best_componentarrays_looping_style.jl:
Calls f(cell_id, resolved_field1, resolved_field2, ...)
"""
@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    local properties
    if G1 <: VirtualFVMArray
        properties = fieldnames(G1.parameters[2])
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        properties = keys(Ax.parameters[1])
    else
        properties = ()
    end

    exprs = []
    for n in properties
        view_exprs = []
        for i in 1:N
            if groups[i] <: VirtualFVMArray
                push!(view_exprs, :(_resolve_unified_val(groups[$i], Val{$(QuoteNode(n))}())))
            else
                push!(view_exprs, :(getproperty(groups[$i], $(QuoteNode(n)))))
            end
        end
        # Matches best_componentarrays_looping_style.jl: No name passed in cell_id variant
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

"""
    foreach_field_at!(f, groups...)

Matches the signature in best_componentarrays_looping_style.jl:
Calls f(name, groups...)
"""
@generated function foreach_field_at!(f, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    local layout_keys
    if G1 <: VirtualFVMArray
        layout_keys = fieldnames(G1.parameters[2])
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        layout_keys = keys(Ax.parameters[1])
    else
        layout_keys = ()
    end

    exprs = []
    for n in layout_keys
        # Matches best_componentarrays_looping_style.jl: Passes name and the raw groups Tuple
        push!(exprs, :(f($(QuoteNode(n)), groups...)))
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

# ============================================================================
# Test Bench
# ============================================================================

n_cells = 1000
n_faces = 6
reaction_names = (:WGS_rxn, :MD_rxn)
N = 2

# Prototypes for memory layout
du_proto = ComponentVector(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_proto = ComponentVector(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

du_cache_proto = ComponentVector(
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), 
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    rho = zeros(n_cells)
)

u_cache_proto = ComponentVector(
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ),
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    rho = zeros(n_cells)
)

p_proto = ComponentVector(
    diffusion_pre_exponential_factor = [1e-10],
    diffusion_activation_energy = [10000.0]
)

properties_proto = ComponentVector(
    viscosity = zeros(n_cells),
)

# Allocate backing storage
du_vec = Vector(du_proto)
u_vec = Vector(u_proto)
du_cache_vec = Vector(du_cache_proto)
u_cache_vec = Vector(u_cache_proto)
p_vec = Vector(p_proto)
properties_vec = Vector(properties_proto)

# DiffCaches
du_diff_cache_vec = DiffCache(du_cache_vec, N)
u_diff_cache_vec = DiffCache(u_cache_vec, N)

cell_volumes = ones(n_cells)

# Virtual Axes
virtual_du_axes = virtual_merge_axes((du_proto, du_cache_proto))
virtual_u_axes = virtual_merge_axes((u_proto, u_cache_proto, p_proto, properties_proto))

function ode_for_testing_f!(
    du_vec, u_vec, p_vec, t,
    v_du_axes, v_u_axes,
    du_diff_cache, u_diff_cache,
    props_vec,
    cell_vols,
)
    du_vec .= 0.0
    
    # Resolve merged views into VirtualFVMArrays
    if (first(u_vec) + first(p_vec)) isa SparseConnectivityTracer.Dual{Float64}
        u = VirtualFVMArray((u_vec, (get_tmp(u_diff_cache, first(u_vec) + first(p_vec)) .= 0.0), p_vec, props_vec), v_u_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache, first(u_vec) + first(p_vec)) .= 0.0)), v_du_axes)
    else
        u = VirtualFVMArray((u_vec, get_tmp(u_diff_cache, first(u_vec) + first(p_vec)), p_vec, props_vec), v_u_axes)
        du = VirtualFVMArray((du_vec, get_tmp(du_diff_cache, first(u_vec) + first(p_vec))), v_du_axes)
    end

    @batch for cell_id in 1:length(cell_vols)
        # Direct dot syntax
        du.mass_fractions.methylene_blue[cell_id] += 1.0 
        
        # Unrolled iteration over multiple logical arrays
        # Closure args: cid, field1, field2 (Matching best style)
        foreach_field_at!(cell_id, u.mass_fractions, u.molar_concentrations) do cell_id, mass_fractions, molar_concentrations
            # mf and mc are already resolved views
            mass_fractions.methylene_blue[cell_id] += 1.0
            molar_concentrations.methylene_blue[cell_id] += 1.0
        end
        
        # Unrolled iteration over scalars/sub-caches
        # Closure args: name, groups... (Matching best style)
        foreach_field_at!(u.net_rates.reforming_reactions) do reaction_name, reforming_net_rates
            # Match best_componentarrays_looping_style.jl: Use indexing via the symbol
            reforming_net_rates[reaction_name] += 1.0
        end
        
        # Cross-buffer access
        du.mass[cell_id] += du.rho[cell_id] * u.viscosity[cell_id]
    end
end

# Benchmark
du_view = du_vec
u_view = u_vec

println("Benchmarking ode_for_testing_f! with 'Best' Looping Style...")
@btime ode_for_testing_f!(
    $du_view, $u_view, $p_vec, 0.0,
    $virtual_du_axes, $virtual_u_axes,
    $du_diff_cache_vec, $u_diff_cache_vec, 
    $properties_vec,
    $cell_volumes,
)

using OrdinaryDiffEq

f_closure = (du, u, p, t) -> ode_for_testing_f!(
    du, u, p, t,

    virtual_du_axes, virtual_u_axes,

    du_diff_cache_vec, u_diff_cache_vec, 

    properties_vec,
    
    cell_volumes
)

length(du_vec)
length(u_vec)

jac_sparsity = SparseConnectivityTracer.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_vec, 0.0), du_vec, u_vec, SparseConnectivityTracer.TracerLocalSparsityDetector()
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 1.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u_vec, tspan, p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
VSCodeServer.@profview sol = solve(implicit_prob, FBDF())

u_vec .= 0.0
explicit_prob = ODEProblem(f_closure, u_vec, tspan, p_vec)
VSCodeServer.@profview sol = solve(explicit_prob, Tsit5())