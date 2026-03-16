#=
    FVMArray.jl — Custom array type for Finite Volume Method simulations
    
    Replaces ComponentArrays with a type that:
      - Returns views (not copies) for symbolic indexing → mutations propagate
      - Supports nested named fields, face-indexed sub-arrays, and scalar fields
      - Eagerly merges multiple vectors into a single pre-allocated dense buffer
      - Is compatible with DiffCache, SparseConnectivityTracer, @batch, and ODE solvers
=#

# ============================================================================
# FaceVectorView — view into a flat vector as [cell_id][face_idx]
# ============================================================================

struct FaceVectorView{V<:AbstractVector} <: AbstractVector{SubArray}
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
# FVMArray — the main type
# ============================================================================

struct FVMArray{T, D<:AbstractVector{T}, A<:NamedTuple} <: AbstractVector{T}
    data::D
    axes::A
end

# ---------- AbstractVector interface ----------

@inline Base.size(A::FVMArray) = size(getfield(A, :data))
@inline Base.length(A::FVMArray) = length(getfield(A, :data))

@inline function Base.getindex(A::FVMArray, i::Int)
    @boundscheck checkbounds(getfield(A, :data), i)
    return @inbounds getfield(A, :data)[i]
end

@inline function Base.setindex!(A::FVMArray, v, i::Int)
    @boundscheck checkbounds(getfield(A, :data), i)
    @inbounds getfield(A, :data)[i] = v
    return v
end

Base.IndexStyle(::Type{<:FVMArray}) = IndexLinear()

# similar returns a plain Vector so ODE solvers, broadcast, etc. don't need to carry axes
Base.similar(A::FVMArray) = similar(getfield(A, :data))
Base.similar(A::FVMArray, ::Type{S}) where {S} = similar(getfield(A, :data), S)
Base.similar(A::FVMArray, ::Type{S}, dims::Dims) where {S} = similar(getfield(A, :data), S, dims)

# ---------- Internal resolution: axes entry → view/sub-FVMArray ----------

# Cell-indexed field (UnitRange) → view
@inline _resolve(data::AbstractVector, ax::UnitRange{Int}) = view(data, ax)

# Nested group (NamedTuple) → new FVMArray sharing the same backing data
@inline _resolve(data::AbstractVector, ax::NamedTuple) = FVMArray(data, ax)

# Face-indexed field (Tuple{Int,Int}) → FaceVectorView
@inline _resolve(data::AbstractVector, ax::Tuple{Int,Int}) = FaceVectorView(data, ax[2], ax[1])

# Scalar field (single-element range stored as Tuple{Int}) → view of 1 element so += works
@inline _resolve(data::AbstractVector, ax::Tuple{Int}) = view(data, ax[1]:ax[1])

# ---------- Symbolic access (dot notation and [:symbol]) ----------

@inline function Base.getproperty(A::FVMArray, s::Symbol)
    return _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
end

@inline Base.propertynames(A::FVMArray) = keys(getfield(A, :axes))

@inline function Base.getindex(A::FVMArray, s::Symbol)
    return _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
end

# ---------- dotview / view for Symbol indexing (so .+= works) ----------

@inline Base.dotview(A::FVMArray, s::Symbol) = _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
@inline Base.view(A::FVMArray, s::Symbol) = _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))

# ---------- keys for allocation-free iteration ----------

@inline Base.keys(A::FVMArray) = keys(getfield(A, :axes))

# ---------- Broadcasting ----------

Base.BroadcastStyle(::Type{<:FVMArray}) = Broadcast.DefaultArrayStyle{1}()

# dataids so broadcast can detect aliasing with the underlying data
Base.dataids(A::FVMArray) = Base.dataids(getfield(A, :data))

# ============================================================================
# @batch-compatible field iteration (compile-time unrolled)
# ============================================================================

# Val-based resolve for compile-time specialization
@inline function _resolve_val(A::FVMArray, ::Val{s}) where {s}
    return _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
end

"""
    field_views(A::FVMArray) -> Tuple

Return a tuple of views for all fields in `A`, fully typed at compile time.
Useful for pre-computing views outside `@batch` loops.

# Example
```julia
mf_views = field_views(u.mass_fractions)
mc_views = field_views(u.molar_concentrations)
@batch for cell_id in 1:n_cells
    for (mf, mc) in zip(mf_views, mc_views)
        mf[cell_id] += 1.0
        mc[cell_id] += 1.0
    end
end
```
"""
@generated function field_views(A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(_resolve_val(A, Val{$(QuoteNode(n))}())) for n in field_names]
    return quote
        Base.@_inline_meta
        ($(exprs...),)
    end
end

"""
    foreach_field_at!(f, cell_id::Int, groups::FVMArray...)

Call `f(cell_id, view1, view2, ...)` for each field, where `view1, view2, ...`
are the resolved views from each group. All groups must share the same field names.
Fully unrolled at compile time — safe to use inside `@batch`.

# Example
```julia
@batch for cell_id in 1:n_cells
    foreach_field_at!(cell_id, u.mass_fractions, u.molar_concentrations) do cid, mf, mc
        mf[cid] += 1.0
        mc[cid] += 1.0
    end
end
```
"""
@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{FVMArray, N}) where {N}
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

"""
    batch_foreach_field(f, A::FVMArray)

Call `f(name::Symbol, view)` for each field in `A`.
Fully unrolled at compile time — safe to use inside `@batch`.

# Example
```julia
@batch for cell_id in 1:n_cells
    batch_foreach_field(u.mass_fractions) do name, view
        view[cell_id] += 1.0
    end
end
```
"""
@generated function batch_foreach_field(f, A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(f($(QuoteNode(n)), _resolve_val(A, Val{$(QuoteNode(n))}()))) for n in field_names]
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

# ============================================================================
# Axes builder — creates the compile-time NamedTuple of indices from a prototype
# ============================================================================

"""
    create_axes(prototype::NamedTuple; offset::Int=0) -> (axes::NamedTuple, total_length::Int)

Walk a prototype NamedTuple (same shape you'd pass to ComponentArray) and produce
an axes NamedTuple mapping each field to its index range in a flat vector.

Supported leaf types:
  - `Vector{<:Number}` of length n → UnitRange (cell-indexed, length n)
  - `Vector{Vector{<:Number}}` of length n_cells, inner length n_faces → Tuple{Int,Int} (face-indexed)
  - `Number` or length-1 Vector  → Tuple{Int} (scalar)
  - `NamedTuple` → recurse
"""
function create_axes(prototype::NamedTuple; offset::Int=0)
    entries = Symbol[]
    values = Any[]
    current_offset = offset

    for name in keys(prototype)
        val = prototype[name]
        ax_entry, current_offset = _axes_entry(val, current_offset)
        push!(entries, name)
        push!(values, ax_entry)
    end

    axes_nt = NamedTuple{Tuple(entries)}(Tuple(values))
    return axes_nt, current_offset
end

# Dispatch helpers for create_axes

# NamedTuple → recurse
function _axes_entry(val::NamedTuple, offset::Int)
    sub_axes, new_offset = create_axes(val; offset=offset)
    return sub_axes, new_offset
end

# Vector of numbers → cell-indexed UnitRange
function _axes_entry(val::AbstractVector{<:Number}, offset::Int)
    n = length(val)
    if n == 1
        # Single-element vector treated as scalar
        return (offset + 1,), offset + 1
    else
        return (offset + 1):(offset + n), offset + n
    end
end

# Vector of vectors → face-indexed
function _axes_entry(val::AbstractVector{<:AbstractVector}, offset::Int)
    n_cells = length(val)
    n_faces = length(val[1])
    start = offset + 1
    total = n_cells * n_faces
    return (start, n_faces), offset + total
end

# Scalar number
function _axes_entry(val::Number, offset::Int)
    return (offset + 1,), offset + 1
end

# ============================================================================
# Eager merge — pre-allocated in-place copyto!
# ============================================================================

"""
    create_merged_buffer(total_length::Int, T::Type=Float64) -> Vector{T}

Pre-allocate a dense buffer for merging multiple source vectors.
"""
create_merged_buffer(total_length::Int, ::Type{T}=Float64) where {T} = zeros(T, total_length)

"""
    merge_into!(buffer, vecs::Tuple)

Copy each vector in `vecs` contiguously into `buffer` (in order).
Zero-allocation after the first call (no resizing).
"""
@inline function merge_into!(buffer::AbstractVector, vecs::Tuple)
    offset = 0
    for v in vecs
        n = length(v)
        @inbounds copyto!(buffer, offset + 1, v, 1, n)
        offset += n
    end
    return buffer
end

"""
    merge_axes(axes_tuple::Tuple, lengths::Tuple) -> NamedTuple

Merge multiple axes NamedTuples, shifting each by the cumulative offset of
the previous source vectors' lengths. Returns a single combined axes NamedTuple.
"""
function merge_axes(axes_tuple::Tuple, lengths::Tuple)
    @assert length(axes_tuple) == length(lengths)
    entries = Symbol[]
    values = Any[]
    offset = 0

    for (i, ax) in enumerate(axes_tuple)
        _collect_shifted!(entries, values, ax, offset)
        offset += lengths[i]
    end

    return NamedTuple{Tuple(entries)}(Tuple(values))
end

function _collect_shifted!(entries, values, ax::NamedTuple, offset::Int)
    for name in keys(ax)
        if name in entries
            error("Duplicate field name during merge: $name. Rename one of the fields.")
        end
        push!(entries, name)
        push!(values, _shift_ax(ax[name], offset))
    end
end

# Shift each axes leaf type by an offset
_shift_ax(ax::UnitRange{Int}, offset::Int) = (ax.start + offset):(ax.stop + offset)
_shift_ax(ax::Tuple{Int,Int}, offset::Int) = (ax[1] + offset, ax[2])      # face: shift start, keep n_faces
_shift_ax(ax::Tuple{Int}, offset::Int) = (ax[1] + offset,)                # scalar: shift index
function _shift_ax(ax::NamedTuple, offset::Int)
    shifted_entries = Symbol[]
    shifted_values = Any[]
    for name in keys(ax)
        push!(shifted_entries, name)
        push!(shifted_values, _shift_ax(ax[name], offset))
    end
    return NamedTuple{Tuple(shifted_entries)}(Tuple(shifted_values))
end

# ============================================================================
# Convenience constructor
# ============================================================================

"""
    FVMArray(data::AbstractVector, axes::NamedTuple)

Wrap an existing vector with the given axes layout.
"""
FVMArray(data::AbstractVector, axes::NamedTuple) = FVMArray{eltype(data), typeof(data), typeof(axes)}(data, axes)

"""
    FVMArray(prototype::NamedTuple)

Create an FVMArray from a prototype NamedTuple, allocating a new backing vector.
"""
function FVMArray(prototype::NamedTuple)
    axes, total_len = create_axes(prototype)
    data = zeros(Float64, total_len)
    return FVMArray(data, axes)
end

# ============================================================================
# Pretty printing
# ============================================================================

function Base.show(io::IO, A::FVMArray)
    print(io, "FVMArray(")
    print(io, length(A), " elements, fields: ", join(keys(A), ", "))
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", A::FVMArray)
    println(io, "FVMArray with $(length(A)) elements")
    _show_axes(io, getfield(A, :axes), "  ")
end

function _show_axes(io::IO, ax::NamedTuple, indent::String)
    for name in keys(ax)
        val = ax[name]
        if val isa NamedTuple
            println(io, indent, name, " (group):")
            _show_axes(io, val, indent * "  ")
        elseif val isa UnitRange
            println(io, indent, name, " → ", length(val), " cells (indices ", val, ")")
        elseif val isa Tuple{Int,Int}
            n_cells = 0  # can't know without data length, just show the tuple
            println(io, indent, name, " → face-indexed (start=", val[1], ", n_faces=", val[2], ")")
        elseif val isa Tuple{Int}
            println(io, indent, name, " → scalar (index ", val[1], ")")
        else
            println(io, indent, name, " → ", val)
        end
    end
end


# ============================================================================
# Tests
# ============================================================================

function __test__()
    println("=" ^ 60)
    println("FVMArray Tests")
    println("=" ^ 60)

    n_cells = 1000
    n_faces = 6
    passed = 0
    failed = 0

    function check(label, condition)
        if condition
            println("  ✓ $label")
            passed += 1
        else
            println("  ✗ FAIL: $label")
            failed += 1
        end
    end

    # ---- Build prototypes matching requirements_example.jl ----

    du_proto = (
        mass_fractions = (
            methylene_blue = zeros(n_cells),
            water = zeros(n_cells),
        ),
        pressure = zeros(n_cells),
    )

    u_proto = (
        mass_fractions = (
            methylene_blue = zeros(n_cells),
            water = zeros(n_cells),
        ),
        pressure = zeros(n_cells),
    )

    du_cache_proto = (
        mass = zeros(n_cells),
        mass_face = fill(zeros(n_faces), n_cells),
        net_rates = (
            reforming_reactions = (
                WGS_rxn = 0.0,
                MD_rxn = 0.0,
            ),
        ),
        rho = zeros(n_cells),
    )

    u_cache_proto = (
        mass = zeros(n_cells),
        mass_face = fill(zeros(n_faces), n_cells),
        net_rates = (
            reforming_reactions = (
                WGS_rxn = 0.0,
                MD_rxn = 0.0,
            ),
        ),
        rho = zeros(n_cells),
    )

    properties_proto = (
        rho_prop = zeros(n_cells),
        viscosity = zeros(n_cells),
        diffusion_pre_exponential_factor = zeros(n_cells),
        diffusion_activation_energy = zeros(n_cells),
    )

    # ---- Test 1: create_axes ----
    println("\n[1] create_axes")
    du_axes, du_len = create_axes(du_proto)
    check("du_axes has correct keys", keys(du_axes) == (:mass_fractions, :pressure))
    check("du length = 3 * n_cells = $(3 * n_cells)", du_len == 3 * n_cells)

    du_cache_axes, du_cache_len = create_axes(du_cache_proto)
    check("du_cache_axes has correct keys", keys(du_cache_axes) == (:mass, :mass_face, :net_rates, :rho))
    expected_cache_len = n_cells + n_cells * n_faces + 2 + n_cells  # mass + face + 2 scalars + rho
    check("du_cache length = $expected_cache_len", du_cache_len == expected_cache_len)

    # ---- Test 2: FVMArray construction ----
    println("\n[2] FVMArray construction")
    du_data = zeros(du_len)
    du = FVMArray(du_data, du_axes)
    check("FVMArray length matches", length(du) == du_len)
    check("FVMArray integer indexing works", du[1] == 0.0)

    # ---- Test 3: Dot property access returns views (not copies) ----
    println("\n[3] Property access returns views")
    du_data .= 0.0
    du.mass_fractions.methylene_blue[5] = 42.0
    check("dot access mutates backing data", du_data[5] == 42.0)

    du.pressure[10] = 99.0
    check("pressure dot access mutates backing data", du_data[2 * n_cells + 10] == 99.0)

    # ---- Test 4: Symbol indexing returns views (not copies) ----
    println("\n[4] Symbol indexing returns views")
    du_data .= 0.0
    du[:mass_fractions][:methylene_blue][7] = 55.0
    check("symbol indexing mutates backing data", du_data[7] == 55.0)

    # ---- Test 5: keys for looping ----
    println("\n[5] keys iteration")
    species_keys = keys(du.mass_fractions)
    check("keys returns correct names", species_keys == (:methylene_blue, :water))

    du_data .= 0.0
    for species_name in keys(du.mass_fractions)
        du.mass_fractions[species_name][1] += 1.0
    end
    check("looping via keys mutates correctly", du_data[1] == 1.0 && du_data[n_cells + 1] == 1.0)

    # ---- Test 6: FaceVectorView ----
    println("\n[6] FaceVectorView")
    cache_data = zeros(du_cache_len)
    du_cache = FVMArray(cache_data, du_cache_axes)
    du_cache.mass_face[3][2] = 7.7
    check("face indexing mutates backing data", any(x -> x == 7.7, cache_data))
    check("sum of face view works", sum(du_cache.mass_face[3]) == 7.7)

    for face_idx in 1:n_faces
        du_cache.mass_face[1][face_idx] = Float64(face_idx)
    end
    check("face loop writes correctly", sum(du_cache.mass_face[1]) == sum(1:n_faces))

    # ---- Test 7: Nested scalar fields ----
    println("\n[7] Nested scalar fields")
    cache_data .= 0.0
    du_cache.net_rates.reforming_reactions.WGS_rxn .+= 3.14
    check("scalar field += works", du_cache.net_rates.reforming_reactions.WGS_rxn[1] ≈ 3.14)

    # ---- Test 8: Eager merge ----
    println("\n[8] Eager merge (pre-allocated)")
    
    u_axes, u_len = create_axes(u_proto)
    u_cache_axes, u_cache_len = create_axes(u_cache_proto)
    prop_axes, prop_len = create_axes(properties_proto)

    # Create source vectors
    u_src = ones(u_len)
    u_cache_src = 2.0 .* ones(u_cache_len)
    prop_src = 3.0 .* ones(prop_len)

    # Merge axes
    merged_u_axes = merge_axes(
        (prop_axes, u_axes, u_cache_axes),
        (prop_len, u_len, u_cache_len),
    )
    total_u_len = prop_len + u_len + u_cache_len

    # Pre-allocate buffer
    merged_buf = create_merged_buffer(total_u_len)
    merge_into!(merged_buf, (prop_src, u_src, u_cache_src))

    u_merged = FVMArray(merged_buf, merged_u_axes)

    check("merged length correct", length(u_merged) == total_u_len)
    check("properties region has value 3.0", u_merged.rho_prop[1] == 3.0)
    check("u region has value 1.0", u_merged.mass_fractions.methylene_blue[1] == 1.0)
    check("cache region has value 2.0", u_merged.mass[1] == 2.0)

    # ---- Test 9: Nested key iteration ----
    println("\n[9] Nested key iteration")
    for species_name in keys(u_merged.mass_fractions)
        u_merged.mass_fractions[species_name][1] += 10.0
    end
    check("nested key loop mutates merged buffer", merged_buf[prop_len + 1] == 11.0)

    nested_keys = keys(u_merged.net_rates.reforming_reactions)
    check("deep nested keys work", nested_keys == (:WGS_rxn, :MD_rxn))

    for reaction in keys(u_merged.net_rates.reforming_reactions)
        u_merged.net_rates.reforming_reactions[reaction] .+= 1.0
    end
    check("deep nested key loop mutates correctly",
        u_merged.net_rates.reforming_reactions.WGS_rxn[1] ≈ 3.0 &&
        u_merged.net_rates.reforming_reactions.MD_rxn[1] ≈ 3.0)

    # ---- Test 10: merge_into! is allocation-free after first call ----
    println("\n[10] merge_into! allocation check")
    merge_into!(merged_buf, (prop_src, u_src, u_cache_src))  # warm up
    allocs = @allocated merge_into!(merged_buf, (prop_src, u_src, u_cache_src))
    check("merge_into! is allocation-free", allocs == 0)

    # ---- Test 11: keys looping allocation check ---- (GARBAGE PATTERN!)
    println("\n[11] Keys looping allocation check")
    # Warm up
    map(keys(u_merged.mass_fractions)) do species
        u_merged.mass_fractions[species][1] += 0.0
    end
    allocs_keys = @allocated begin
        map(keys(u_merged.mass_fractions)) do species
            u_merged.mass_fractions[species][1] += 0.0
        end
    end
    check("keys looping is allocation-free", allocs_keys == 0)

    # ---- Test 12: .= broadcast ----
    println("\n[12] Broadcast .= ")
    du_data .= 0.0
    du.pressure .= 5.0
    check("broadcast .= on field works", all(du.pressure .== 5.0))

    # ---- Test 13: foreach_field_at! correctness ----
    println("\n[13] foreach_field_at! correctness")
    
    # Build a second prototype with matching keys
    dual_proto = (
        mass_fractions = (
            methylene_blue = zeros(n_cells),
            water = zeros(n_cells),
        ),
    )
    dual_axes, dual_len = create_axes(dual_proto)
    
    data_a = zeros(dual_len)
    data_b = zeros(dual_len)
    arr_a = FVMArray(data_a, dual_axes)
    arr_b = FVMArray(data_b, dual_axes)
    
    for cell_id in 1:n_cells
        foreach_field_at!(cell_id, arr_a.mass_fractions, arr_b.mass_fractions) do cid, va, vb
            va[cid] += 2.0
            vb[cid] += 3.0
        end
    end
    check("foreach_field_at! mutates arr_a", all(data_a .== 2.0))
    check("foreach_field_at! mutates arr_b", all(data_b .== 3.0))
    
    # Test field_views
    println("\n[14] field_views correctness")
    data_a .= 0.0
    data_b .= 0.0
    views_a = field_views(arr_a.mass_fractions)
    views_b = field_views(arr_b.mass_fractions)
    for (va, vb) in zip(views_a, views_b)
        for cell_id in 1:n_cells
            va[cell_id] += 5.0
            vb[cell_id] += 7.0
        end
    end
    check("field_views mutates arr_a", all(data_a .== 5.0))
    check("field_views mutates arr_b", all(data_b .== 7.0))
    
    # Test batch_foreach_field
    println("\n[15] batch_foreach_field correctness")
    data_a .= 0.0
    for cell_id in 1:n_cells
        batch_foreach_field(arr_a.mass_fractions) do name, view
            view[cell_id] += 11.0
        end
    end
    check("batch_foreach_field mutates correctly", all(data_a .== 11.0))

    # ---- Summary ----
    println("\n" * "=" ^ 60)
    println("Results: $passed passed, $failed failed out of $(passed + failed)")
    println("=" ^ 60)

    return failed == 0
end

# Run tests when included directly
if abspath(PROGRAM_FILE) == @__FILE__
    __test__()
end

__test__()