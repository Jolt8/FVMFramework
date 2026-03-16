#=
    Experiments for allocation-free looping over FVMArray fields by key.

    The core problem: `for s in keys(fvmarray.group)` loops over a Tuple of Symbols.
    Inside the loop body, `fvmarray.group[s]` dispatches on `s::Symbol`, which Julia
    can't specialize because the *return type* varies per symbol (SubArray vs FVMArray 
    vs FaceVectorView). This causes dynamic dispatch → allocations.

    Below we try several approaches to eliminate this.
    (No BenchmarkTools dependency — uses @allocated and @elapsed only)
=#

include("FVMArray.jl")

# ---- Setup ----
n_cells = 1000
proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
        ethanol = zeros(n_cells),
    ),
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
        ethanol = zeros(n_cells),
    ),
    pressure = zeros(n_cells),
)

ax, len = create_axes(proto)
data = zeros(len)
u = FVMArray(data, ax)

# Helper to benchmark
function bench(label, f!, u, n_cells; warmup=3, runs=100)
    for _ in 1:warmup; f!(u, n_cells); end
    allocs = @allocated f!(u, n_cells)
    
    total = 0.0
    for _ in 1:runs; total += @elapsed f!(u, n_cells); end
    avg_ns = total / runs * 1e9

    status = allocs == 0 ? "✓ 0 allocs" : "✗ $allocs bytes"
    println("  $label: $status, avg $(round(avg_ns, digits=1)) ns")
    return allocs
end

# ============================================================================
# Approach 0: Naive for loop (BASELINE)
# ============================================================================

function loop_naive!(u, n_cells)
    for species in keys(u.mass_fractions)
        for cell_id in 1:n_cells
            u.mass_fractions[species][cell_id] += 1.0
            u.molar_concentrations[species][cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach 1: map(f, keys(nt))
# ============================================================================

function loop_map!(u, n_cells)
    map(keys(u.mass_fractions)) do species
        for cell_id in 1:n_cells
            u.mass_fractions[species][cell_id] += 1.0
            u.molar_concentrations[species][cell_id] += 1.0
        end
    end
    nothing
end

# ============================================================================
# Approach 2: foreach(f, keys(nt))
# ============================================================================

function loop_foreach!(u, n_cells)
    foreach(keys(u.mass_fractions)) do species
        for cell_id in 1:n_cells
            u.mass_fractions[species][cell_id] += 1.0
            u.molar_concentrations[species][cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach 3: @generated foreach_field — compile-time unrolling via Val dispatch
# ============================================================================

@inline function _resolve_val(A::FVMArray, ::Val{s}) where {s}
    return _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
end

@generated function foreach_field(f, A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(f(_resolve_val(A, Val{$(QuoteNode(n))}()))) for n in field_names]
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

@generated function foreach_field_named(f, A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(f($(QuoteNode(n)), _resolve_val(A, Val{$(QuoteNode(n))}()))) for n in field_names]
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

function loop_generated!(u, n_cells)
    foreach_field(u.mass_fractions) do field_view
        for cell_id in 1:n_cells
            field_view[cell_id] += 1.0
        end
    end
end

function loop_generated_named!(u, n_cells)
    foreach_field_named(u.mass_fractions) do name, field_view
        for cell_id in 1:n_cells
            field_view[cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach 4: Get all views as a tuple, then iterate the tuple
# ============================================================================

@generated function field_views(A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(_resolve_val(A, Val{$(QuoteNode(n))}())) for n in field_names]
    return quote
        Base.@_inline_meta
        ($(exprs...),)
    end
end

function loop_tuple!(u, n_cells)
    views = field_views(u.mass_fractions)
    for v in views
        for cell_id in 1:n_cells
            v[cell_id] += 1.0
        end 
    end
end

# ============================================================================
# Approach 5: @fvm_foreach macro (wraps @generated foreach_field)
# ============================================================================

macro fvm_foreach(group_expr, var_name, body)
    quote
        foreach_field($(esc(group_expr))) do $(esc(var_name))
            $(esc(body))
        end
    end
end

macro fvm_foreach_named(group_expr, name_var, view_var, body)
    quote
        foreach_field_named($(esc(group_expr))) do $(esc(name_var)), $(esc(view_var))
            $(esc(body))
        end
    end
end

function loop_macro!(u, n_cells)
    @fvm_foreach u.mass_fractions species_view begin
        for cell_id in 1:n_cells
            species_view[cell_id] += 1.0
        end
    end
end

function loop_macro_named!(u, n_cells)
    @fvm_foreach_named u.mass_fractions species_name species_view begin
        for cell_id in 1:n_cells
            species_view[cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach 6: Recursive tuple head/tail (no @generated, pure Julia)
# ============================================================================

@inline _foreach_recurse(f, A::FVMArray, ::Tuple{}) = nothing
@inline function _foreach_recurse(f, A::FVMArray, names::Tuple)
    s = first(names)
    f(_resolve(getfield(A, :data), getfield(getfield(A, :axes), s)))
    _foreach_recurse(f, A, Base.tail(names))
end

function foreach_field_recursive(f, A::FVMArray)
    _foreach_recurse(f, A, keys(A))
end

function loop_recursive!(u, n_cells)
    foreach_field_recursive(u.mass_fractions) do field_view
        for cell_id in 1:n_cells
            field_view[cell_id] += 1.0
        end
    end
end

# ============================================================================
# Run all benchmarks
# ============================================================================

println("=" ^ 60)
println("Key Looping Allocation Experiments (n_cells=$n_cells, 3 species)")
println("=" ^ 60)

results = Dict{String,Int}()
results["0-naive"]     = bench("Approach 0 (naive for loop)      ", loop_naive!, u, n_cells)
results["1-map"]       = bench("Approach 1 (map over keys)       ", loop_map!, u, n_cells)
results["2-foreach"]   = bench("Approach 2 (foreach over keys)   ", loop_foreach!, u, n_cells)
results["3-generated"] = bench("Approach 3 (@generated foreach)  ", loop_generated!, u, n_cells)
results["3b-gen-named"]= bench("Approach 3b(@generated named)    ", loop_generated_named!, u, n_cells)
results["4-tuple"]     = bench("Approach 4 (tuple of views)      ", loop_tuple!, u, n_cells)
results["5-macro"]     = bench("Approach 5 (@fvm_foreach macro)  ", loop_macro!, u, n_cells)
results["5b-macro-nm"] = bench("Approach 5b(@fvm_foreach_named)  ", loop_macro_named!, u, n_cells)
results["6-recursive"] = bench("Approach 6 (recursive head/tail) ", loop_recursive!, u, n_cells)

println("\n" * "=" ^ 60)
winners = [k for (k,v) in results if v == 0]
if !isempty(winners)
    println("Zero-allocation winners: ", join(winners, ", "))
else
    println("No zero-allocation approach found :(")
end
println("=" ^ 60)
