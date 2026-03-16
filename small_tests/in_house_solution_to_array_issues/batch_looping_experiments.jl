#=
    Batch-compatible key looping experiments.
    
    Problem: `map`/`foreach` closures + `@batch` = massive allocations.
    Manual unrolling + `@batch` = fast and ~0 allocs.
    Goal: find a way to auto-unroll the inner species loop so @batch works.
=#


#Commenting here: YOU DID IT!

using Polyester
include("FVMArray.jl")

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

# ---- Helpers from FVMArray ----
@inline function _resolve_val(A::FVMArray, ::Val{s}) where {s}
    return _resolve(getfield(A, :data), getfield(getfield(A, :axes), s))
end

# ============================================================================
# BASELINE: manual unroll + @batch (user's known-good reference)
# ============================================================================

function baseline_manual!(u, n_cells)
    @batch for cell_id in 1:n_cells
        u.mass_fractions[:methylene_blue][cell_id] += 1.0
        u.mass_fractions[:water][cell_id] += 1.0
        u.mass_fractions[:ethanol][cell_id] += 1.0
        u.molar_concentrations[:methylene_blue][cell_id] += 1.0
        u.molar_concentrations[:water][cell_id] += 1.0
        u.molar_concentrations[:ethanol][cell_id] += 1.0
    end
end

# ============================================================================
# BASELINE 2: map + @batch (user's known-BAD reference)
# ============================================================================

function baseline_map_batch!(u, n_cells)
    @batch for cell_id in 1:n_cells
        map(keys(u.mass_fractions)) do species
            u.mass_fractions[species][cell_id] += 1.0
            u.molar_concentrations[species][cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach A: Pre-compute views as tuples OUTSIDE @batch, iterate tuples inside
# ============================================================================

@generated function field_views(A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(_resolve_val(A, Val{$(QuoteNode(n))}())) for n in field_names]
    return quote
        Base.@_inline_meta
        ($(exprs...),)
    end
end

function approach_a_precompute!(u, n_cells)
    mf_views = field_views(u.mass_fractions)
    mc_views = field_views(u.molar_concentrations)
    @batch for cell_id in 1:n_cells
        for (mf, mc) in zip(mf_views, mc_views)
            mf[cell_id] += 1.0
            mc[cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach B: @generated inner function that takes cell_id and does all species
# ============================================================================

@generated function _do_all_species!(groupA::FVMArray{T,D,AxA}, groupB::FVMArray{T2,D2,AxB}, cell_id::Int) where {T,D,AxA,T2,D2,AxB}
    # Unroll over all fields in groupA (assumes groupA and groupB have same keys)
    field_names = fieldnames(AxA)
    exprs = []
    for n in field_names
        push!(exprs, quote
            _resolve_val(groupA, Val{$(QuoteNode(n))}())[cell_id] += 1.0
            _resolve_val(groupB, Val{$(QuoteNode(n))}())[cell_id] += 1.0
        end)
    end
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

function approach_b_generated_inner!(u, n_cells)
    @batch for cell_id in 1:n_cells
        _do_all_species!(u.mass_fractions, u.molar_concentrations, cell_id)
    end
end

# ============================================================================
# Approach C: @generated that takes a function + group + cell_id 
#   Most generic — user passes their own body function
# ============================================================================
#WINNER WINNER CHICKEN DINNER!!!

@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{FVMArray, N}) where {N}
    # All groups must have the same keys (uses first group's keys)
    first_ax = fieldnames(groups[1].parameters[3])  # Ax type parameter
    exprs = []
    for n in first_ax
        view_exprs = [:(  _resolve_val(groups[$i], Val{$(QuoteNode(n))}())  ) for i in 1:N]
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

function approach_c_generic!(u, n_cells)
    @batch for cell_id in 1:n_cells
        foreach_field_at!(cell_id, u.mass_fractions, u.molar_concentrations) do cid, mf, mc
            mf[cid] += 1.0
            mc[cid] += 1.0
        end
    end
end

# ============================================================================
# Approach D: Macro that unrolls the species loop at parse time
#   The macro expands BEFORE @batch sees the code, so @batch only sees
#   direct field accesses with no closures.
# ============================================================================

"""
    @unroll_fields group_expr var_name body

Expands at macro-expansion time into one copy of `body` per field in the group,
substituting `var_name` with each field's symbol literal — no closures involved.

NOTE: This requires knowing the field names at macro expansion time,
so we use a @generated approach under the hood instead.
"""

# Actually, the cleanest approach for @batch is a @generated function:
@generated function batch_foreach_field(f, A::FVMArray{T,D,Ax}) where {T,D,Ax}
    field_names = fieldnames(Ax)
    exprs = [:(f($(QuoteNode(n)), _resolve_val(A, Val{$(QuoteNode(n))}()))) for n in field_names]
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

function approach_d_batch_foreach!(u, n_cells)
    @batch for cell_id in 1:n_cells
        batch_foreach_field(u.mass_fractions) do name, mf_view
            mf_view[cell_id] += 1.0
        end
        batch_foreach_field(u.molar_concentrations) do name, mc_view
            mc_view[cell_id] += 1.0
        end
    end
end

# variant: also pass name as a Val for dispatch if needed
function approach_d_combined!(u, n_cells)
    @batch for cell_id in 1:n_cells
        batch_foreach_field(u.mass_fractions) do name, view
            view[cell_id] += 1.0
            # use `name` to index into a parallel group:
            u.molar_concentrations[name][cell_id] += 1.0
        end
    end
end

# ============================================================================
# Approach E: ntuple-based — get views as NTuple, loop the NTuple inside @batch
# ============================================================================

function approach_e_ntuple!(u, n_cells)
    mf = field_views(u.mass_fractions)
    mc = field_views(u.molar_concentrations)
    n_fields = length(mf)
    @batch for cell_id in 1:n_cells
        for i in 1:n_fields
            mf[i][cell_id] += 1.0
            mc[i][cell_id] += 1.0
        end
    end
end

# ============================================================================
# Run benchmarks 
# ============================================================================

function bench(label, f!, u, n; warmup=3, runs=50)
    for _ in 1:warmup; f!(u, n); end
    allocs = @allocated f!(u, n)
    total = 0.0
    for _ in 1:runs; total += @elapsed f!(u, n); end
    avg_us = total / runs * 1e6
    status = allocs <= 128 ? "✓ $allocs bytes" : "✗ $allocs bytes"
    println("  $label: $status, avg $(round(avg_us, digits=2)) μs")
    return allocs
end

println("=" ^ 70)
println("@batch-Compatible Key Looping (n_cells=$n_cells, 3 species × 2 groups)")
println("=" ^ 70)

results = Dict{String,Int}()
results["baseline-manual"]  = bench("BASELINE manual unroll         ", baseline_manual!, u, n_cells)
results["baseline-map-BAD"] = bench("BASELINE map+@batch (BAD)      ", baseline_map_batch!, u, n_cells)
results["A-precompute"]     = bench("A: pre-compute views tuple     ", approach_a_precompute!, u, n_cells)
results["B-generated-inner"]= bench("B: @generated inner function   ", approach_b_generated_inner!, u, n_cells)
results["C-generic"]        = bench("C: foreach_field_at! (generic) ", approach_c_generic!, u, n_cells)
results["D-batch-foreach"]  = bench("D: batch_foreach_field         ", approach_d_batch_foreach!, u, n_cells)
results["D-combined"]       = bench("D: batch_foreach combined      ", approach_d_combined!, u, n_cells)
results["E-ntuple"]         = bench("E: ntuple views + index loop   ", approach_e_ntuple!, u, n_cells)

println("\n" * "=" ^ 70)
good = [(k,v) for (k,v) in results if v <= 128]  # 128 bytes = @batch's own allocation
sort!(good, by=x->x[2])
println("Low-allocation approaches (≤128 bytes = @batch overhead):")
for (k,v) in good
    println("  $k: $v bytes")
end
println("=" ^ 70)

@btime approach_a_precompute!($u, $n_cells) #1.570 μs (0 allocations: 0 bytes)
@btime approach_b_generated_inner!($u, $n_cells) #1.780 μs (1 allocation: 128 bytes)
@btime approach_c_generic!($u, $n_cells) #1.970 μs (1 allocation: 128 bytes)
@btime approach_d_batch_foreach!($u, $n_cells) #6.875 μs (0 allocations: 0 bytes)
@btime approach_d_combined!($u, $n_cells) #788.500 μs (24936 allocations: 1.28 MiB)
@btime approach_e_ntuple!($u, $n_cells) #20.800 μs (0 allocations: 0 bytes)
