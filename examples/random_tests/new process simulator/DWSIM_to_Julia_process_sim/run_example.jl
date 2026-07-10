"""
    run_example.jl

Example driver that demonstrates the full DWSIM → Julia pipeline:

  1. Runs `extract_dwsim_flowsheet.py` on a `.dwxmz` file (producing JSON).
  2. Calls `load_dwsim_flowsheet(...)` to get Julia data structures.
  3. Inspects the resulting `topo`, `streams`, `components`, and compound data.

Prerequisites:
  - Python with `pythonnet` / DWSIM installed (for step 1)
  - JSON3.jl, ComponentArrays.jl (for step 2)
"""

using Revise

# Include the loader (adjust path if needed)
Revise.includet(joinpath(@__DIR__, "load_dwsim_flowsheet.jl"))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Step 1 — Run the Python extractor                                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# Path to the DWSIM flowsheet
dwxmz_path = joinpath(@__DIR__, "..", "small_recycle_test.dwxmz")
json_output = joinpath(@__DIR__, "flowsheet_data.json")

# Path to the Python extraction script
extractor_script = joinpath(@__DIR__, "extract_dwsim_flowsheet.py")

println("═" ^ 72)
println("  Step 1: Running Python extractor")
println("═" ^ 72)

# Run the Python script. Adjust the python executable if needed.
# If you have `pythonnet` installed in a specific env, set PYTHON_PATH accordingly.
python_exe = "python"

try
    run(`$python_exe $extractor_script $dwxmz_path $json_output`)
    println("✓ JSON written to: $json_output")
catch e
    @warn "Python extraction failed. If you've already generated the JSON manually, " *
          "you can skip this step and proceed to Step 2." exception = e
end

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Step 2 — Load into Julia                                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

println("\n" * "═" ^ 72)
println("  Step 2: Loading JSON into Julia data structures")
println("═" ^ 72)

result = load_dwsim_flowsheet(json_output)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Step 3 — Inspect results                                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

println("\n" * "═" ^ 72)
println("  Step 3: Inspection")
println("═" ^ 72)

println("\n── Property Packages ────────────────────────────────────────────────")
for pp in result.property_packages
    println("  • $pp")
end

println("\n── Species ──────────────────────────────────────────────────────────")
for sp in result.species_symbols
    cd = result.compound_data[sp]
    println("  • $sp  (MW = $(cd.molecular_weight) g/mol)")
end

println("\n── Material Stream Names ────────────────────────────────────────────")
for sn in result.stream_names
    println("  • $sn")
end

println("\n── Energy Stream Names ──────────────────────────────────────────────")
for en in result.energy_names
    println("  • $en")
end

println("\n── Component (Unit Op) Names ────────────────────────────────────────")
for cn in result.component_names
    println("  • $cn")
end

println("\n── Topo ─────────────────────────────────────────────────────────────")
for (name, comp) in pairs(result.topo)
    println("  $name => $(typeof(comp))")
    println("    $(comp)")
end

println("\n── Sample Stream Data ───────────────────────────────────────────────")
if !isempty(result.stream_names)
    first_stream = first(result.stream_names)
    s = getproperty(result.streams, first_stream)
    println("  Stream :$first_stream")
    println("    mass_flow     = $(s.mass_flow) kg/s")
    println("    temp          = $(s.temp) K")
    println("    pressure      = $(s.pressure) Pa")
    println("    vapor_fraction = $(s.vapor_fraction)")
    println("    mass_fractions = $(s.mass_fractions)")
end

println("\n── Components (Unit Op Properties) ─────────────────────────────────")
for cn in result.component_names
    println("  :$cn => $(getproperty(result.components, cn))")
end

println("\n✓ Done! All data structures are ready for the process simulator.")
