"""
    load_dwsim_flowsheet.jl

Reads the JSON produced by `extract_dwsim_flowsheet.py` and reconstructs
the Julia data structures expected by the process simulator:

  - `species_symbols`  — Tuple of Symbols for compounds
  - `stream_names`     — Tuple of Symbols for material streams
  - `energy_names`     — Tuple of Symbols for energy streams
  - `component_names`  — Tuple of Symbols for unit operations
  - `streams`          — ComponentVector of all material + energy streams
  - `components`       — ComponentVector of all unit-op properties
  - `topo`             — NamedTuple mapping unit-op name => Julia struct
  - `property_packages`— Vector of property-package name strings
  - `compound_data`    — Dict{Symbol, NamedTuple} of per-species constant properties
"""

using JSON3
using ComponentArrays
using Unitful

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1 — Abstract types & unit-operation structs                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

abstract type AbstractComponent{CalcMode} end

# ── Existing types (matching small_recycle_test.jl) ────────────────────────────

struct Compressor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::Compressor) = (c.stream_in,)
get_inlet_energy_streams(c::Compressor) = (c.energy_in,)
get_outlet_material_streams(c::Compressor) = (c.stream_out,)
get_outlet_energy_streams(::Compressor) = ()

struct TwoMixer{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in_1::Symbol
    stream_in_2::Symbol
    stream_out::Symbol
end

get_inlet_material_streams(c::TwoMixer) = (c.stream_in_1, c.stream_in_2)
get_inlet_energy_streams(::TwoMixer) = ()
get_outlet_material_streams(c::TwoMixer) = (c.stream_out,)
get_outlet_energy_streams(::TwoMixer) = ()

struct ThreeMixer{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in_1::Symbol
    stream_in_2::Symbol
    stream_in_3::Symbol
    stream_out::Symbol
end

get_inlet_material_streams(c::ThreeMixer) = (c.stream_in_1, c.stream_in_2, c.stream_in_3)
get_inlet_energy_streams(::ThreeMixer) = ()
get_outlet_material_streams(c::ThreeMixer) = (c.stream_out,)
get_outlet_energy_streams(::ThreeMixer) = ()

struct TwoSplitter{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out_1::Symbol
    stream_out_2::Symbol
end

get_inlet_material_streams(c::TwoSplitter) = (c.stream_in,)
get_inlet_energy_streams(::TwoSplitter) = ()
get_outlet_material_streams(c::TwoSplitter) = (c.stream_out_1, c.stream_out_2)
get_outlet_energy_streams(::TwoSplitter) = ()

struct ThreeSplitter{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out_1::Symbol
    stream_out_2::Symbol
    stream_out_3::Symbol
end

get_inlet_material_streams(c::ThreeSplitter) = (c.stream_in,)
get_inlet_energy_streams(::ThreeSplitter) = ()
get_outlet_material_streams(c::ThreeSplitter) = (c.stream_out_1, c.stream_out_2, c.stream_out_3)
get_outlet_energy_streams(::ThreeSplitter) = ()

# ── New DWSIM unit-op types ───────────────────────────────────────────────────

struct Heater{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::Heater) = (c.stream_in,)
get_inlet_energy_streams(c::Heater) = (c.energy_in,)
get_outlet_material_streams(c::Heater) = (c.stream_out,)
get_outlet_energy_streams(::Heater) = ()

struct Cooler{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_out::Symbol
end

get_inlet_material_streams(c::Cooler) = (c.stream_in,)
get_inlet_energy_streams(::Cooler) = ()
get_outlet_material_streams(c::Cooler) = (c.stream_out,)
get_outlet_energy_streams(c::Cooler) = (c.energy_out,)

struct Pump{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::Pump) = (c.stream_in,)
get_inlet_energy_streams(c::Pump) = (c.energy_in,)
get_outlet_material_streams(c::Pump) = (c.stream_out,)
get_outlet_energy_streams(::Pump) = ()

struct Valve{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
end

get_inlet_material_streams(c::Valve) = (c.stream_in,)
get_inlet_energy_streams(::Valve) = ()
get_outlet_material_streams(c::Valve) = (c.stream_out,)
get_outlet_energy_streams(::Valve) = ()

struct Expander{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_out::Symbol
end

get_inlet_material_streams(c::Expander) = (c.stream_in,)
get_inlet_energy_streams(::Expander) = ()
get_outlet_material_streams(c::Expander) = (c.stream_out,)
get_outlet_energy_streams(c::Expander) = (c.energy_out,)

struct HeatExchanger{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    hot_stream_in::Symbol
    hot_stream_out::Symbol
    cold_stream_in::Symbol
    cold_stream_out::Symbol
end

get_inlet_material_streams(c::HeatExchanger) = (c.hot_stream_in, c.cold_stream_in)
get_inlet_energy_streams(::HeatExchanger) = ()
get_outlet_material_streams(c::HeatExchanger) = (c.hot_stream_out, c.cold_stream_out)
get_outlet_energy_streams(::HeatExchanger) = ()

struct ConversionReactor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::ConversionReactor) = (c.stream_in,)
get_inlet_energy_streams(c::ConversionReactor) = (c.energy_in,)
get_outlet_material_streams(c::ConversionReactor) = (c.stream_out,)
get_outlet_energy_streams(::ConversionReactor) = ()

struct EquilibriumReactor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::EquilibriumReactor) = (c.stream_in,)
get_inlet_energy_streams(c::EquilibriumReactor) = (c.energy_in,)
get_outlet_material_streams(c::EquilibriumReactor) = (c.stream_out,)
get_outlet_energy_streams(::EquilibriumReactor) = ()

struct GibbsReactor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::GibbsReactor) = (c.stream_in,)
get_inlet_energy_streams(c::GibbsReactor) = (c.energy_in,)
get_outlet_material_streams(c::GibbsReactor) = (c.stream_out,)
get_outlet_energy_streams(::GibbsReactor) = ()

struct CSTR{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::CSTR) = (c.stream_in,)
get_inlet_energy_streams(c::CSTR) = (c.energy_in,)
get_outlet_material_streams(c::CSTR) = (c.stream_out,)
get_outlet_energy_streams(::CSTR) = ()

struct PFR{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::PFR) = (c.stream_in,)
get_inlet_energy_streams(c::PFR) = (c.energy_in,)
get_outlet_material_streams(c::PFR) = (c.stream_out,)
get_outlet_energy_streams(::PFR) = ()

struct ShortcutColumn{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    distillate_out::Symbol
    bottoms_out::Symbol
    condenser_energy::Symbol
    reboiler_energy::Symbol
end

get_inlet_material_streams(c::ShortcutColumn) = (c.stream_in,)
get_inlet_energy_streams(::ShortcutColumn) = ()
get_outlet_material_streams(c::ShortcutColumn) = (c.distillate_out, c.bottoms_out)
get_outlet_energy_streams(c::ShortcutColumn) = (c.condenser_energy, c.reboiler_energy)

struct ComponentSeparator{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out_1::Symbol
    stream_out_2::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::ComponentSeparator) = (c.stream_in,)
get_inlet_energy_streams(c::ComponentSeparator) = (c.energy_in,)
get_outlet_material_streams(c::ComponentSeparator) = (c.stream_out_1, c.stream_out_2)
get_outlet_energy_streams(::ComponentSeparator) = ()

struct Recycle{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
end

get_inlet_material_streams(c::Recycle) = (c.stream_in,)
get_inlet_energy_streams(::Recycle) = ()
get_outlet_material_streams(c::Recycle) = (c.stream_out,)
get_outlet_energy_streams(::Recycle) = ()

# Catch-all for any DWSIM type we haven't explicitly modeled
struct GenericUnitOp{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    material_streams_in::Tuple{Vararg{Symbol}}
    material_streams_out::Tuple{Vararg{Symbol}}
    energy_streams_in::Tuple{Vararg{Symbol}}
    energy_streams_out::Tuple{Vararg{Symbol}}
end

get_inlet_material_streams(c::GenericUnitOp) = c.material_streams_in
get_inlet_energy_streams(c::GenericUnitOp) = c.energy_streams_in
get_outlet_material_streams(c::GenericUnitOp) = c.material_streams_out
get_outlet_energy_streams(c::GenericUnitOp) = c.energy_streams_out

# ── Stream graph-node types (for MetaGraphsNext) ──────────────────────────────

abstract type AbstractStream end

struct MaterialStreamNode <: AbstractStream
    id::Symbol
end

struct EnergyStreamNode <: AbstractStream
    id::Symbol
end

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2 — Name sanitization & property-name overrides                   ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
    sanitize_name(s::AbstractString) -> Symbol

Convert a free-form DWSIM tag (e.g. `"Compressor 1"`, `"My Feed Stream"`) into
a Julia-safe Symbol (`:compressor_1`, `:my_feed_stream`).
"""
function sanitize_name(s::AbstractString)::Symbol
    cleaned = lowercase(strip(s))
    cleaned = replace(cleaned, r"[^a-z0-9]+" => "_")
    cleaned = strip(cleaned, '_')
    # Ensure it doesn't start with a digit
    if !isempty(cleaned) && isdigit(first(cleaned))
        cleaned = "_" * cleaned
    end
    return Symbol(cleaned)
end

"""
    PROPERTY_NAME_OVERRIDES

User-editable dictionary that maps cryptic DWSIM property keys
(e.g. `"MSTR-001"`) to cleaner Symbol names (e.g. `:mass_flow`).

Only applies to **unit-operation property keys**. Stream tags and unit-op tags
are always sanitized via `sanitize_name`.

Add your own overrides here as needed.
"""
const PROPERTY_NAME_OVERRIDES = Dict{String, Symbol}(
    # DWSIM compressor properties
    "PROP_CO_0" => :pressure_out,
    "PROP_CO_1" => :delta_t,
    "PROP_CO_2" => :power_required,
    "PROP_CO_3" => :delta_p,
    "PROP_CO_4" => :efficiency_adiabatic,
    "PROP_CO_5" => :efficiency_polytropic,

    # DWSIM expander / turbine properties
    "PROP_EX_0" => :pressure_out,
    "PROP_EX_1" => :delta_t,
    "PROP_EX_2" => :power_generated,
    "PROP_EX_3" => :delta_p,
    "PROP_EX_4" => :efficiency_adiabatic,
    "PROP_EX_5" => :efficiency_polytropic,

    # DWSIM heater properties
    "PROP_HT_0" => :outlet_temperature,
    "PROP_HT_1" => :delta_t,
    "PROP_HT_2" => :heat_added,
    "PROP_HT_3" => :delta_p,
    "PROP_HT_4" => :efficiency,

    # DWSIM cooler properties
    "PROP_CL_0" => :outlet_temperature,
    "PROP_CL_1" => :delta_t,
    "PROP_CL_2" => :heat_removed,
    "PROP_CL_3" => :delta_p,
    "PROP_CL_4" => :efficiency,

    # DWSIM pump properties
    "PROP_PU_0" => :delta_p,
    "PROP_PU_1" => :delta_t,
    "PROP_PU_2" => :power_required,
    "PROP_PU_3" => :efficiency,
    "PROP_PU_4" => :npsh_available,

    # DWSIM valve properties
    "PROP_VA_0" => :pressure_out,
    "PROP_VA_1" => :delta_t,
    "PROP_VA_2" => :delta_p,

    # DWSIM heat exchanger properties
    "PROP_HX_0" => :overall_htc,
    "PROP_HX_1" => :area,
    "PROP_HX_2" => :heat_duty,
    "PROP_HX_3" => :hot_side_delta_t,
    "PROP_HX_4" => :cold_side_delta_t,
    "PROP_HX_5" => :hot_side_delta_p,
    "PROP_HX_6" => :cold_side_delta_p,

    # DWSIM energy stream
    "PROP_ES_0" => :energy_flow,

    # DWSIM mixer
    "PROP_MX_0" => :pressure_calculation,

    # DWSIM splitter
    #"PROP_SP_0" => :split_ratio_stream_1,
    #"PROP_SP_1" => :split_ratio_stream_2,
    #"PROP_SP_2" => :split_ratio_stream_3,
    "SR1" => :split_ratio_stream_1,
    "SR2" => :split_ratio_stream_2,
    "SR3" => :split_ratio_stream_3
)

"""
    resolve_property_name(raw_key::String) -> Symbol

Look up `raw_key` in PROPERTY_NAME_OVERRIDES; if absent, fall through to
`sanitize_name`.
"""
function resolve_property_name(raw_key::String)::Symbol
    return get(PROPERTY_NAME_OVERRIDES, raw_key, sanitize_name(raw_key))
end


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3 — DWSIM type-name → Julia struct mapping                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
Map from DWSIM `GetType().Name` string to the Julia struct constructor.
Energy-port conventions:
  - Compressor / Heater / Pump / reactors → energy_in
  - Cooler / Expander                     → energy_out
  - Valve                                 → no energy port
  - HeatExchanger                         → two material in, two material out
  - Mixer / Splitter                      → N material in or out (2 or 3)
"""
const DWSIM_TYPE_MAP = Dict{String, Type}(
    "Compressor"          => Compressor,
    "Mixer"               => TwoMixer,       # may upgrade to ThreeMixer dynamically
    "Splitter"            => TwoSplitter,    # may upgrade to ThreeSplitter dynamically
    "Heater"              => Heater,
    "Cooler"              => Cooler,
    "Pump"                => Pump,
    "Valve"               => Valve,
    "Expander"            => Expander,
    "HeatExchanger"       => HeatExchanger,
    "ConversionReactor"   => ConversionReactor,
    "EquilibriumReactor"  => EquilibriumReactor,
    "GibbsReactor"        => GibbsReactor,
    "Reactor_CSTR"        => CSTR,
    "Reactor_PFR"         => PFR,
    "ShortcutColumn"      => ShortcutColumn,
    "ComponentSeparator"  => ComponentSeparator,
    "Recycle"             => Recycle,
)


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4 — Build topo from edges                                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
    classify_streams(data)

Given the parsed JSON, return two Sets of sanitized-Symbol stream tags:
one for material streams, one for energy streams.
"""
function classify_streams(data)
    material_tags = Set{Symbol}()
    energy_tags = Set{Symbol}()

    for material_stream in data["material_streams"]
        push!(material_tags, sanitize_name(material_stream["tag"]))
    end
    for energy_stream in data["energy_streams"]
        push!(energy_tags, sanitize_name(energy_stream["tag"]))
    end

    return material_tags, energy_tags
end

"""
    build_connectivity(data, material_tags, energy_tags)

From the JSON edges, build per-unit-op connectivity:
    unit_op_sym => (mat_in = Symbol[], mat_out = Symbol[],
                    ene_in = Symbol[], ene_out = Symbol[])
"""
function build_connectivity(data, material_tags, energy_tags)
    # Set of unit-op sanitized tags
    unit_op_tags = Set{Symbol}(sanitize_name(unit_op["tag"]) for unit_op in data["unit_operations"])

    # Adjacency
    connectivity = Dict{Symbol, NamedTuple{(:mat_in, :mat_out, :ene_in, :ene_out),
                                            NTuple{4, Vector{Symbol}}}}()

    for unit_op in data["unit_operations"]
        sym = sanitize_name(unit_op["tag"])
        connectivity[sym] = (mat_in = Symbol[], mat_out = Symbol[],
                             ene_in = Symbol[], ene_out = Symbol[])
    end

    for edge in data["edges"]
        source = sanitize_name(edge["source"])
        target = sanitize_name(edge["target"])

        # Edge from stream → unit_op  ⟹  stream is an inlet to that unit op
        if (source ∈ material_tags) && (target ∈ unit_op_tags)
            push!(connectivity[target].mat_in, source)
        elseif (source ∈ energy_tags) && (target ∈ unit_op_tags)
            push!(connectivity[target].ene_in, source)

        # Edge from unit_op → stream  ⟹  stream is an outlet from that unit op
        elseif (source ∈ unit_op_tags) && (target ∈ material_tags)
            push!(connectivity[source].mat_out, target)
        elseif (source ∈ unit_op_tags) && (target ∈ energy_tags)
            push!(connectivity[source].ene_out, target)
        end
    end

    return connectivity
end

"""
    make_component_struct(type_name, calc_mode_str, id, connection)

Instantiate the correct Julia struct for this unit operation.
"""
function make_component_struct(type_name::String, calc_mode_str, id::Symbol, connection)
    if calc_mode_str === nothing
        calc_mode = :Default
    else
        calc_mode = Symbol(calc_mode_str)
    end

    base_type = get(DWSIM_TYPE_MAP, type_name, nothing)

    # ── Mixer: pick TwoMixer or ThreeMixer based on actual inlet count ─────
    if type_name == "Mixer"
        streams_in = length(connection.mat_in)
        if streams_in == 3
            return ThreeMixer{calc_mode}(id,
                connection.mat_in[1], connection.mat_in[2], connection.mat_in[3],
                first(connection.mat_out))
        elseif streams_in >= 2
            return TwoMixer{calc_mode}(id,
                connection.mat_in[1], connection.mat_in[2],
                first(connection.mat_out))
        else
            # single-inlet mixer — treat as pass-through, use GenericUnitOp
            return GenericUnitOp{calc_mode}(id,
                Tuple(connection.mat_in), Tuple(connection.mat_out),
                Tuple(connection.ene_in), Tuple(connection.ene_out))
        end
    end

    # ── Splitter: pick TwoSplitter or ThreeSplitter ────────────────────────
    if type_name == "Splitter"
        streams_out = length(connection.mat_out)
        if streams_out == 3
            return ThreeSplitter{calc_mode}(id,
                first(connection.mat_in),
                connection.mat_out[1], connection.mat_out[2], connection.mat_out[3])
        elseif streams_out >= 2
            return TwoSplitter{calc_mode}(id,
                first(connection.mat_in),
                connection.mat_out[1], connection.mat_out[2])
        else
            return GenericUnitOp{calc_mode}(id,
                Tuple(connection.mat_in), Tuple(connection.mat_out),
                Tuple(connection.ene_in), Tuple(connection.ene_out))
        end
    end

    # ── HeatExchanger: 2 material in, 2 material out ──────────────────────
    if type_name == "HeatExchanger"
        return HeatExchanger{calc_mode}(id,
            connection.mat_in[1], connection.mat_out[1],
            connection.mat_in[2], connection.mat_out[2])
    end

    # ── ShortcutColumn: 1 mat in, 2 mat out, 2 energy out ─────────────────
    if type_name == "ShortcutColumn"
        return ShortcutColumn{calc_mode}(id,
            first(connection.mat_in),
            connection.mat_out[1], connection.mat_out[2],
            connection.ene_out[1], connection.ene_out[2])
    end

    # ── ComponentSeparator: 1 mat in, 2 mat out, 1 energy in ──────────────
    if type_name == "ComponentSeparator"
        if isempty(connection.ene_in)
            energy_in = :none
        else
            energy_in = first(connection.ene_in)
        end
        return ComponentSeparator{calc_mode}(id,
            first(connection.mat_in),
            connection.mat_out[1], connection.mat_out[2],
            energy_in)
    end

    if type_name == "Recycle"
        return Recycle{calc_mode}(id,
            first(connection.mat_in),
            connection.mat_out[1])
    end

    # ── Standard 1-in / 1-out with energy port ─────────────────────────────
    if base_type !== nothing
        # Types with energy_in: Compressor, Heater, Pump, reactors
        if base_type in (Compressor, Heater, Pump,
                         ConversionReactor, EquilibriumReactor, GibbsReactor,
                         CSTR, PFR)
            if isempty(connection.ene_in)
                energy_in = :none
            else
                energy_in = first(connection.ene_in)
            end
            return base_type{calc_mode}(id, first(connection.mat_in), first(connection.mat_out), energy_in)
        end

        # Types with energy_out: Cooler, Expander
        if base_type in (Cooler, Expander)
            if isempty(connection.ene_out)
                energy_out = :none
            else
                energy_out = first(connection.ene_out)
            end
            return base_type{calc_mode}(id, first(connection.mat_in), first(connection.mat_out), energy_out)
        end

        # No energy port: Valve
        if base_type === Valve
            return Valve{calc_mode}(id, first(connection.mat_in), first(connection.mat_out))
        end
    end

    # ── Fallback: GenericUnitOp ────────────────────────────────────────────
    @warn "Unknown DWSIM type '$type_name' for '$id' — wrapping in GenericUnitOp"
    return GenericUnitOp{calc_mode}(id,
        Tuple(connection.mat_in), Tuple(connection.mat_out),
        Tuple(connection.ene_in), Tuple(connection.ene_out))
end

"""
    build_topo(data, material_tags, energy_tags)

Return a NamedTuple of `AbstractComponent` structs, keyed by sanitized unit-op name.
"""
function build_topo(data, material_tags, energy_tags)
    connectivity = build_connectivity(data, material_tags, energy_tags)

    topo_pairs = Pair{Symbol, AbstractComponent}[]

    for unit_op in data["unit_operations"]
        id = sanitize_name(unit_op["tag"])
        type_name = unit_op["type_name"]
        calc_mode = get(unit_op, "calc_mode", nothing)

        connection = connectivity[id]
        comp = make_component_struct(type_name, calc_mode, id, connection)
        push!(topo_pairs, id => comp)
    end

    return (; (p.first => p.second for p in topo_pairs)...)
end


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5 — Build streams & components ComponentVectors                   ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
    build_compound_data(data) -> (species_symbols, compound_data)

Returns `species_symbols::Tuple{Vararg{Symbol}}` and a
`Dict{Symbol, NamedTuple}` with molecular_weight, ig_cp constants, etc.
"""
function build_compound_data(data)
    compounds = data["compounds"]
    species = Symbol[]
    compound_data = Dict{Symbol, NamedTuple}()

    for compound in compounds
        sym = sanitize_name(compound["name"])
        push!(species, sym)

        mw = compound["molecular_weight"]
        if mw === nothing
            mw = 0.0
        end

        tc = compound["critical_temperature"]
        if tc === nothing
            tc = 0.0
        end

        pc = compound["critical_pressure"]
        if pc === nothing
            pc = 0.0
        end

        af = compound["acentric_factor"]
        if af === nothing
            af = 0.0
        end

        ig_a = get(compound, "ig_cp_a", nothing)
        if ig_a === nothing
            ig_a = 0.0
        end

        ig_b = get(compound, "ig_cp_b", nothing)
        if ig_b === nothing
            ig_b = 0.0
        end

        ig_c = get(compound, "ig_cp_c", nothing)
        if ig_c === nothing
            ig_c = 0.0
        end

        ig_d = get(compound, "ig_cp_d", nothing)
        if ig_d === nothing
            ig_d = 0.0
        end

        ig_e = get(compound, "ig_cp_e", nothing)
        if ig_e === nothing
            ig_e = 0.0
        end

        compound_data[sym] = (
            molecular_weight = mw,
            critical_temperature = tc,
            critical_pressure = pc,
            acentric_factor = af,
            ig_cp_a = ig_a,
            ig_cp_b = ig_b,
            ig_cp_c = ig_c,
            ig_cp_d = ig_d,
            ig_cp_e = ig_e,
            ig_cp_equation = string(get(compound, "ig_cp_equation", "")),
        )
    end

    return Tuple(species), compound_data
end

"""
    make_material_stream_template(species_symbols)

Build a ComponentVector template for a material stream, matching the layout
of `small_recycle_test.jl` but parameterized on `species_symbols`.
Includes `vapor_fraction` as requested.
"""
function make_material_stream_template(species_symbols)
    zeros_unitless = NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])

    return ComponentVector(
        mass_flow = 0.0,            # kg/s  (base SI, unitless)
        mass_flows = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        molar_flow = 0.0,           # mol/s
        molar_flows = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        volumetric_flow = 0.0,      # m^3/s
        volumetric_flows = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        temp = 0.0,                 # K
        pressure = 0.0,             # Pa
        vapor_fraction = 0.0,       # molar vapor fraction (0–1)
        mass_fractions = (ComponentVector(; zeros_unitless...),),
        molar_fractions = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        molar_concentrations = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        specific_heat_ratio = 0.0,
        molecular_weights = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        mw_avg = 0.0,               # g/mol
        rho = 0.0,                  # kg/m^3
        species_cps = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        cp_avg = 0.0,               # J/(kg·K)
        species_cvs = (ComponentVector(;
            NamedTuple{species_symbols}([0.0 for _ in eachindex(species_symbols)])...),),
        cv_avg = 0.0,               # J/(kg·K)
        R_gas = 8.314,              # J/(mol·K)
    )
end

"""
    make_energy_stream_template()

Build a ComponentVector template for an energy stream.
"""
function make_energy_stream_template()
    return ComponentVector(
        electrical_wattage = 0.0,   # W
        heat_wattage = 0.0,         # W
    )
end

"""
    build_streams(data, species_symbols, compound_data)

Create and populate the combined `streams` ComponentVector.
"""
function build_streams(data, species_symbols, compound_data)
    material_template = make_material_stream_template(species_symbols)
    energy_template = make_energy_stream_template()

    # ── Collect names ──────────────────────────────────────────────────────
    stream_names = Tuple(sanitize_name(material_stream["tag"]) for material_stream in data["material_streams"])
    energy_names = Tuple(sanitize_name(energy_stream["tag"]) for energy_stream in data["energy_streams"])

    # ── Build the combined ComponentVector ─────────────────────────────────
    streams = ComponentVector(;
        NamedTuple{stream_names}([deepcopy(material_template) for _ in eachindex(stream_names)])...,
        NamedTuple{energy_names}([deepcopy(energy_template) for _ in eachindex(energy_names)])...,
    )

    # ── Populate material stream values ────────────────────────────────────
    for material_stream in data["material_streams"]
        sym = sanitize_name(material_stream["tag"])
        stream_ref = getproperty(streams, sym)

        if material_stream["mass_flow"] !== nothing
            stream_ref.mass_flow = material_stream["mass_flow"]
        end
        if material_stream["temperature"] !== nothing
            stream_ref.temp = material_stream["temperature"]
        end
        if material_stream["pressure"] !== nothing
            stream_ref.pressure = material_stream["pressure"]
        end
        if material_stream["vapor_fraction"] !== nothing
            stream_ref.vapor_fraction = material_stream["vapor_fraction"]
        end

        # Mass fractions
        if material_stream["mass_fractions"] !== nothing
            i = 1
            for (comp_name, frac) in pairs(material_stream["mass_fractions"])
                species_sym = sanitize_name(String(comp_name))
                if species_sym in species_symbols
                    stream_mass_fraction = getproperty(getproperty(stream_ref, :mass_fractions), i)
                    stream_mass_fraction = frac
                    i += 1
                end
            end
        end

        # Fill in compound constant properties (molecular_weights from database)
        i = 1
        for species_symbol in species_symbols
            if haskey(compound_data, species_symbol)
                this_compound_data = compound_data[species_symbol]
                stream_molecular_weight = getproperty(getproperty(stream_ref, :molecular_weights), i)
                stream_molecular_weight = this_compound_data.molecular_weight  # g/mol
                i += 1
            end
        end
    end

    # ── Populate energy stream values ──────────────────────────────────────
    for energy_stream in data["energy_streams"]
        sym = sanitize_name(energy_stream["tag"])
        stream_ref = getproperty(streams, sym)

        if energy_stream["power"] !== nothing
            stream_ref.electrical_wattage = energy_stream["power"]  # W
        end
    end

    return streams, stream_names, energy_names
end

"""
    build_components(data)

Create the `components` ComponentVector with all DWSIM properties dumped
as flat numeric values per unit op. Property names are resolved through
`PROPERTY_NAME_OVERRIDES` → `sanitize_name` fallback.
"""
function build_components(data)
    unit_ops = data["unit_operations"]

    if isempty(unit_ops)
        return ComponentVector(), ()
    end

    component_names = Tuple(sanitize_name(unit_op["tag"]) for unit_op in unit_ops)

    # Build per-unit-op NamedTuples of properties
    component_pairs = Pair{Symbol, Any}[]

    for unit_op in unit_ops
        sym = sanitize_name(unit_op["tag"])
        props = get(unit_op, "properties", nothing)

        if props === nothing || isempty(props)
            # Unit op with no extractable numeric properties → scalar placeholder
            push!(component_pairs, sym => 0.0)
            continue
        end

        # Collect numeric property pairs
        prop_pairs = Pair{Symbol, Float64}[]
        for (raw_key, pval) in pairs(props)
            pname = resolve_property_name(String(raw_key))
            val = pval["value"]
            if val !== nothing
                push!(prop_pairs, pname => Float64(val))
            end
        end

        if isempty(prop_pairs)
            push!(component_pairs, sym => 0.0)
        else
            cv = ComponentVector(; NamedTuple(prop_pairs)...)
            push!(component_pairs, sym => cv)
        end
    end

    components = ComponentVector(; NamedTuple(component_pairs)...)
    return components, component_names
end


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6 — Top-level entry point                                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
    load_dwsim_flowsheet(json_path::String)

Parse a JSON file produced by `extract_dwsim_flowsheet.py` and return
a NamedTuple with all the Julia objects the process simulator needs:

```julia
result = load_dwsim_flowsheet("flowsheet_data.json")
result.topo              # NamedTuple of AbstractComponent structs
result.streams           # ComponentVector
result.components        # ComponentVector
result.species_symbols   # Tuple of Symbols
result.stream_names      # Tuple of Symbols
result.energy_names      # Tuple of Symbols
result.component_names   # Tuple of Symbols
result.property_packages # Vector{String}
result.compound_data     # Dict{Symbol, NamedTuple}
```
"""
function load_dwsim_flowsheet(json_path::String)
    raw = read(json_path, String)
    data = JSON3.read(raw)

    # Property packages
    property_packages = String[string(pp) for pp in data["property_packages"]]

    # Compounds
    species_symbols, compound_data = build_compound_data(data)

    # Classify stream tags
    material_tags, energy_tags = classify_streams(data)

    # Streams
    streams, stream_names, energy_names = build_streams(data, species_symbols, compound_data)

    # Topo
    topo = build_topo(data, material_tags, energy_tags)

    # Components
    components, component_names = build_components(data)

    return (
        topo = topo,
        streams = streams,
        components = components,
        species_symbols = species_symbols,
        stream_names = stream_names,
        energy_names = energy_names,
        component_names = component_names,
        property_packages = property_packages,
        compound_data = compound_data,
    )
end
