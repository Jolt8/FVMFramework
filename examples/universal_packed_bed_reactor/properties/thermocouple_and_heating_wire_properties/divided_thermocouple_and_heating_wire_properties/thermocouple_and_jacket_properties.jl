
function get_cell_lengths_along_property_group(grid, cell_set_name, cell_lengths_along_pipe)
    return cell_lengths_along_pipe[grid.cellsets[cell_set_name]]
end 

function get_thermocouple_closest_cell_id(grid, cell_lengths_along_pipe, axial_position_of_thermocouple)
    cell_ids = grid.cellsets["thermocouple_and_heating_wire"]
    return argmin(id -> abs(cell_lengths_along_pipe[id] - axial_position_of_thermocouple), cell_ids)
end


function get_thermocouple_and_jacket_properties(grid, pipe_length, n_cells_axial, cell_lengths_along_pipe)
    # ── Thermocouple axial positions (measured from inlet face) ──────────────
    # ⚠️ TODO: confirm exact axial positions of each TC bead from your physical build notes

    #note that cell_lengths_along_pipe contains the cell lenghts along the reactor for all 400 cells
    #thus, we have to index it with only the cells in this region to get the correct cell id
    TC1_closest_cell_id = get_thermocouple_closest_cell_id(grid, cell_lengths_along_pipe, 1.5u"inch" |> u"m")
    TC2_closest_cell_id = get_thermocouple_closest_cell_id(grid, cell_lengths_along_pipe, 3.0u"inch" |> u"m")
    TC3_closest_cell_id = get_thermocouple_closest_cell_id(grid, cell_lengths_along_pipe, 5.0u"inch" |> u"m")
    TC4_closest_cell_id = get_thermocouple_closest_cell_id(grid, cell_lengths_along_pipe, 7.5u"inch" |> u"m")
    TC5_closest_cell_id = get_thermocouple_closest_cell_id(grid, cell_lengths_along_pipe, 10.0u"inch" |> u"m")
    #note that cell_lenghts_along_pipe[region_cells] doesn't work here and returns the wrong index

    # ── Effective lumped properties for the TC+tape annular region ───────────
    # Physical layers in this annular shell (from pipe OD outward):
    #   inner mica tape:   ~1–2 mm,  k ≈ 0.50 W/(m·K),  ρ ≈ 2900 kg/m³,  cp ≈ 880 J/(kg·K)
    #   inner fiberglass:  ~3–5 mm,  k ≈ 0.07 W/(m·K),  ρ ≈ 200  kg/m³,  cp ≈ 840 J/(kg·K)
    #   Type-K TC bead:    negligible thermal mass
    #   outer fiberglass:  ~3–5 mm,  k ≈ 0.07 W/(m·K),  ρ ≈ 200  kg/m³,  cp ≈ 840 J/(kg·K)
    #   outer mica tape:   ~1–2 mm,  k ≈ 0.50 W/(m·K),  ρ ≈ 2900 kg/m³,  cp ≈ 880 J/(kg·K)
    #   mica wrap:         ~1–2 mm,  k ≈ 0.50 W/(m·K),  ρ ≈ 2900 kg/m³,  cp ≈ 880 J/(kg·K)
    #   heating wire:      thin nichrome/kanthal wire — heat source only, negligible mass
    #
    # The effective k below is a harmonic-mean (series resistance) dominated
    # by the low-k fiberglass layers. For a first pass we use the fiberglass k
    # since it is the dominant resistor; ρ and cp are volume-fraction weighted.
    #
    # Approximate total annular thickness: ~10–14 mm
    # ⚠️ TODO: confirm exact tape thicknesses from physical build

    # Volume fractions (rough, adjust if you measure thicknesses more precisely):
    # mica total ~4 mm, fiberglass total ~8 mm, out of ~12 mm total
    vf_mica      = 4.0 / 12.0  # ~1/3
    vf_fiberglass = 8.0 / 12.0  # ~2/3

    rho_mica         = 2900.0   # kg/m³
    cp_mica          = 880.0    # J/(kg·K)
    rho_fiberglass   = 200.0    # kg/m³  (loose woven tape, not dense board)
    cp_fiberglass    = 840.0    # J/(kg·K)
    k_fiberglass     = 0.07     # W/(m·K) — dominates series resistance

    rho_eff = vf_mica * rho_mica + vf_fiberglass * rho_fiberglass  # ≈ 1097 kg/m³
    cp_eff  = (vf_mica * rho_mica * cp_mica + vf_fiberglass * rho_fiberglass * cp_fiberglass) / rho_eff # ≈ 873 J/(kg·K)
    k_eff   = k_fiberglass  # series-resistance dominated by fiberglass

    thermocouple_and_heating_wire_properties = ComponentVector(
        k = k_eff * u"W/(m*K)",
        cp = cp_eff * u"J/(kg*K)",
        rho = rho_eff * u"kg/m^3",

        TC1_closest_cell_id = TC1_closest_cell_id,
        TC2_closest_cell_id = TC2_closest_cell_id,
        TC3_closest_cell_id = TC3_closest_cell_id,
        TC4_closest_cell_id = TC4_closest_cell_id,
        TC5_closest_cell_id = TC5_closest_cell_id,
    )

    return thermocouple_and_heating_wire_properties
end