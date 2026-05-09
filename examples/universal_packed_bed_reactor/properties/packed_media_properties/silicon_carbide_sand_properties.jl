function get_silicon_carbide_sand_properties(pipe_length, n_cells, cell_lengths_along_pipe)
    # ── Silicon Carbide (SiC) packing / preheater bed properties ─────────────
    # SiC is used here purely as a thermal pre-heater / evaporator bed,
    # not as a catalyst. No reactions occur in this zone.
    #
    # Bulk / packed-bed values (solid + void):
    #   SiC solid:  k ≈ 120 W/(m·K), ρ ≈ 3210 kg/m³, cp ≈ 750 J/(kg·K)  at ~300°C
    #   Bed void fraction ε ≈ 0.38–0.42 for random packing of equal spheres
    #
    # The effective bed k (Zehner-Schlünder or simple mixing rule) is dominated
    # by the gas phase in a packed bed, but for a hot SiC bed the solid conduction
    # is high enough that the effective k is much higher than pure air.
    # A conservative estimate using the simple parallel-model lower bound:
    #   k_eff ≈ (1 - ε) * k_solid + ε * k_gas ≈ 0.60 * 120 + 0.40 * 0.05 ≈ 72 W/(m·K)
    # However for thin beds and coarse particles a value of 1–5 W/(m·K) is
    # commonly used in literature for gas-filled packed beds. 
    # ⚠️ TODO: clarify particle size and packing density of the SiC bed

    bed_void_fraction = 0.40  # dimensionless

    # Effective packed-bed properties (volume-fraction weighted):
    rho_SiC = 3210.0   # kg/m³  solid SiC density
    cp_SiC  = 750.0    # J/(kg·K) at ~300°C; use 670 at room temp if preferred
    k_SiC   = 120.0    # W/(m·K) solid SiC

    rho_gas = 1.0      # kg/m³  (air/steam approximate)
    cp_gas  = 1050.0   # J/(kg·K)
    k_gas   = 0.05     # W/(m·K)

    rho_eff = (1.0 - bed_void_fraction) * rho_SiC + bed_void_fraction * rho_gas   # ≈ 1926 kg/m³
    cp_eff  = ((1.0 - bed_void_fraction) * rho_SiC * cp_SiC + bed_void_fraction * rho_gas * cp_gas) / rho_eff
    k_eff   = (1.0 - bed_void_fraction) * k_SiC + bed_void_fraction * k_gas  # parallel upper bound ≈ 72 W/(m·K)
    # ⚠️ TODO: k_eff is highly uncertain for packed beds — consider using 1–10 W/(m·K)
    # from literature if you don't have a better estimate

    particle_diameter = 1.0u"mm"  # ⚠️ TODO: confirm actual SiC particle/grain size

    silicon_carbide_sand_properties = ComponentVector(
        k = k_eff * u"W/(m*K)",
        cp = cp_eff * u"J/(kg*K)",
        rho = rho_eff * u"kg/m^3",

        #pipe_inside_diameter = 16.0u"mm",
        #pipe_length = pipe_length,
        #per_cell_pipe_length = pipe_length / n_cells,
        #cell_lengths_along_pipe = cell_lengths_along_pipe,

        bed_void_fraction = bed_void_fraction,
        particle_diameter = particle_diameter,

        diffusion_coefficients = (
            methanol = 1e-5u"m^2/s",
            water = 1e-5u"m^2/s",
            carbon_monoxide = 1e-5u"m^2/s",
            hydrogen = 1e-5u"m^2/s",
            carbon_dioxide = 1e-5u"m^2/s",
            air = 1e-5u"m^2/s"
        ),
        molecular_weights = (
            methanol = 32.04u"g/mol",
            water = 18.02u"g/mol",
            carbon_monoxide = 28.01u"g/mol",
            hydrogen = 2.02u"g/mol",
            carbon_dioxide = 44.01u"g/mol",
            air = 28.97u"g/mol"
        ),
    )

    return silicon_carbide_sand_properties
end
