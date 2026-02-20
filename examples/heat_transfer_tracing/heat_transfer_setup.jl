using ComponentArrays
using Unitful

# ── Initial Conditions (State Variables) ─────────────────────────────

# Region: copper (24 cells)
state_copper = ComponentVector(
    temp = 300.0
)

# Region: steel (24 cells)
state_steel = ComponentVector(
    temp = 350.0
)


# ── Fixed Properties ─────────────────────────────────────────────────

# Region: copper (24 cells)
fixed_copper = ComponentVector(
    k = 401.0,
    rho = 8960.0,
    cp = 385.0
)

# Region: steel (24 cells)
fixed_steel = ComponentVector(
    k = 15.0,
    rho = 7850.0,
    cp = 450.0
)


# ── Do not modify below this line ────────────────────────────────────

initial_state = Dict{Symbol, Any}(
    :copper => state_copper,
    :steel => state_steel
)

fixed_properties = Dict{Symbol, Any}(
    :copper => fixed_copper,
    :steel => fixed_steel
)
