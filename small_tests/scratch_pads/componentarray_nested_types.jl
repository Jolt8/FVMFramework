using ComponentArrays

u = (
    mass_fractions = (
        methanol = [0.0, 0.0],
        water = [0.0, 0.0]
    ),
    net_rates = (
        WGS_rxn = 0.0,
        MD_rxn = 0.0,
        MSR_rxn = 0.0
    ),
    temp = [0.0, 0.0],
    R_gas = 8.314
)

u = ComponentVector(u)


typeof(u.mass_fractions)
typeof(u.net_rates)
typeof(u.temp)
typeof(u.R_gas)

u.mass_fractions = ComponentArray(; u.mass_fractions..., test = 0.0)


test = Dict(
    Symbol("a") => 1
)

ComponentArray(test)