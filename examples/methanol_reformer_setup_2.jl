using ComponentArrays
using Unitful

# ── Initial Conditions (State Variables) ─────────────────────────────

# Region: reforming_area (1424 cells)
state_reforming_area = ComponentVector(
        mass_fractions = ComponentVector(
        water = 1.0,
        methanol = 1.3,
        carbon_dioxide = 0.0001,
        carbon_monoxide = 0.0001,
        hydrogen = 0.0001
    ),
    pressure = ustrip(1.0u"bar" |> u"Pa"),
    temp = ustrip(250.0u"°C" |> u"K")
)

# Region: inlet (24 cells)
state_inlet = ComponentVector(
        mass_fractions = ComponentVector(
        water = 0.0,
        methanol = 0.0,
        carbon_dioxide = 0.0,
        carbon_monoxide = 0.0,
        hydrogen = 0.0
    ),
    pressure = ustrip(1.1u"bar" |> u"Pa"),
    temp = ustrip(250.0u"°C" |> u"K")
)

# Region: outlet (24 cells)
state_outlet = ComponentVector(
        mass_fractions = ComponentVector(
        methanol = 0.0,
        water = 0.0,
        carbon_dioxide = 0.0,
        carbon_monoxide = 0.0,
        hydrogen = 0.0
    ),
    pressure = ustrip(1.1u"bar" |> u"Pa"),
    temp = ustrip(250.0u"°C" |> u"K")
)

# Region: wall (6382 cells)
state_wall = ComponentVector(
    temp = ustrip(250.0u"°C" |> u"K")
)

# Region: heating_areas (178 cells)
state_heating_areas = ComponentVector(
    temp = ustrip(250.0u"°C" |> u"K")
)

# Controller: temp_controller
state_temp_controller = ComponentVector(
    integral_error = 0.0
)


# ── Fixed Properties ─────────────────────────────────────────────────

# Region: reforming_area (1424 cells)
fixed_reforming_area = ComponentVector(
    permeability = 6.0e-12,
    mu = 1.0e-5,
    k = 237.0,
    cp = 3000.0
    #=species_diffusion_coeffs = ComponentVector(
        water = 1e-5,
        methanol = 1e-5,
        carbon_dioxide = 1e-5,
        carbon_monoxide = 1e-5,
        hydrogen = 1e-5
    )=#
)

# Region: inlet (24 cells)
fixed_inlet = ComponentVector(
    permeability = 6.0e-12,
    mu = 1.0e-5,
    k = 237.0,
    cp = 3000.0
    #=species_diffusion_coeffs = ComponentVector(
        water = 1e-5,
        methanol = 1e-5,
        carbon_dioxide = 1e-5,
        carbon_monoxide = 1e-5,
        hydrogen = 1e-5
    )=#
)

# Region: outlet (24 cells)
fixed_outlet = ComponentVector(
    permeability = 6.0e-12,
    mu = 1.0e-5,
    k = 237.0,
    cp = 3000.0
    #=species_diffusion_coeffs = ComponentVector(
        water = 1e-5,
        methanol = 1e-5,
        carbon_dioxide = 1e-5,
        carbon_monoxide = 1e-5,
        hydrogen = 1e-5
    )=#
)

# Region: wall (6382 cells)
fixed_wall = ComponentVector(
    k = 237.0,
    rho = 2700.0,
    cp = 921.0
)

# Region: heating_areas (178 cells)
fixed_heating_areas = ComponentVector(
    k = 237.0,
    rho = 2700.0,
    cp = 921.0
)

# Controller: temp_controller
fixed_temp_controller = ComponentVector(
        controllers = ComponentVector(
        integral_time = 1.0,
        initial_volumetric_input = 1e7,
        desired_value = ustrip(250.0u"°C" |> u"K"),
        max_volumetric_input = 1e8,
        min_volumetric_input = 0.0,
        proportional_gain = 1.0,
        derivative_time = 1.0
    )
)


# ── Do not modify below this line ────────────────────────────────────

initial_state = Dict{Symbol, Any}(
    :reforming_area => state_reforming_area,
    :inlet => state_inlet,
    :outlet => state_outlet,
    :wall => state_wall,
    :heating_areas => state_heating_areas,
    :temp_controller => state_temp_controller
)

fixed_properties = Dict{Symbol, Any}(
    :reforming_area => fixed_reforming_area,
    :inlet => fixed_inlet,
    :outlet => fixed_outlet,
    :wall => fixed_wall,
    :heating_areas => fixed_heating_areas,
    :temp_controller => fixed_temp_controller
)
