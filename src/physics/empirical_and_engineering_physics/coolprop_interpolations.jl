temp_range = range(270.0, 570.0, length = 100)
pressure_range = range(100000.0, 1000000.0, length = 100)

temp_and_pressure_range = [(xi, yi) for xi in temp_range, yi in pressure_range]

function generate_coolprop_interpolation(
    interpolated_prop_name,
    prop_1_name::String, prop_1_min, prop_1_max, prop_1_n_samples,
    prop_2_name::String, prop_2_min, prop_2_max, prop_2_n_samples,
    mixture::String
)
    prop_1_range = range(prop_1_min, prop_1_max, length = prop_1_n_samples)
    prop_2_range = range(prop_2_min, prop_2_max, length = prop_2_n_samples)

    interpolated_prop = linear_interpolation(
        (prop_1_range, prop_2_range), 
        [PropsSI(interpolated_prop_name, prop_1_name, T, prop_2_name, P, mixture) for T in prop_1_range, P in prop_2_range]
    )

    return interpolated_prop
end
#=
prop_1_name = "T"
min_temp = 270.0
max_temp = 570.0
temp_n_samples = 100

prop_2_name = "P"
min_pressure = 100000.0
max_pressure = 1000000.0
pressure_n_samples = 100

mixture = "water"

dynamic_viscosity_interp_TP = generate_coolprop_interpolation(
    "V", 
    prop_1_name, min_temp, max_temp, temp_n_samples,
    prop_2_name, min_pressure, max_pressure, pressure_n_samples,
    mixture
)

k_interp_TP = generate_coolprop_interpolation(
    "conductivity", 
    prop_1_name, min_temp, max_temp, temp_n_samples,
    prop_2_name, min_pressure, max_pressure, pressure_n_samples,
    mixture
)

cp_interp_TP = generate_coolprop_interpolation(
    "C", 
    prop_1_name, min_temp, max_temp, temp_n_samples,
    prop_2_name, min_pressure, max_pressure, pressure_n_samples,
    mixture
)

cp_interp_TP(273, 100000)

function prandtl_number!(du, u, cell_id, vol)
    dynamic_viscosity = u.viscosity(u.temp, u.pressure)[1] #not sure why this returns a matrix 
    k = u.k(u.temp, u.pressure)[1]
    cp = u.cp(u.temp, u.pressure)[1]

    u.Pr[cell_id] = (dynamic_viscosity * cp) / (k)
end

u = (temp = [270.0], pressure = [100000.0], Pr = [0.0], viscosity = dynamic_viscosity_interp_TP, k = k_interp_TP, cp = cp_interp_TP)
du = (temp = [0.0], pressure = [0.0], Pr = [0.0])

prandtl_number!(du, u, 1, 1)
#@benchmark prandtl_number!($du, $u, 1, 1)
=#