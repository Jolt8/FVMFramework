using Interpolations
using CoolProp
using BenchmarkTools

temp_range = range(270.0, 570.0, length = 100)
pressure_range = range(100000.0, 1000000.0, length = 100)

temp_and_pressure_range = [(xi, yi) for xi in temp_range, yi in pressure_range]
corresponding_viscosity = [PropsSI("V", "T", T, "P", P, "water") for T in temp_range, P in pressure_range]

viscosity = linear_interpolation((temp_range, pressure_range), corresponding_viscosity)

function coolprop_test!(temp, pressure)
    PropsSI("V", "T", temp, "P", pressure, "water")
end

function interpolated_test!(temp::Float64, pressure::Float64)
    viscosity(temp, pressure)
end

temp = 332.0
pressure = 121938.0

@btime for i in 1:10 coolprop_test!($temp, $pressure) end #4.864 ms
@btime for i in 1:10 interpolated_test!($temp, $pressure) end #1.030 μs