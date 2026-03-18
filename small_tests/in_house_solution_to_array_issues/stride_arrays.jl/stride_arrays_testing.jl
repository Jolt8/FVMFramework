using StrideArraysCore
using BenchmarkTools
using Polyester
using ComponentArrays
include("C://Users//wille//Desktop//FVMFramework//small_tests//in_house_solution_to_array_issues//FVMArray.jl")

n_cells = 1000

u = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_axes, u_offset = create_axes(u)
u_vec = StrideArray{Float64}(zero, u_offset)

u = FVMArray(u_vec, u_axes)

u.ptr

function test_strided_array(u, n_cells)
    @batch for i in 1:n_cells
        u.mass_fractions.methylene_blue[i] += 1.0
    end
end

u = ComponentArray(mass_fractions = (methylene_blue = zeros(n_cells),),)

@btime test_strided_array($u, $n_cells)
