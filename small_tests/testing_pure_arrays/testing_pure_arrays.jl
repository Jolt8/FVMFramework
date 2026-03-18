using ComponentArrays
using Polyester
using BenchmarkTools
using StrideArraysCore

include("C://Users//wille//Desktop//FVMFramework//small_tests//in_house_solution_to_array_issues//FVMArray.jl")

StrideArray(undef, 2, 2)

n_cells = 100000

du_test = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

s, du_test_offset = create_axes(du_test)
du = zeros(du_test_offset)

Base.summarysize(du) #scales as you would expect
Base.summarysize(s) #very tiny, only scales with fields

function test_regular_array(du, s, n_cells)
    @batch for cell_id in 1:n_cells
        du[s.mass_fractions.methylene_blue][cell_id] += 1.0
    end
end

@btime test_regular_array($du, $s, $n_cells)

du_fvm = FVMArray(du, s)

function test_FVMArray(du_fvm, n_cells)
    @batch for cell_id in 1:n_cells
        du_fvm.mass_fractions.methylene_blue[cell_id] += 1.0
    end
end

@btime test_FVMArray($du_fvm, $n_cells)

function test_raw_array(du, n_cells)
    @batch for cell_id in 1:n_cells
        du[cell_id] += 1.0
    end
end

@btime test_raw_array($du, $n_cells)

#=
#things I have found out from this
    # - using the regular array pattern is much faster than using the FVMArray pattern and allows for allocation-free SIMD
    # - I don't know why the FVMArray pattern is slower, but it is
    # - I also don't know if there's a way to fix this
    # - based on the fact that FVM requires maximum performance, and the fact that the raw pattern of du[s.stuff...] isn't that ugly, 
    # - I think we're probably going to have to use the regular array pattern
    # - However, if we could overload something in StrideArraysCore to make it so that the FVMArray pattern works as intended, that would be awesome