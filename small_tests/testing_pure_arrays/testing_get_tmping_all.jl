using ComponentArrays
using Polyester
using BenchmarkTools
using PreallocationTools
include("C://Users//wille//Desktop//FVMFramework//small_tests//in_house_solution_to_array_issues//FVMArray.jl")

n_cells = 10000000

du_test = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

s, du_test_offset = create_axes(du_test)
du = zeros(du_test_offset)

du_diff_cache = DiffCache(du, 2)

Base.summarysize(du) #scales as you would expect
Base.summarysize(s) #very tiny, only scales with fields

function test_diff_cached_array(du_diff_cache, s, n_cells)
    du = get_tmp(du_diff_cache, 1.0)
    #du = view(get_tmp(du_diff_cache, 1.0), 1:n_cells) #somehow this is faster
    @batch for cell_id in 1:n_cells
        du[s.mass_fractions.methylene_blue][cell_id] += 1.0
    end
end

@time test_diff_cached_array(du_diff_cache, s, n_cells)
@btime test_diff_cached_array($du_diff_cache, $s, $n_cells) 
##somehow slightly faster than raw

function test_raw_array(du, n_cells)
    @batch for cell_id in 1:n_cells
        du[cell_id] += 1.0
    end
end

@time test_raw_array(du, n_cells)
@btime test_raw_array($du, $n_cells)