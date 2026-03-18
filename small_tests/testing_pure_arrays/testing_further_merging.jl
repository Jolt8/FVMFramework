using ComponentArrays
using Polyester
using BenchmarkTools
using PreallocationTools
include("C://Users//wille//Desktop//FVMFramework//small_tests//in_house_solution_to_array_issues//FVMArray.jl")

n_cells = 100000

u_state = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

reaction_names = (:WGS_rxn, :MD_rxn)

u_cached = (
    rho = zeros(n_cells),
    mass_face = fill(zeros(6), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ),
)

u_properties = (
    viscosity = zeros(n_cells),
)

s, u_offset = create_axes((; u_state..., u_cached..., u_properties...))
u = zeros(u_offset)

u_diff_cache = DiffCache(u, 2)

_, state_offset = create_axes(u_state)
_, cached_offset = create_axes(u_cached)
_, properties_offset = create_axes(u_properties)

properties_unit_range = (state_offset + cached_offset + 1) : u_offset
properties_vec = zeros(length(properties_unit_range))

Base.summarysize(u) #scales as you would expect
Base.summarysize(s) #very tiny, only scales with fields

function test_diff_cached_array(u_diff_cache, s, properties_vec, properties_unit_range, n_cells)
    u = get_tmp(u_diff_cache, 1.0)
    #u[1000003:1100002] = properties_vec #this is what's causing the slowdown
    @batch for cell_id in 1:n_cells
    end
    return 
end

@time test_diff_cached_array(u_diff_cache, s, properties_vec, properties_unit_range, n_cells)
@btime test_diff_cached_array($u_diff_cache, $s, $properties_vec, $properties_unit_range, $n_cells) #45.7 μs

function test_raw_array(u, n_cells)
    @batch for cell_id in 1:n_cells
    end
end

@time test_raw_array(u, n_cells)
@btime test_raw_array($u, $n_cells) #12.7 μs