using ComponentArrays
using PreallocationTools
using Polyester
using BenchmarkTools
using SparseConnectivityTracer
import ADTypes
using OrdinaryDiffEq
using FVMFramework

n_cells = 100000
n_faces = 6
reaction_names = (:WGS_rxn, :MD_rxn)
N = 2

du_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

du_axes = create_axes(du_proto, n_cells)
du_vec = Vector(ComponentVector(du_proto))

u_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_axes = create_axes(u_proto, n_cells)
u_vec = Vector(ComponentVector(u_proto))

du_cache_proto = ComponentVector(
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), 
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    rho = zeros(n_cells)
)

du_cache_axes = create_axes(du_cache_proto, n_cells)
du_cache_vec = Vector(ComponentVector(du_cache_proto))

u_cache_proto = ComponentVector(
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ),
    molar_concentrations = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    rho = zeros(n_cells)
)

u_cache_axes = create_axes(u_cache_proto, n_cells)
u_cache_vec = Vector(ComponentVector(u_cache_proto))

p_proto = ComponentVector(
    diffusion_pre_exponential_factor = [1e-10],
    diffusion_activation_energy = [10000.0]
)

p_axes = create_axes(p_proto, n_cells)
p_vec = Vector(ComponentVector(p_proto))

properties_proto = ComponentVector(
    viscosity = zeros(n_cells),
    rho = zeros(n_cells)
)

properties_axes = create_axes(properties_proto, n_cells)
properties_vec = Vector(ComponentVector(properties_proto))

# DiffCaches
du_diff_cache = DiffCache(du_cache_vec, N)
u_diff_cache = DiffCache(u_cache_vec, N)

cell_volumes = ones(n_cells)

function ode_for_testing_f!(
    du_vec, u_vec, p_vec, t,
    du_axes, u_axes,
    du_diff_cache, u_diff_cache,
    du_cache_axes, u_cache_axes,
    properties,
    cell_volumes
)
    du_cache_vec = get_tmp(du_diff_cache, first(u_vec) + first(p_vec))
    u_cache_vec = get_tmp(u_diff_cache, first(u_vec) + first(p_vec))

    #println(typeof(du_cache_vec))

    #do you know what is fucking crazy?
    #the only reason this worked in the past was because I zeroed out the u_cache_vec
    #otherwise, it just returns #undef for everything
    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = create_views_inline(du_cache_vec, du_cache_axes)
    u_cache_nt = create_views_inline(u_cache_vec, u_cache_axes)

    du_vec .= 0.0
    du_nt = create_views_inline(du_vec, du_axes)
    u_nt = create_views_inline(u_vec, u_axes)

    du = (; du_nt..., du_cache_nt...)
    u = (; properties..., u_nt..., u_cache_nt...) 

    u.rho .= (u_cache_nt.rho .+ properties.rho)

    for cell_id in eachindex(cell_volumes)
        # Direct dot syntax
        du.mass_fractions.methylene_blue[cell_id] += 1.0 
        
        #map(keys(du.mass_fractions)) do species_name
        #    du.mass_fractions[species_name][cell_id] += 1.0
        #    du.mass_fractions[species_name][cell_id] += 1.0
        #end

        du.mass_fractions.methylene_blue[cell_id] += 1.0
        du.mass_fractions.water[cell_id] += 1.0
        du.mass_fractions.methylene_blue[cell_id] += 1.0
        du.mass_fractions.water[cell_id] += 1.0
        
        #map(keys(du.reforming_reactions.net_rates.reforming_reactions)) do reaction_name
        #    du.reforming_reactions.net_rates.reforming_reactions[reaction_name][cell_id] += 1.0
        #end

        du.mass[cell_id] += du.rho[cell_id] * u.viscosity[cell_id]
    end
end

# Benchmark
println("Benchmarking ode_for_testing_f! with 'Best' Looping Style...")
@btime ode_for_testing_f!(
    $du_vec, $u_vec, $p_vec, 0.0,
    $du_axes, $u_axes,
    $du_diff_cache, $u_diff_cache,
    $du_cache_axes, $u_cache_axes,
    $properties_proto,
    $cell_volumes
)
#1.134 ms (6 allocations, 1.53 MiB)

f_closure = (du, u, p, t) -> ode_for_testing_f!(
    du, u, p, t,

    du_axes, u_axes,
    du_diff_cache, u_diff_cache,
    du_cache_axes, u_cache_axes,
    properties_proto,
    cell_volumes
)

length(du_vec)
length(u_vec)

jac_sparsity = SparseConnectivityTracer.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_vec, 0.0), du_vec, u_vec, SparseConnectivityTracer.TracerLocalSparsityDetector()
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 1.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u_vec, tspan, p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#@btime sol = solve(implicit_prob, FBDF())
#846.439 ms (505760 allocations, 903.6 MiB)

u_vec .= 0.0
explicit_prob = ODEProblem(f_closure, u_vec, tspan, p_vec)
#@btime sol = solve(explicit_prob, Tsit5())
#131.480 ms (41128 allocations: 185.25 MiB)