using ComponentArrays
using PreallocationTools
using Polyester
using BenchmarkTools
using SparseConnectivityTracer
import ADTypes
using OrdinaryDiffEq
using NonlinearSolve
using FVMFramework

n_cells = 100000
n_faces = 6
reaction_names = (:WGS_rxn, :MD_rxn)
N = 2

du_proto = ComponentVector(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_proto = ComponentVector(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

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

p_proto = ComponentVector(
    diffusion_pre_exponential_factor = [1e-10],
    diffusion_activation_energy = [10000.0]
)

properties_proto = ComponentVector(
    viscosity = zeros(n_cells),
)

# Allocate backing storage
du_vec = Vector(du_proto)
u_vec = Vector(u_proto)
du_cache_vec = Vector(du_cache_proto)
u_cache_vec = Vector(u_cache_proto)
p_vec = Vector(p_proto)
properties_vec = Vector(properties_proto)

# DiffCaches
du_diff_cache_vec = DiffCache(du_cache_vec, N)
u_diff_cache_vec = DiffCache(u_cache_vec, N)

cell_volumes = ones(n_cells)

# Virtual Axes
virtual_du_axes = virtual_merge_axes((du_proto, du_cache_proto))
virtual_u_axes = virtual_merge_axes((u_proto, u_cache_proto, p_proto, properties_proto))

function ode_for_testing_f!(
    du_vec, u_vec, p_vec, t,
    v_du_axes, v_u_axes,
    du_diff_cache, u_diff_cache,
    props_vec,
    cell_vols,
)
    du_vec .= 0.0

    # Resolve merged views into VirtualFVMArrays
    if (first(u_vec) + first(p_vec)) isa SparseConnectivityTracer.Dual{Float64}
        u = VirtualFVMArray((u_vec, (get_tmp(u_diff_cache, first(u_vec) + first(p_vec)) .= 0.0), p_vec, props_vec), v_u_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache, first(u_vec) + first(p_vec)) .= 0.0)), v_du_axes)
    else
        u = VirtualFVMArray((u_vec, get_tmp(u_diff_cache, first(u_vec) + first(p_vec)), p_vec, props_vec), v_u_axes)
        du = VirtualFVMArray((du_vec, (get_tmp(du_diff_cache, first(u_vec) + first(p_vec)) .= 0.0)), v_du_axes)
    end

    @batch for cell_id in eachindex(cell_vols)
        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] += 1.0
        end
        du.pressure[cell_id] += u.pressure[cell_id]
    end
    return 
end

@btime ode_for_testing_f!(
    $du_vec, $u_vec, $p_vec, 0.0,
    $virtual_du_axes, $virtual_u_axes,
    $du_diff_cache_vec, $u_diff_cache_vec, 
    $properties_vec,
    $cell_volumes,
)

VSCodeServer.@profview ode_for_testing_f!(
    du_vec, u_vec, p_vec, 0.0,
    virtual_du_axes, virtual_u_axes,
    du_diff_cache_vec, u_diff_cache_vec, 
    properties_vec,
    cell_volumes,
)

f_closure = (du, u, p, t) -> ode_for_testing_f!(
    du, u, p, t,
    virtual_du_axes, virtual_u_axes,
    du_diff_cache_vec, u_diff_cache_vec, 
    properties_vec,
    cell_volumes
)

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

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
#@btime sol = solve(implicit_prob, FBDF())
#841.673 ms (260641 allocations: 863.35 MiB)
#1.127 s (394500 allocations: 870.32 MiB) (multithreaded)

u_vec .= 0.0
explicit_prob = ODEProblem(f_closure, u_vec, tspan, p_vec)
@btime sol = solve(explicit_prob, Tsit5())
#82.043 ms (32370 allocations: 134.37 MiB)
#183.951 ms (84321 allocations: 136.97 MiB) (multithreaded)

@btime sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
#719.297 ms (270879 allocations: 1.07 GiB)
#582.135 ms (458871 allocations: 1.09 GiB) (multithreaded)
