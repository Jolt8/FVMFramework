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

function ode_for_testing_f!(du_vec, u_vec, p, t)
    du_vec .= 0.0

    @batch for cell_id in 1:100000
        du_vec[cell_id + 100] += 1.0
        du_vec[cell_id + 200] += 1.0
        du_vec[cell_id] += u_vec[cell_id]
    end
    return 
end

jac_sparsity = SparseConnectivityTracer.jacobian_sparsity(
    (du, u) -> ode_for_testing_f!(du, u, p_vec, 0.0), du_vec, u_vec, SparseConnectivityTracer.TracerLocalSparsityDetector()
)

ode_func = ODEFunction(ode_for_testing_f!, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 1.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u_vec, tspan, p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
#@time sol = solve(implicit_prob, FBDF())
#819.045 ms (255688 allocations: 863.02 MiB)
#1.104 s (400281 allocations: 870.58 MiB) (multithreaded)

u_vec .= 0.0
explicit_prob = ODEProblem(ode_for_testing_f!, u_vec, tspan, p_vec)
@btime sol = solve(explicit_prob, Tsit5())
#81.686 ms (33877 allocations: 171.06 MiB) (non-multithreaded)
#70.306 ms (265 allocations: 169.38 MiB) (multithreaded)

@btime sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
#1.642 s (4181 allocations: 2.40 GiB) (non-multithreaded)
#1.648 s (4181 allocations: 2.40 GiB) (multithreaded)