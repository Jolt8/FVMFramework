using ComponentArrays 
using PreallocationTools
using Polyester
using LazyArrays
using BenchmarkTools
using SparseConnectivityTracer
import ADTypes
using OrdinaryDiffEq
using NonlinearSolve
include("FVMArray.jl")

function ode_for_testing_f!(
    du, u, p_vec, t,

    du_merged_axes, u_merged_axes,
    
    cell_volumes,
)
    

    #=
    du_cache_vec = get_tmp(du_diff_cache_vec, u_vec)
    u_cache_vec = get_tmp(u_diff_cache_vec, u_vec) 

    du_cache_vec .= 0.0 
    u_cache_vec .= 0.0

    du_vec .= 0.0

    du_cache_fvm = FVMArray(du_cache_vec, du_caches_axes)
    u_cache_fvm = FVMArray(u_cache_vec, u_caches_axes)
    
    properties_fvm = FVMArray(properties_vec, properties_axes)
    p_fvm = FVMArray(p_vec, p_axes)

    #merge_into!(du_merged_buffer, (du_cache_fvm,))
    #merge_into!(u_merged_buffer, (u_cache_fvm, properties_fvm, p_fvm))

    du = FVMArray(du_merged_buffer, du_merged_axes)
    u = FVMArray(u_merged_buffer, u_merged_axes)

    #u.diffusion_pre_exponential_factor .= p_fvm.diffusion_pre_exponential_factor[1]
    #u.diffusion_activation_energy .= p_fvm.diffusion_activation_energy[1]
    =#

    du .= 0.0
    du = FVMArray(du, du_merged_axes)
    u = FVMArray(u, u_merged_axes)

    @batch for cell_id in 1:length(cell_volumes)
        du.mass_fractions[:methylene_blue][cell_id] += 1.0 
        du.mass_fractions.methylene_blue[cell_id] += 1.0 

        du.net_rates.reforming_reactions.WGS_rxn[1] += 1.0 

        foreach_field_at!(cell_id, du.mass_fractions) do species_name, mass_fraction
            #mass_fraction[species_name] += 1.0 
            #foreach_field_at!(1, du.net_rates.reforming_reactions) do reaction_name, net_rate #this allocates a ton when doing @batch
                #net_rate[reaction_name] += 1.0 
            #end
        end

        #foreach_field_at!(1, du.net_rates.reforming_reactions) do reaction_name, net_rate #oh, this is fine when used outside another foreach_field_at!
            #net_rate[reaction_name] += 1.0 
        #end

        for face_idx in 1:6
            du.mass_face[cell_id][face_idx] += 1.0 
        end

        du.mass[cell_id] += sum(du.mass_face[cell_id]) 
    end
end

n_cells = 10000
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

du_axes, du_offset = create_axes(du_proto)
du_vec = zeros(du_offset)

u_proto = (
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_axes, u_offset = create_axes(u_proto)
u_vec = zeros(u_offset)

du_caches = (
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

du_caches_axes, du_caches_offset = create_axes(du_caches)
du_caches_vec = zeros(du_caches_offset)
du_diff_cache_vec = DiffCache(du_caches_vec, N)

u_caches = (
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

u_caches_axes, u_caches_offset = create_axes(u_caches)
u_caches_vec = zeros(u_caches_offset)
u_diff_cache_vec = DiffCache(u_caches_vec, N)

p = (
    diffusion_pre_exponential_factor = [1e-10],
    diffusion_activation_energy = [10000.0]
)

p_axes, p_offset = create_axes(p)
p_vec = zeros(p_offset)

properties = (
    viscosity = zeros(n_cells),
    #diffusion_pre_exponential_factor = zeros(n_cells),
    #diffusion_activation_energy = zeros(n_cells)
)

properties_axes, properties_offset = create_axes(properties)
properties_vec = zeros(properties_offset)

du_merged_buffer = create_merged_buffer(length(du_vec) + length(du_caches_vec) + length(properties_vec) + length(p_vec))
du_merged_axes = merge_axes((du_axes, du_caches_axes, properties_axes, p_axes), (du_offset, du_caches_offset, properties_offset, p_offset))
u_merged_buffer = create_merged_buffer(length(u_vec) + length(u_caches_vec) + length(properties_vec) + length(p_vec))
u_merged_axes = merge_axes((u_axes, u_caches_axes, properties_axes, p_axes), (u_offset, u_caches_offset, properties_offset, p_offset))

cell_volumes = ones(n_cells)

t = 0.0

du_view = view(du_merged_buffer, 1:length(du_vec))
u_view = view(u_merged_buffer, 1:length(u_vec))

du_merged_fvm = FVMArray(du_merged_buffer, du_merged_axes)
u_merged_fvm = FVMArray(u_merged_buffer, u_merged_axes)

@btime ode_for_testing_f!(
    $du_merged_buffer, $u_merged_buffer, $p_vec, t,

    $du_merged_axes, $u_merged_axes,
    
    $cell_volumes,
) 
#1.185 ms (10 allocations: 736 bytes) (with @batch)

function test_raw_loop_performance(du_merged_buffer, u_merged_buffer)
    @batch for i in eachindex(du_merged_buffer)
        du_merged_buffer[i] += 1.0
        u_merged_buffer[i] += 1.0
    end
end

@btime test_raw_loop_performance(du_merged_buffer, u_merged_buffer)
#1.308 ms (0 allocations: 0 bytes) for 100000 cells (with @batch)
#1.086 ms (0 allocations: 0 bytes) for 100000 cells (without @batch)
#1.102 ms (58 allocations: 3.50 KiB) for 100000 cells (with #threads)


f_closure = (du, u, p, t) -> ode_for_testing_f!(
    du, u, p, t,

    du_merged_axes, u_merged_axes,
    
    cell_volumes,
)

detector = SparseConnectivityTracer.TracerSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_vec, 0.0), du_merged_buffer, u_merged_buffer, detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 1.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u_merged_buffer, tspan, p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

using ILUZero

function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end

VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
@btime sol = solve(implicit_prob, FBDF())
#=
explicit_prob = ODEProblem(f_closure, u_merged_buffer, tspan, p_vec)
@btime sol = solve(explicit_prob, Tsit5())

findall(x -> x != 0.0, sol.u[1])

u_named = FVMArray(sol.u[1], u_axes)

u_named.mass_fractions.methylene_blue
