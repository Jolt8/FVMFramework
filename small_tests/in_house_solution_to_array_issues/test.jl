#Ok, so, I've been struggling for a while now to get a data structure that can support all the things shown off in this file
#thus, I think I'm going to have to create my own data structure 
#

using ComponentArrays 
using PreallocationTools
using Polyester
using LazyArrays
using BenchmarkTools
using SparseConnectivityTracer
import ADTypes
using OrdinaryDiffEq

function ode_for_testing_f!(
    du_vec, u_vec, p, t, 

    properties_vec,
    
    cell_volumes,

    p_axes,
    du_master_diff_cache_vec, u_master_diff_cache_vec,
    du_total_axes, u_total_axes
)
    #ABSOLUTELY REQUIRED


    #START OF SETUP
    #get_tmp cached variables (ABSOLUTELY REQUIRED)
    du_vec = get_tmp(du_master_diff_cache_vec, first(u_vec) + first(p))
    u_vec = get_tmp(u_master_diff_cache_vec, first(u_vec) + first(p))

    du_vec .= 0.0
    u_vec .= 0.0
    
    du = ComponentArray(du_vec, du_total_axes)
    u = ComponentArray(u_vec, u_total_axes)

    p_named = ComponentVector(p, p_axes)

    u.diffusion_pre_exponential_factor .= p_named.diffusion_pre_exponential_factor[1]
    u.diffusion_activation_energy .= p_named.diffusion_activation_energy[1]
    #END OF SETUP

    #CAPABILITIES THAT REQUIRE 0 ALLOCATIONS
    @batch for cell_id in 1:length(cell_volumes) #supports @batch
        #du.mass_fractions[:methylene_blue][cell_id] += 1.0 #symbolic indexing
        du.mass_fractions.methylene_blue[cell_id] += 1.0 #dot indexing

        du[1] #looks or is an array under the hood to @batch can use it

        u.net_rates.reforming_reactions.WGS_rxn += 1.0 #supports non cell_id indexed fields

        #=for species_name in keys(u.mass_fractions)
            du.mass_fractions[species_name][cell_id] += 1.0 #supports fast looping through fields without allocations
            for reaction in keys(u.net_rates.reforming_reactions) #supports nested field looping
                du.net_rates.reforming_reactions[reaction] += 1.0 #supports fast looping through single value fields without allocations
            end
        end=#

        for face_idx in 1:6
            du.mass_face[cell_id][face_idx] += 1.0 #supports fast looping through face indexed fields without allocations
        end

        du.mass[cell_id] += sum(du.mass_face[cell_id]) #supports basic array operations on fields
    end
    #END OF CAPABILITIES THAT REQUIRE 0 ALLOCATIONS
end

#=
other things to keep in mind
    - offset indexing like du.mass_face[cell_id + 1][5]
    - LazyMerging of du_in, cache, and properties 
    - ability to store interpolations like where u.viscosity would call u.viscosity(u.temp, u.pressure) under the hood
    - it's probably going to be called FVMArray or maybe something even more generic because I could see other people using this
=#

n_cells = 1000
n_faces = 6
reaction_names = (:WGS_rxn, :MD_rxn)

du_proto = ComponentArray(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u_proto = ComponentArray(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

du_caches = ComponentArray(
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

u_caches = ComponentArray(
    mass = zeros(n_cells),
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
    ), #it would also be nice if commans were not required for fields with 1 field within them
    rho = zeros(n_cells)
)

p = ComponentArray(
    diffusion_pre_exponential_factor = [1e-10],
    diffusion_activation_energy = [10000.0]
)

p_vec = Vector(p)
p_axes = getaxes(p)

properties = ComponentArray(
    rho = zeros(n_cells),
    viscosity = zeros(n_cells),
    diffusion_pre_exponential_factor = zeros(n_cells),
    diffusion_activation_energy = zeros(n_cells)
)

du_vec = Vector(du_proto)
u_vec = Vector(u_proto)
du_cache_vec = Vector(du_caches)
u_cache_vec = Vector(u_caches)
properties_vec = Vector(properties)

du_master_diff_cache_vec = DiffCache(vcat(du_vec, du_cache_vec))
u_master_diff_cache_vec = DiffCache(vcat(u_vec, u_cache_vec, properties_vec))

du_combined = ComponentArray(; NamedTuple(du_proto)..., NamedTuple(du_caches)...)
u_combined = ComponentArray(; NamedTuple(u_proto)..., NamedTuple(u_caches)..., properties...)

du_combined_axes = getaxes(du_combined)
u_combined_axes = getaxes(u_combined)

cell_volumes = ones(n_cells)

t = 0.0

du_view = view(get_tmp(du_master_diff_cache_vec, du_vec), 1:length(du_vec))
u_view = view(get_tmp(u_master_diff_cache_vec, u_vec), 1:length(u_vec))

@btime ode_for_testing_f!(
    du_view, u_view, p_vec, t,

    properties_vec, 
    
    cell_volumes,
    
    p_axes,
    du_master_diff_cache_vec, u_master_diff_cache_vec,
    du_combined_axes, u_combined_axes
)

f_closure = (du, u, p, t) -> ode_for_testing_f!(
    du, u, p, t,

    properties_vec, 
    
    cell_volumes,
    
    p_axes,
    du_master_diff_cache_vec, u_master_diff_cache_vec,
    du_combined_axes, u_combined_axes
)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_vec, 0.0), du_vec, u_vec, detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 100.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u_vec, tspan, p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF())