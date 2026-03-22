using BenchmarkTools
using Polyester
using ComponentArrays
using PreallocationTools
using SparseConnectivityTracer

include("virtual_ca_for_looping_base_methods.jl")
include("looping_methods.jl")

n_cells = 1000
n_faces = 6
reaction_names = (:WGS_rxn, :MD_rxn)
N = 2

# Prototypes for memory layout
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
        WGS_rxn = 0.0,
        MD_rxn = 0.0,
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
        WGS_rxn = 0.0,
        MD_rxn = 0.0,
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
    molecular_weights = ComponentVector(
        methylene_blue = 271.18,
        water = 18.015
    )
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
    u = VirtualFVMArray((u_vec, get_tmp(u_diff_cache, first(u_vec) + first(p_vec)), p_vec, props_vec), v_u_axes)
    du = VirtualFVMArray((du_vec, get_tmp(du_diff_cache, first(u_vec) + first(p_vec))), v_du_axes)

    for cell_id in eachindex(cell_vols)
        # Direct dot syntax
        du.mass_fractions.methylene_blue[cell_id] += 1.0 
        
        for_fields!(du.mass_fractions, du.molar_concentrations) do species, mass_fractions, molar_concentrations
            mass_fractions[species[cell_id]] += 1.0
            molar_concentrations[species[cell_id]] += 1.0
        end
        
        # Cross-buffer access
        du.mass[cell_id] += du.rho[cell_id] * u.viscosity[cell_id]

        du.mass_face[cell_id, 2] += 1.0
    end
    
    foreach_field_at!(du.net_rates) do rxn_idx, net_rates
        net_rates[rxn_idx] += 1.0
    end
end

# Benchmark
println("Benchmarking ode_for_testing_f! with 'Best' Looping Style...")
@btime ode_for_testing_f!(
    $du_vec, $u_vec, $p_vec, 0.0,
    $virtual_du_axes, $virtual_u_axes,
    $du_diff_cache_vec, $u_diff_cache_vec, 
    $properties_vec,
    $cell_volumes,
)