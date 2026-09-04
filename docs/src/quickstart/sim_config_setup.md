

# Getting Started with FVMFramework

This tutorial will teach you how to use the tools provided by FVMFramework to solve a 3D heat transfer problem

In this tutorial, we will be solving the following heat transfer equation:

```math
Q_{a, in} = -k A \frac{(T_{a}-T_{b})}{\Delta d}
```
where:
- $$Q_{a,in}$$ is the heat flux going into cell a in [W]
- $$k$$ is the thermal conductivity of the material in [W/(m*K)]
- $$A$$ is the area of the face shared by boths cells in [m^2]
- $$T_{a}$$ and $$T_{b}$$ are the temperatures of cell a and cell b respectively in [K]
- $$d$$ is the distance between the cell centers of each cell in [m]


```@docs
create_fvm_config
``` 

```julia
using Revise

using FVMFramework

using Ferrite
using FerriteGmsh
using OrdinaryDiffEq
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays

using Unitful
```

The general workflow of this package is to:
- Define a simple mesh using Ferrite.jl or import a mesh using FerriteGmsh.jl
- Create a config struct
- add cached fields, special caches, and optimized parameters
- add regions (groups of volumes with similar properties) and patches (groups of faces with similar properties)
- define the physics functions including physics functions for calculating fluxes between faces and fluxes between


```julia
grid_dimensions = (10, 10, 10)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Tetrahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))
```julia

```julia
u_proto = ComponentVector(
    temp = zeros(length(grid.cells))u"K",
)
```

config = create_fvm_config(grid, u_proto)

n_faces = length(config.geo.cell_neighbor_areas[1])

struct Solid <: AbstractPhysics end

#property updating/retrieval

#variable summation

#internal physics

#sources

#boundary conditions

#capacities

add_setup_syms!(
    config;
    cache_syms_and_units = (heat = u"J",),
    special_caches = ComponentVector(),
    second_order_syms = [],
    optimized_parameters = ComponentVector()
)

function cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    # J/s /= m^3 * kg*m^3 * J/(kg*K)
    # = K/s
    #@show du.heat[cell_id]
    #@show vol
    #@show u.rho[cell_id]
    #@show u.cp[cell_id]
    du.temp[cell_id] += du.heat[cell_id] / (vol * u.rho[cell_id] * u.cp[cell_id])
    #@show du.temp[cell_id]
end

function cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    # kg/s /= (m^3 / (J/(mol*K) * K))
    #remember: J = Pa*m^3
    # = Pa/s
    du_moles = du.mass[cell_id] / u.mw_avg[cell_id]
    du.pressure[cell_id] += (du_moles * u.R_gas[cell_id] * u.temp[cell_id]) / vol
end

function cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    total_mass = vol * u.rho[cell_id]

    for_fields!(du.mass_fractions, u.mass_fractions, du.species_mass_flows) do species, du_mass_fractions, u_mass_fractions, species_mass_flows
        du_mass_fractions[species[cell_id]] += (species_mass_flows[species[cell_id]] - u_mass_fractions[species[cell_id]] * du.mass[cell_id]) / total_mass
    end
end



add_region!(
    config, "copper";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 270.0u"°C",
    ),
    properties = ComponentVector(
        k = 237.0u"W/(m*K)", 
        rho = 2700.0u"kg/m^3",
        cp = 921.0u"J/(kg*K)",
    ),
    property_update_function = function copper_property_update!(properties, u)
        
    end,
    region_function =
    function heat_transfer!(du, u, cell_id, vol)
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "steel";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 350.0u"°C",
    ),
    properties = ComponentVector(
        k = 123.0u"W/(m*K)",
        rho = 7800.0u"kg/m^3",
        cp = 450.0u"J/(kg*K)",
    ),
    property_update_function = function steel_property_update!(properties, u)
        
    end,
    region_function =
    function heat_transfer!(du, u, p, t, cell_id, vol)
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

function heat_diffusion!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist

    du.heat[idx_a] -= -k_effective * grad_T * area
end

function solid_solid_flux!(
    du, u, p, t, 
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    heat_diffusion!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Solid && typeof(type_b) <: Solid && return solid_solid_flux!
end

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

function solve_system!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, p, t, geo, system)
    solve_controller_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = 0.0

test_prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1000.0), p_guess)
sol = solve(test_prob, Tsit5(), tspan = (0.0, 10.0))

t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#728.605 ms (341178 allocations: 1.08 GiB) (non-multithreaded)
#787.708 ms (202070 allocations: 1.07 GiB) (only connections multithreading)
#799.862 ms (219386 allocations: 1.07 GiB) (everything multithreaded)
#864.517 ms (211023 allocations: 1.07 GiB) (only regions multithreading)


VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess)

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end