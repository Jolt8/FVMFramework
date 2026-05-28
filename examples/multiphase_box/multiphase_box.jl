using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
# Overload power_by_squaring to support tracer types inside Clapeyron's complex root solver
#Base.power_by_squaring(x::Any, p::SparseConnectivityTracer.Dual) = p
Base.power_by_squaring(x::SparseConnectivityTracer.Dual, p::SparseConnectivityTracer.Dual) = x+p
#Base.power_by_squaring(x::Any, p::SparseConnectivityTracer.Dual) = x^p #this is the only thing left that fails, the above works though

# Overload ldexp for SparseConnectivityTracer.Dual
Base.ldexp(x::SparseConnectivityTracer.Dual, i::Int) = x

# Define explicit promotion rules to resolve ambiguity
Base.promote_rule(::Type{SparseConnectivityTracer.Dual{P, T}}, ::Type{ForwardDiff.Dual{Tx, Vx, Nx}}) where {P, T, Tx, Vx, Nx} = 
    ForwardDiff.Dual{Tx, promote_type(SparseConnectivityTracer.Dual{P, T}, Vx), Nx}

Base.promote_rule(::Type{ForwardDiff.Dual{Tx, Vx, Nx}}, ::Type{SparseConnectivityTracer.Dual{P, T}}) where {P, T, Tx, Vx, Nx} = 
    ForwardDiff.Dual{Tx, promote_type(SparseConnectivityTracer.Dual{P, T}, Vx), Nx}

# Define convert overload to resolve conversion ambiguity
Base.convert(::Type{ForwardDiff.Dual{Tx, Vx, Nx}}, y::D) where {Tx, Vx, Nx, P, T, D <: SparseConnectivityTracer.Dual{P, T}} = 
    ForwardDiff.Dual{Tx, Vx, Nx}(convert(Vx, y), ForwardDiff.Partials{Nx, Vx}(ntuple(i -> zero(Vx), Nx)))

# Resolve ambiguities for common binary operators and comparisons (+, -, *, /, ^, ==, <, <=)
for TracerType in (SparseConnectivityTracer.GradientTracer, SparseConnectivityTracer.HessianTracer)
    # Define ^(::D, ::D) using the exp(y * log(x)) identity
    @eval begin
        function Base.:^(x::D, y::D) where {P, T <: $TracerType, D <: SparseConnectivityTracer.Dual{P, T}}
            return exp(y * log(x))
        end
    end
    for op in (:+, :-, :*, :/, :^, Symbol("=="), Symbol("<"), Symbol("<="))
        @eval begin
            function Base.$op(x::ForwardDiff.Dual{Tx, Vx, Nx}, y::D) where {Tx, Vx, Nx, P, T <: $TracerType, D <: SparseConnectivityTracer.Dual{P, T}}
                x_prom, y_prom = promote(x, y)
                return Base.$op(x_prom, y_prom)
            end
            function Base.$op(x::D, y::ForwardDiff.Dual{Tx, Vx, Nx}) where {Tx, Vx, Nx, P, T <: $TracerType, D <: SparseConnectivityTracer.Dual{P, T}}
                x_prom, y_prom = promote(x, y)
                return Base.$op(x_prom, y_prom)
            end
        end
    end
end

using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak

#=
using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using ForwardDiff
using DataFrames
using CSV
using DataInterpolations
=#
using Dates

using FVMFramework

grid_dimensions = (1, 1, 100)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 20.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)
#=
addcellset!(grid, "wall", xyz -> 
    xyz[1] <= 1.0 &&
    xyz[1] >= 99.0 &&
    xyz[2] <= 1.0 &&
    xyz[2] >= 99.0 &&
    xyz[3] <= 1.0 &&
    xyz[3] >= 99.0
)
getcellset(grid, "wall")

addcellset!(grid, "fluid", xyz -> 
    xyz[1] < 99.0 &&
    xyz[1] > 1.0 &&
    xyz[2] < 99.0 &&
    xyz[2] > 1.0 &&
    xyz[3] < 99.0 &&
    xyz[3] > 1.0
)
getcellset(grid, "fluid")
=#

addcellset!(grid, "fluid", xyz -> true)

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

#properties
Revise.includet(joinpath(@__DIR__, "properties", "water_methanol_properties.jl"))
water_methanol_properties = get_water_methanol_properties()

#physics
Revise.includet(joinpath(@__DIR__, "physics", "drift_flux_vaporization.jl"))

function get_mole_fractions_vec(du, u, cell_id, vol)
    gas_mole_fractions_vec = [
        max(u.gas_mole_fractions.methanol[cell_id], 1e-15),
        max(u.gas_mole_fractions.water[cell_id], 1e-15)
    ]
    liquid_mole_fractions_vec = [
        max(u.liquid_mole_fractions.methanol[cell_id], 1e-15),
        max(u.liquid_mole_fractions.water[cell_id], 1e-15)
    ]
    # Normalize so they sum to 1
    gas_mole_fractions_vec ./= sum(gas_mole_fractions_vec)
    liquid_mole_fractions_vec ./= sum(liquid_mole_fractions_vec)
    return gas_mole_fractions_vec, liquid_mole_fractions_vec
end

function update_solid_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id] * u.steel_thermal_mass_multiplier[cell_id]
    u.rho[cell_id] = properties.solid_rho[cell_id]
end

#clapeyron_model = CPA(["methanol", "water"])
clapeyron_model = PR(["methanol", "water"])

Clapeyron.volume(clapeyron_model, 100000, 273.13, [0.5, 0.5], phase = :liquid)
Clapeyron.volume(clapeyron_model, 100000, 273.13, [1.0, 1.0], phase = :liquid)

function update_fluid_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model)
    
    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = u.bed_porosity[cell_id] * u.fluid_rho[cell_id] + (1.0 - u.bed_porosity[cell_id]) * properties.solid_rho[cell_id]
end

function fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    #it's definitely one of these
    sum_mass_flux_face_to_cell!(du, u, cell_id) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
end

function solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

u_proto = ComponentVector(
    gas_densities = (
        methanol = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    liquid_densities = (
        methanol = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    #I think overall mass fractions would be a cached variable then
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
#reaction_names = keys(water_methanol_properties.reactions.reforming_reactions)
species_names = keys(u_proto.gas_densities)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)", #k and cp cause a dimsnion error for some reason
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        mass = u"kg",
        gas_holdup = u"1",
        liquid_holdup = u"1",
        gas_enthalpy = u"J/kg",
        liquid_enthalpy = u"J/kg",
        overall_gas_density = u"kg/m^3",
        overall_liquid_density = u"kg/m^3",
        fluid_rho = u"kg/m^3",

        gas_mw_avg = u"kg/mol",
        liquid_mw_avg = u"kg/mol",

        eos_gas_density = u"kg/m^3", #these are the intrinsic density of the current phase found with the eos model
        eos_liquid_density = u"kg/m^3",
        
        mass_face = u"kg",
        
        viscous_drag = u"kg/(m*s)",
        inertial_drag = u"kg/(m^3)",
        driving_force = u"N/m^3",
        mixture_superficial_velocity = u"m/s",
        mixture_pore_velocity = u"m/s",

        gas_velocity_face = u"m/s",
        liquid_velocity_face = u"m/s",

        gas_mass_fractions = u"1",
        liquid_mass_fractions = u"1",

        gas_mole_fractions = u"1",
        liquid_mole_fractions = u"1",

        gas_generation = u"kg/(m^3*s)",

        K_vle = u"1",
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",

        viscous_drag = zeros(n_cells, n_faces)u"kg/(m*s)",
        inertial_drag = zeros(n_cells, n_faces)u"kg/(m^3)",
        driving_force = zeros(n_cells, n_faces)u"N/m^3",
        mixture_superficial_velocity = zeros(n_cells, n_faces)u"m/s",
        mixture_pore_velocity = zeros(n_cells, n_faces)u"m/s",

        gas_velocity_face = zeros(n_cells, n_faces)u"m/s",
        liquid_velocity_face = zeros(n_cells, n_faces)u"m/s",

        gas_mass_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),
        liquid_mass_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),

        gas_mole_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),
        liquid_mole_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),

        gas_generation = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"kg/(m^3*s)" for _ in 1:length(species_names))
        ),

        K_vle = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector()
)

add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = ComponentVector(
        gas_densities = ComponentVector(
            methanol = 0.1u"kg/m^3",
            water = 0.2u"kg/m^3",
        ),
        liquid_densities = ComponentVector(
            methanol = 80.0u"kg/m^3",
            water = 100.0u"kg/m^3",
        ),
        pressure = 1.0u"atm",
        temp = 60.0u"°C",
    ),
    properties = water_methanol_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        #update_fluid_properties!(du, u, cell_id, vol, system)

        #du.heat[cell_id] += 10.0

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)
#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    overall_drift_flux_model!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )

    #this is the only equation that's still needed
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function fluid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function solid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_solid_flux!

    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return fluid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!
end

function solve_system!(du, u, p_vec, t, geo, system)
    for cell_id in grid.cellsets["fluid"]
        update_fluid_properties!(du, u, cell_id, geo.cell_volumes[cell_id], system)
    end

    solve_connection_groups!(du, u, geo, system)
    solve_controller_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system) #this seems to be the culprit
end

du0_vec, u0_vec, state_axes, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = 0.0

#test_prob = ODEProblem(f_closure, u0_vec, (0.0, 0.01), p_guess)
#@time sol = solve(test_prob, Tsit5(), callback = approximate_time_to_finish_cb)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 30000.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

@time sol = solve(implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)];
u_named[33].pressure
u_named[33].gas_densities.methanol
u_named[33].liquid_densities.water

root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

sol_to_vtk(sol, u_named, grid, @__FILE__, root_dir)

hi = 1

#=

using Plots

gas_inlet_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.methanol[1] for i in 1:length(sol.u)], label = "methanol", title = "gas_densities inlet")
plot!(gas_inlet_plot, [sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.water[1] for i in 1:length(sol.u)], label = "water")

gas_middle_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.methanol[50] for i in 1:length(sol.u)], label = "methanol", title = "gas_densities middle")
plot!(gas_middle_plot, [sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.water[50] for i in 1:length(sol.u)], label = "water")

gas_outlet_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.methanol[100] for i in 1:length(sol.u)], label = "methanol", title = "gas_densities outlet")
plot!(gas_outlet_plot, [sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.water[100] for i in 1:length(sol.u)], label = "water")

liquid_inlet_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.methanol[1] for i in 1:length(sol.u)], label = "methanol", title = "liquid_densities inlet")
plot!(liquid_inlet_plot, [sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.water[1] for i in 1:length(sol.u)], label = "water")

liquid_middle_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.methanol[50] for i in 1:length(sol.u)], label = "methanol", title = "liquid_densities middle")
plot!(liquid_middle_plot, [sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.water[50] for i in 1:length(sol.u)], label = "water")

liquid_outlet_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.methanol[100] for i in 1:length(sol.u)], label = "methanol", title = "liquid_densities outlet")
plot!(liquid_outlet_plot, [sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.water[100] for i in 1:length(sol.u)], label = "water")

methanol_liquid_cell_before_outlet_plot = plot([sol.t[i] for i in 270:length(sol.u)], [u_named[i].liquid_densities.methanol[99] for i in 270:length(sol.u)], label = "methanol liquid density", title = "liquid and gas densities before outlet") #this is the strange one
plot!([sol.t[i] for i in 270:length(sol.u)], [u_named[i].gas_densities.methanol[99] for i in 270:length(sol.u)], label = "methanol gas density") #this is also a strange one
water_liquid_cell_before_outlet_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].liquid_densities.water[99] for i in 1:length(sol.u)], label = "water", title = "liquid_densities before outlet")
water_gas_cell_before_outlet_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].gas_densities.water[99] for i in 1:length(sol.u)], label = "water", title = "liquid_densities before outlet")

inlet_pressure_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].pressure[1] for i in 1:length(sol.u)], label = "inlet", title = "pressure inlet")
middle_pressure_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].pressure[50] for i in 1:length(sol.u)], label = "pressure", title = "pressure middle")
outlet_pressure_plot = plot([sol.t[i] for i in 1:length(sol.u)], [u_named[i].pressure[100] for i in 1:length(sol.u)], label = "pressure", title = "pressure outlet")

plot(gas_inlet_plot, gas_middle_plot, gas_outlet_plot, liquid_inlet_plot, liquid_middle_plot, liquid_outlet_plot, methanol_liquid_cell_before_outlet_plot, water_liquid_cell_before_outlet_plot, inlet_pressure_plot, middle_pressure_plot, outlet_pressure_plot, layout = (11, 1), legend = false, size = (3000, 3000))
=#