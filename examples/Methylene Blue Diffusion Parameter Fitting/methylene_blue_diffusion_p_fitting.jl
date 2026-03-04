using Revise

using FVMFramework

#I think this take a long time because Enzyme (despite not being installed) throws a warning from DiffEqBaseEnzyme.jl

using Ferrite
using FerriteGmsh
using OrdinaryDiffEq
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays
using PreallocationTools
using ForwardDiff
using Polyester
using DataInterpolations

using Unitful

mesh_path = joinpath(@__DIR__, "cone_end/dialysis_tubing_cone_output.msh")

grid = togrid(mesh_path)

grid.cellsets
grid.facetsets

n_cells = length(grid.cells)
u_proto = (
    mass_fractions = (methylene_blue = zeros(n_cells), water = zeros(n_cells)),
)

config = create_fvm_config(grid, u_proto)

#=
for (cell_id, facet_idx) in grid.facetsets["dialysis_tubing_surface"]
    println(cell_id, " ", facet_idx)
end
=#

check_cellset_connectivity(config.grid, "dialysis_tubing_interior")



#since each mass fraction is modified, they have to be vectors


function normalize_mass_fractions(mass_fractions)
    total_mass_fractions = 0.0

    for (species_name, mass_fraction) in pairs(mass_fractions)
        total_mass_fractions += mass_fractions[species_name][1]
    end

    for species_name in keys(mass_fractions)
        mass_fractions[species_name][1] /= total_mass_fractions
    end

    return NamedTuple{keys(mass_fractions)}(first.(values(mass_fractions)))
end

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

dialysis_tubing_initial_mass_fractions = (
    methylene_blue = [0.0003846],
    water = [0.9996154]
)

dialysis_tubing_initial_mass_fractions = normalize_mass_fractions(dialysis_tubing_initial_mass_fractions)

add_region!(
    config, "dialysis_tubing_interior";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = dialysis_tubing_initial_mass_fractions,
    ),
    properties = (
        temp = ustrip(21.13u"°C" |> u"K"),
        k = 0.6, # k (W/(m*K))
        cp = 4184, # cp (J/(kg*K))
        rho = 1000, # rho (kg/m^3)
        mu = 1e-3, # mu (Pa*s)
        species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
        diffusion_coefficients = (
            methylene_blue = 1e-5,
            water = 1e-5
        ), #diffusion coefficients (m^2/s)
        molecular_weights = (
            methylene_blue = 0.31985, 
            water = 0.01802
        ), #species_molecular_weights [kg/mol]
    ), 
    optimized_syms = [],
    state_syms = [:mass_fractions],
    cache_syms = [:heat, :molar_concentrations, :mw_avg, :rho], 
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

surrounding_fluid_initial_mass_fractions = (
    methylene_blue = [0.0],
    water = [1.0]
)

surrounding_fluid_initial_mass_fractions = normalize_mass_fractions(surrounding_fluid_initial_mass_fractions)

add_region!(
    config, "surrounding_fluid";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = surrounding_fluid_initial_mass_fractions,
    ),
    properties = (
        temp = ustrip(21.13u"°C" |> u"K"),
        k = 0.6, # k (W/(m*K))
        cp = 4184, # cp (J/(kg*K))
        rho = 1000, # rho (kg/m^3)
        mu = 1e-3, # mu (Pa*s)
        species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
        diffusion_coefficients = (
            methylene_blue = 1e-5,
            water = 1e-5
        ), #diffusion coefficients (m^2/s)
        molecular_weights = (
            methylene_blue = 0.31985, 
            water = 0.01802
        ), #species_molecular_weights [kg/mol]
    ), 
    optimized_syms = [],
    state_syms = [:mass_fractions],
    cache_syms = [:heat, :molar_concentrations, :mw_avg, :rho], 
    region_function =
    function surrounding_fluid!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

include("specific_physics/arrenhius_mass_fraction_diffusion_meth_blue_and_water.jl") 
#for arrenhius_mass_fraction_diffusion_meth_blue_and_water!

#this does nothing for now, but it will probably be used to determine the length of the p_vector in the future

#we'll come back to this later, but I'm thinking that instead I'll just add optimized_ in front of 
#the property name to denote that it's an optimized parameter
#=
macro optimize(key_value_pair) 
    param_name = key_value_pair.args[1] #this seem to always return the correct name
    param_value = key_value_pair.args[2]

    #=
    if typeof(key_value_pair.args[2]) <: Expr
        param_value = key_value_pair.args[2].args
    else
        param_value = key_value_pair.args[2]
    end
    =#

    #println(param_name)
    println(param_value)
    
    push!(p_tracker, OptimizedParameterTracker(param_name, param_value))
    
    return :($param_name => $param_value[1])
end

properties = (
    @optimize :test = 1,
    @optimize :diffusion_pre_exponential_factor = 1e-5,
    @optimize :diffusion_activation_energy = 1000.0
)
=#

add_patch!(
    config, "dialysis_tubing_surface";
    properties = (
        diffusion_pre_exponential_factor = 1e-5,
        diffusion_activation_energy = 1000.0
    ), 
    optimized_syms = [
        :diffusion_pre_exponential_factor,
        :diffusion_activation_energy
    ],
    patch_function =
    function membrane_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        
        #zeroing out regular diffusion that would happen across this face
        #du.mass_fractions[:methylene_blue][idx_a] *= 0.0
        #du.mass_fractions[:water][idx_a] *= 0.0
        #du.mass_fractions[:methylene_blue][idx_b] *= 0.0
        #du.mass_fractions[:water][idx_b] *= 0.0

        arrenhius_mass_fraction_diffusion_meth_blue_and_water!(
            du, u,
            idx_a, idx_b, face_idx,
            cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
            cell_volumes[idx_a], cell_volumes[idx_b]
        )
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    #=
    temp_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

#this is the smallest I could make this function
#I like using <: here because it makes it look nice with syntax highlighting
function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
)

#START OF OPTIMIZATION SETUP (just to help differentiate between what is and isn't optimization related)
#after doing some work on this, it might be best to make it so that custom methods are used every single time
#I think making something that is generalizable to any problem would be almost impossible to code
#So, we'll just have some basic methods that help with setup but then everything else requires editing the operator itself 

#=stuff we still have to add
    - methods for assigning experimental data only to certain cellsets 
        - I wonder if the memory cost from passing in experiemental data of n_fields * n_timestamps * n_cells is too much
        - For example, if we had 1 field, 4000 timesteps, and 10000 cells (pretty typical), that's 40,000,000 data points
        - However if we assume the field is constant for each field, then we only need n_fields * n_timestamps data points
        - We could also reduce the resolution of the timesteps, but that's not ideal
    - 
=#
using XLSX

experimental_data_path = joinpath(@__DIR__, "dialysis_tubing_diffusion_data.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

#start of logic that doesn't need to be seen by the user

struct TrialData
    state_data::NamedTuple
    state_time::NamedTuple
    compared_data::NamedTuple
    compared_time::NamedTuple
end

function initialize_trials()
    return Dict{String, TrialData}()
end

function input_experimental_data(trials, xf, name::String; 
    data::NamedTuple, 
    trial_processor_function::Function
)  
    trial_processor_function(trials, data, name, xf)
end

function dv(data) #column_to_vector/row_to_vector, called it dv because it's data_to_vector
    return float.(vec(data))
end
#end of logic that doesn't need to be seen by the user

function moles_to_mass_fraction(species_moles, species_molecular_weight, rho, mixture_volume)
    return ((species_moles / mixture_volume) * species_molecular_weight) / rho
end

function process_experimental_data!(trials, data::NamedTuple, name, xf)
    #state data (affects the conditions of the simulation over time)
    temperatures = data.temp_data.temp .* u"°C"
    state_data = (
        temp = CubicSpline(ustrip.((temperatures .|> u"K")), data.temp_data.timestamps),
    )

    state_time = (
        temp = data.temp_data.timestamps, #this is only needed for defining tMax
    )



    #compared data (what the simulated data is compared against to get loss)
    mixture_volume = 500u"ml"
    mixture_rho = 1000.0u"kg/m^3"

    calibration_concentration = 1.6u"g/L"
    calibration_volume_fractions = dv(xf["calibration"]["A2:A6"]) .* u"ml/ml"
    stock_solution_density = 1015.0u"kg/m^3"

    calibration_mass_fractions = (calibration_volume_fractions .* calibration_concentration) / stock_solution_density
    
    corresonding_absorbances = dv(xf["calibration"]["B2:B6"])

    calibration_absorbances_to_mass_fractions_interp = CubicSpline(ustrip.(calibration_mass_fractions .|> u"kg/kg"), corresonding_absorbances)

    methylene_blue_mass_fractions = calibration_absorbances_to_mass_fractions_interp.(data.absorbance_data.methylene_blue_absorbance)
    water_mass_fractions = 1.0 .- methylene_blue_mass_fractions

    compared_data = (
        mass_fractions = (
            methylene_blue = ustrip.((methylene_blue_mass_fractions .|> u"kg/kg")),
            water = ustrip.((water_mass_fractions .|> u"kg/kg")),
        ),
    )

    compared_time = (
        mass_fractions = data.absorbance_data.timestamps,
    )

    trial = TrialData(
        state_data,
        state_time,
        compared_data,
        compared_time
    )

    trials[name] = trial
end

trials = initialize_trials()

input_experimental_data(trials, xf, "T1"; 
    data = (
        temp_data = (
            timestamps = dv(xf["T1 Temps"]["A2:A805"]),
            temp = dv(xf["T1 Temps"]["B2:B805"]),
        ),
        absorbance_data = ( #I think it would be safe to assume that values that are measured in the same field have the same timestamps
            timestamps = dv(xf["T1"]["A2:A6"]),
            methylene_blue_absorbance = dv(xf["T1"]["B2:B6"]),
        )
    ),
    trial_processor_function = process_experimental_data!
)

#END OF OPTIMIZATION SETUP

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, special_caches)

u_test = (; create_views_inline(u0_vec, system.u_proto_axes)..., create_views_inline(get_tmp(system.u_diff_cache_vec, 0.0), system.u_cache_axes)...
)

f_closure_implicit_pre = (du, u, p, t, state_data, state_time) -> 
methylene_blue_diffuion_parameter_fitting_f!(
    du, u, p, t, 

    state_data, state_time,
    system.p_axes,
    
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals, 

    system.connection_groups, system.controller_groups, system.patch_groups, system.region_groups,

    system.merged_properties, system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system
f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trials["T1"].state_data, trials["T1"].state_time)

using BenchmarkTools

tMax = trials["T1"].state_time.temp[end]

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, tMax), system.p_vec)
#@time sol = solve(prob, Tsit5(), saveat = tMax / 100)#, callback = approximate_time_to_finish_cb)

#VSCodeServer.@profview sol = solve(prob, Tsit5())

#sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)

#sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

#sol_u_named_0.mass_fractions.methylene_blue == sol_u_named_end.mass_fractions.methylene_blue

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, system.p_vec, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
#tMax is already defined above based on experimental data
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, system.p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization()))

sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)

sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

sol_u_named_0.mass_fractions.methylene_blue == sol_u_named_end.mass_fractions.methylene_blue

VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization()))
#VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

#observed_cell_id = grid.cellsets["molar_concentrations_observation_point"][1] #this is what we should do
observed_cell_id = grid.cellsets["surrounding_fluid"][1] #this is just temporary

using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ForwardDiff

#START OF OPTIMIZATION SOLVING
function loss(θ)
    total_loss = 0.0
    for trial_name in keys(trials)
        trial = trials[trial_name]

        tspan = (trial.state_time.temp[1], trial.state_time.temp[end])
        p_adjusted = [exp(θ[1]), θ[2]]

        f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.state_data, trial.state_time)

        ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
        #since this is just a diffusion problem, I doubt implicit solving is necessary

        prob_trial = ODEProblem(ode_func, u0_vec, tspan, system)

        sol = solve(
            prob_trial, FBDF(linsolve = KLUFactorization()), p = p_adjusted,
            sensealg = InterpolatingAdjoint(autodiff = AutoForwardDiff()),
            saveat = trial.compared_time.mass_fractions
        )

        u_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

        trial_error = 0.0
        for time_idx in eachindex(sol.u)
            pred = u_named[time_idx].mass_fractions.methylene_blue[observed_cell_id]
            obs = trial.compared_data.mass_fractions.methylene_blue[time_idx]
            trial_error += sum(abs2, pred .- obs)
        end

        total_loss += trial_error
    end
    return total_loss
end

loss([log(1516.363624655201), 50000])

loss_history = [] #loss accumulator
parameter_history = [] #parameter accumulator

optimization_callback_tracker = function (state, l)
    println("loss: ", l)
    append!(loss_history, l)
    append!(parameter_history, [create_views_inline(state.u, system.u_proto_axes)])
end

N = ForwardDiff.pickchunksize(length(u0_vec))

adtype = Optimization.AutoForwardDiff(chunksize = N)
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

#diffusion pre-exponential factor, diffusion activation energy
guess_params = Float64[log(1516.363624655201), 50000]
loss(guess_params)

lower_bounds = [log(300.0), 1000.0]
upper_bounds = [log(1e12), 100000.0]

optprob = Optimization.OptimizationProblem(optf, guess_params, lb=lower_bounds, ub=upper_bounds)

@time res = Optimization.solve(
    optprob,
    #callback=cb,
    OptimizationOptimJL.LBFGS(),
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    f_abstol=1e-4,
    g_abstol=1e-4,
)

optimized_diffusion_pre_exponential_factor = exp(res.u[1])
optimized_diffusion_activation_energy = res.u[2]
#END OF OPTIMIZATION SOLVING

#START OF INITIAL SEARCHING

guess_lb = 1000
guess_ub = 1000000

Ea_range = collect(range(guess_lb, guess_ub, 20))

function loss_for_fixed_Ea(ln_A_val, fixed_Ea)
    return loss([ln_A_val[1], fixed_Ea])
end

function profile_Ea_scan(Ea_range)
    results = []

    for Ea_val in Ea_range
        prob_1D = Optimization.OptimizationProblem(
            Optimization.OptimizationFunction((x, p) -> loss_for_fixed_Ea(x, Ea_val), Optimization.AutoForwardDiff()),
            [log(1000.0)], # Initial guess for A
            lb = [log(100.0)], ub = [log(1e12)]
        )

        sol = Optimization.solve(prob_1D, OptimizationOptimJL.BFGS()) # Brent is optimal for 1D

        push!(results, (Ea=Ea_val, A=exp(sol.u[1]), Loss=sol.objective))
        println("Ea: $Ea_val, Opt_A: $(exp(sol.u[1])), Loss: $(sol.objective)")
    end
    return results
end

results = profile_Ea_scan(Ea_range)

A_range = [results[i].Ea for i in eachindex(results)]
loss_range = [results[i].Loss for i in eachindex(results)]

min_idx = argmin(loss_range)

results[min_idx]
results[min_idx].Ea
results[min_idx].A
results[min_idx].Loss

#FINAL PARAMETERS (DO NOT TOUCH)
#61.14478947368421
#650.6316641398093
#0.03102647047431776

#using Plots
plot(Ea_range, A_range)
plot(Ea_range, loss_range)
#END OF INITIAL SEARCHING

record_sol = true

sim_file = @__FILE__

u_proto_named = [(; create_views_inline(sol.u[i], system.u_proto_axes)...) for i in eachindex(sol.u)]

#u_named = rebuild_u_named(sol.u, u_proto_named)

cell = grid.cellsets["dialysis_tubing_interior"][1]

mass_fractions_dialysis_tubing_interior = []

for step in eachindex(sol.u)
    push!(mass_fractions_dialysis_tubing_interior, u_proto_named[step].mass_fractions.methylene_blue[cell])
end

mass_fractions_dialysis_tubing_interior

sol.u[1] == sol.u[end]

mass_fractions_beginning = u_proto_named[1].mass_fractions.methylene_blue[cell]

mass_fractions_end = u_proto_named[end].mass_fractions.methylene_blue[cell]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

if record_sol == true
    sol_to_vtk(sol, u_proto_named, grid, sim_file)
end