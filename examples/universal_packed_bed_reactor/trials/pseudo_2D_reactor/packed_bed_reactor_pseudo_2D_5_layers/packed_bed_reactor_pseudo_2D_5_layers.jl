using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak

using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using ForwardDiff
using DataFrames
using CSV
using DataInterpolations
using Dates

using FVMFramework

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "multiphase.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "energy.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "momentum.jl"))
Revise.includet(joinpath(@__DIR__, "..", "packing_and_fluid_property_merging.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "convective_heat_transfer.jl"))

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "5_layer_common_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "packed_media_properties", "copper_mesh_reformer_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "packed_media_properties", "silicon_carbide_sand_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "steel_pipe_wall_properties", "steel_pipe_wall_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "thermocouple_and_heating_wire_properties", "divided_thermocouple_and_heating_wire_properties", "thermocouple_and_jacket_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "thermocouple_and_heating_wire_properties", "divided_thermocouple_and_heating_wire_properties", "heating_wire_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "insulation_properties", "insulation_properties.jl"))
Revise.includet(joinpath(@__DIR__, "pseudo_2D_geometry_editing_5_layers.jl"))

Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_dry", "empty_trial_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_dry", "air_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "thermocouple_data_processing", "thermocouple_data_with_wattage.jl"))

Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "hot_water_trial_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "water_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "recorded_data_processing", "inlet_and_outlet_temperatures.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "recorded_data_processing", "values_of_note.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "thermocouple_data_processing", "thermocouple_data.jl"))

include("src/physics.jl")

include("src/setup_fvm.jl")

include("src/dry_run_trial.jl")
include("src/hot_water_trial.jl")

include("src/ode_system.jl")

include("src/evaluation.jl")
dry_run_config, dry_run_properties, dry_run_thermocouple_data, 
dry_run_heater_1_wattage_per_cell, dry_run_heater_2_wattage_per_cell, dry_run_heater_3_wattage_per_cell, 
dry_run_heater_4_wattage_per_cell, dry_run_heater_5_wattage_per_cell = dry_run_trial();
dry_run_du0_vec, dry_run_u0_vec, dry_run_state_axes, dry_run_geo, dry_run_system = finish_fvm_config(dry_run_config, connection_map_function, check_units = false);

hot_water_config, hot_water_properties, hot_water_pump_shut_off, hot_water_thermocouple_data, 
hot_water_heater_1_wattage_per_cell, hot_water_heater_2_wattage_per_cell, hot_water_heater_3_wattage_per_cell, 
hot_water_heater_4_wattage_per_cell, hot_water_heater_5_wattage_per_cell = hot_water_trial();
hot_water_du0_vec, hot_water_u0_vec, hot_water_state_axes, hot_water_geo, hot_water_system = finish_fvm_config(hot_water_config, connection_map_function, check_units = false);

p_axes = hot_water_system.p_axes

f_closure_dry_run = (du, u, p, t) -> fvm_operator!(du, u, p, t, dry_run_solve_system!, dry_run_geo, dry_run_system)
f_closure_hot_water = (du, u, p, t) -> fvm_operator!(du, u, p, t, hot_water_solve_system!, hot_water_geo, hot_water_system)

first_p_guess_init = ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 3.9215249544953634u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.21843778568123667u"W/K",
    heater_weight_1 = 0.04233643564131222,
    heater_weight_2 = 0.6599848372698764,
    heater_weight_3 = 0.4480116180759421,
    heater_weight_4 = 0.007962092313654773,
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 61.49737866475787u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 3.7753230723789706,
    TC1_thermal_resistance = 28.129313413301414u"K/W",
    TC2_thermal_resistance = 97.1810904949918u"K/W",
    TC3_thermal_resistance = 91.60608551471157u"K/W",
    TC4_thermal_resistance = 57.90696805772799u"K/W",
    TC5_thermal_resistance = 67.15056955503321u"K/W",
)

dry_run_best_params = ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 1.3613733508973453u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.046224916624426536u"W/K",
    heater_weight_1 = 0.04233643564131222,
    heater_weight_2 = 0.6599848372698764,
    heater_weight_3 = 0.4480116180759421,
    heater_weight_4 = 0.007962092313654773,
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 1194.6515661535152u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 2.1703947052768693,
    TC1_thermal_resistance = 28.129313413301414u"K/W",
    TC2_thermal_resistance = 97.1810904949918u"K/W",
    TC3_thermal_resistance = 91.60608551471157u"K/W",
    TC4_thermal_resistance = 57.90696805772799u"K/W",
    TC5_thermal_resistance = 67.15056955503321u"K/W",
)

hot_water_best_params = ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 3.9215249544953634u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.21843778568123667u"W/K",
    heater_weight_1 = 0.3211331590638462,
    heater_weight_2 = 0.14535035339598748,
    heater_weight_3 = 0.6840217991331898,
    heater_weight_4 = 0.9984762651258827,
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 61.49737866475787u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 3.7753230723789706,
    TC1_thermal_resistance = 31.29201644638937u"K/W",
    TC2_thermal_resistance = 2.204398102203747u"K/W",
    TC3_thermal_resistance = 4.116298695255978u"K/W",
    TC4_thermal_resistance = 89.06507428995272u"K/W",
    TC5_thermal_resistance = 185.05463278797103u"K/W"
)

p_axes = getaxes(first_p_guess_init)
p_guess = ustrip.(upreferred.(Vector(first_p_guess_init)))

dry_run_p_best = ustrip.(upreferred.(Vector(dry_run_best_params)))
hot_water_p_best = ustrip.(upreferred.(Vector(hot_water_best_params)))

dry_run_implicit_prob = built_trial_implicit_prob(f_closure_dry_run, dry_run_du0_vec, dry_run_u0_vec, dry_run_thermocouple_data, dry_run_p_best);
hot_water_implicit_prob = built_trial_implicit_prob(f_closure_hot_water, hot_water_du0_vec, hot_water_u0_vec, hot_water_thermocouple_data, hot_water_p_best);

desired_steps = 100
dry_run_saveat = ustrip(upreferred(dry_run_thermocouple_data.timestamps[end])) / desired_steps
hot_water_saveat = ustrip(upreferred(hot_water_thermocouple_data.timestamps[end])) / desired_steps

# Executing test solves (uncomment if you want to run single simulations in the REPL)
# @time dry_run_test_sol = solve(dry_run_implicit_prob, FBDF(linsolve = SparspakFactorization()), dtmax = dry_run_saveat, callback = approximate_time_to_finish_cb)
# @time hot_water_test_sol = solve(hot_water_implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)

# dry_run_u_named = [ComponentVector(dry_run_test_sol.u[i], dry_run_state_axes) for i in 1:length(dry_run_test_sol.u)];
# hot_water_u_named = [ComponentVector(hot_water_test_sol.u[i], hot_water_state_axes) for i in 1:length(hot_water_test_sol.u)];

sim_file = @__FILE__
# dry_run_root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files\\dry_run_trial_output"
# hot_water_root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files\\hot_water_flow_trial_output"
# sol_to_vtk(dry_run_test_sol, dry_run_u_named, grid, sim_file, dry_run_root_dir)
# sol_to_vtk(hot_water_test_sol, hot_water_u_named, grid, sim_file, hot_water_root_dir)

#generate_validation_plots(dry_run_p_best, hot_water_p_best)

dry_run_timestamps = ustrip.(dry_run_thermocouple_data.timestamps)
hot_water_timestamps = ustrip.(hot_water_thermocouple_data.timestamps)

TC1_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC1_closest_cell_id))
TC2_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC2_closest_cell_id))
TC3_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC3_closest_cell_id))
TC4_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC4_closest_cell_id))
TC5_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC5_closest_cell_id))

n_saves = 100
dry_run_save_interval = ustrip(upreferred(dry_run_thermocouple_data.timestamps[end])) / n_saves
hot_water_save_interval = ustrip(upreferred(hot_water_thermocouple_data.timestamps[end])) / n_saves

dry_run_heater_powers = dry_run_thermocouple_data.heater_power_interp.(ustrip.(dry_run_thermocouple_data.timestamps))
dry_run_tstops = ustrip.(upreferred.(dry_run_thermocouple_data.timestamps[findall(!iszero, dry_run_heater_powers)]))

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> pure_dry_run_loss(x), adtype)

p_lower_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 0.01u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.01u"W/K",
    heater_weight_1 = 0.05u"1",
    heater_weight_2 = 0.12u"1",
    heater_weight_3 = 0.20u"1",
    heater_weight_4 = 0.30u"1",
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 1.0u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 0.1,
    TC1_thermal_resistance = 0.1u"K/W",
    TC2_thermal_resistance = 0.1u"K/W",
    TC3_thermal_resistance = 0.1u"K/W",
    TC4_thermal_resistance = 0.1u"K/W",
    TC5_thermal_resistance = 0.1u"K/W",
))))

p_upper_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 100.0u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 5.0u"W/K",
    heater_weight_1 = 0.3u"1",
    heater_weight_2 = 0.35u"1",
    heater_weight_3 = 0.50u"1",
    heater_weight_4 = 0.80u"1",
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 2000.0u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 10.0,
    TC1_thermal_resistance = 200.0u"K/W",
    TC2_thermal_resistance = 200.0u"K/W",
    TC3_thermal_resistance = 200.0u"K/W",
    TC4_thermal_resistance = 200.0u"K/W",
    TC5_thermal_resistance = 200.0u"K/W",
))))

optprob = Optimization.OptimizationProblem(optf, p_guess, lb = p_lower_bounds, ub = p_upper_bounds)

function randomize(lower, upper)
    return lower + (upper - lower) * rand()
end

p_ensemble = [[randomize(p_lower_bounds[i], p_upper_bounds[i]) for i in eachindex(p_lower_bounds)] for _ in 1:Sys.CPU_THREADS]

function prob_func(prob, i, repeat)
    return remake(prob, u0 = p_ensemble[i])
end

ensembleprob = EnsembleProblem(optprob; prob_func)

LOSS = Float64[]
PARS = []

mkpath(joinpath(@__DIR__, "optimization_results"))
results_path = joinpath(@__DIR__, "optimization_results", "optimization_results_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).csv")

open(results_path, "w") do io
    header_str = "loss," * join(string.(propertynames(first_p_guess_init)), ",")
    println(io, header_str)
end

const cb_lock = ReentrantLock()

cb = function (state, l)
    display(l)
    display(state.u)
    
    lock(cb_lock) do
        push!(LOSS, l)
        push!(PARS, state.u)
        
        p_named = NamedTuple(ComponentVector(state.u, p_axes))
        row = merge((loss = l, ), p_named)
        
        CSV.write(results_path, DataFrame([row]), append=true)
    end
    
    false
end

#=
@time res = Optimization.solve(
    ensembleprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    EnsembleThreads(),
    trajectories = 50,
    #Sys.CPU_THREADS,
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    #f_abstol=1e-8,
    #g_abstol=1e-8,
)=#

#=
@time res = Optimization.solve(
    optprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    #EnsembleThreads(),
    #trajectories = 50,
    #Sys.CPU_THREADS,
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    #f_abstol=1e-8,
    #g_abstol=1e-8,
)=#

@time res = solve(
    ensembleprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    EnsembleThreads(),
    trajectories = Sys.CPU_THREADS,
    #PopulationSize = 100,
    #Method = :RandomSearcher,
    callback = cb,
)

#=
@time res = solve(
    optprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    callback = cb,
    PoulationSize = 100,
    Method = :RandomSearcher
    #Method = :SepReal
)
    =#

#=
dry_run_optimization_results_path = joinpath(@__DIR__, "optimization_results", "dry_run_optimization_results.csv")
dry_run_df = CSV.read(dry_run_optimization_results_path, DataFrame)

min_dry_run_loss_idx = argmin(dry_run_df.loss)
min_dry_run_loss = dry_run_df.loss[min_dry_run_loss_idx]
dry_run_best_params = Vector(dry_run_df[min_dry_run_loss_idx, 2:end])

hot_water_optimization_results_path = joinpath(@__DIR__, "optimization_results", "hot_water_optimization_results.csv")
hot_water_df = CSV.read(hot_water_optimization_results_path, DataFrame)

min_hot_water_loss_idx = argmin(hot_water_df.loss)
min_hot_water_loss = hot_water_df.loss[min_hot_water_loss_idx]
hot_water_best_params = Vector(hot_water_df[min_hot_water_loss_idx, 2:end])
=#