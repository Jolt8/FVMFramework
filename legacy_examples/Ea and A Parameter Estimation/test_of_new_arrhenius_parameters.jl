

using Ferrite
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using SciMLSensitivity
using Optimization, OptimizationPolyalgorithms, Zygote
#using Enzyme
using RecursiveArrayTools
using OptimizationOptimJL
using ILUZero
import AlgebraicMultigrid
#import SparseConnectivityTracer, ADTypes
using NonlinearSolve
import Logging
using ComponentArrays
using StaticArrays
using ProfileView
using FVMFramework
using Unitful

abstract type AbstractPhysics end

struct ChemicalReaction
    heat_of_reaction::Float64
    delta_gibbs_free_energy::Float64
    K_gibbs_free_ref_temp::Float64
    reactants::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    reactant_stoich_coeffs::Vector{Int} #(stoich_coeff_for_reactant_1, stoich_coeff_for_reactant_2, etc...)
    products::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    product_stoich_coeffs::Vector{Int} #[stoich_coeff_for_product_1, stoich_coeff_for_product_2, etc...]
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoich_for_chemical_id_2, etc...)
    #=
    commenting these out because there are our optimized_parameters

    kf_A # Pre-exponential factor for forward reaction
    kf_Ea # Activation energy for forward reaction
    kr_A # Pre-exponential factor for reverse reaction
    kr_Ea  # Activation energy for reverse reaction
    =#
end

struct ChemPhysics <: AbstractPhysics
    k::Float64
    rho::Float64
    cp::Float64
    chemical_reactions::Vector{ChemicalReaction}
end

abstract type AbstractBC end

struct HeatBC <: AbstractBC
    initial_temp::Float64
end

struct ChemBC <: AbstractBC
    initial_mass_fractions::Vector{Float64}
end

struct MultiPhysicsBCs
    chem_bcs::Vector{ChemBC}
    temp_bcs::Vector{HeatBC}
end

struct BoundarySystem
    boundary_map::MultiPhysicsBCs
    free_idxs::Vector{Int}
    dirichlet_idxs::Vector{Int}
end

function numerical_flux(k_avg, T_L, T_R, area, dist)
    grad_T = (T_R - T_L) / dist
    q = -k_avg * grad_T
    return q * area
end

function net_reaction_rate(chemical_reaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    R_gas_kj = 0.008314
    kf = (kf_A * exp(-kf_Ea / (R_gas_kj * T)))
    kr = (kr_A * exp(-kr_Ea / (R_gas_kj * T)))

    forward_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.reactants)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.reactant_stoich_coeffs[i]
        forward_term *= concentration^stoich_coeff
    end

    reverse_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.products)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.product_stoich_coeffs[i]
        reverse_term *= concentration^stoich_coeff
    end

    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))

    return net_reaction_rate
end

function K_gibbs_free_2(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    K_ref = exp(-ΔG_rxn_ref / (8.314e-3 * T_ref)) #R is in kJ

    ln_K_ratio = (-ΔH_rxn_ref / 8.314e-3) * (1 / T_actual - 1 / T_ref)

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

function FVM_iter_f!(
    du, u, p, t, timestamps, temperatures_vec,
    cell_volumes,
    cell_centroids, connection_areas, connection_normals,
    #cell volumes and cell centroids are accessed at the id of the cell
    connection_distances, unconnected_areas,
    #connection areas, normals, and distances are  accessed by their location in the 
    #list which corresponds to the respective connection in cell_neighbor_map
    species_molecular_weights,
    cell_props_id_map, bc_sys::BoundarySystem, chem_phys::Vector{ChemPhysics}, ax, n_reactions, n_species
)

    A_Ea_pairs = eachcol(reshape(p, :, n_reactions))
    #unflattened_p would be [[reaction_1_kf_A, reaction_1_kf_Ea], [reaction_2_kf_A, reaction_2_kf_Ea], etc..] 

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)

    du .= 0.0

    #=
    for (i, (idx_a, idx_b)) in enumerate(cell_neighbor_map)
        #diffusion here later
    end
    =#

    # Source and Capacity Loop

    time_idx = argmin(abs.(timestamps .- t))

    for cell_id in bc_sys.free_idxs
        vol = cell_volumes[cell_id]
        props = cell_props_id_map[cell_id]

        k = chem_phys[props].k
        rho = chem_phys[props].rho
        cp = chem_phys[props].cp
        cell_chemical_reactions_vec = chem_phys[props].chemical_reactions

        #S = chem_phys[props].source_term * vol 
        # we should probably create separate containers in chem_phys for both source terms on a per area and per cell basis

        species_mass_fractions = view(u.mass_fractions, :, cell_id)

        species_molar_concentrations = [
            (rho * species_mass_fraction) / species_molecular_weights[species_id]
            for (species_id, species_mass_fraction) in enumerate(species_mass_fractions)
        ]
        #this above causes a lot of GC

        for (reaction_id, reaction) in enumerate(cell_chemical_reactions_vec)
            kf_A = A_Ea_pairs[reaction_id][1]
            kf_Ea = A_Ea_pairs[reaction_id][2]

            #find reverse pre exponential_factor
            K_ref = K_gibbs_free_2(reaction.K_gibbs_free_ref_temp, temperatures_vec[time_idx, cell_id], reaction.delta_gibbs_free_energy, reaction.heat_of_reaction)

            kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (8.314e-3 * temperatures_vec[time_idx, cell_id]))

            #find reverse Ea
            kr_Ea = kf_Ea - reaction.heat_of_reaction

            net_rates = [net_reaction_rate(reaction, species_molar_concentrations, temperatures_vec[time_idx, cell_id], kf_A, kf_Ea, kr_A, kr_Ea) for reaction in cell_chemical_reactions_vec]

            for (species_id, species_molar_concentration) in enumerate(species_molar_concentrations)
                change_in_species_molar_concentration = 0.0

                for (reaction_idx, reaction) in enumerate(cell_chemical_reactions_vec)
                    stoich = reaction.all_stoich_coeffs[species_id]
                    change_in_species_molar_concentration += net_rates[reaction_idx] * stoich
                end

                du.mass_fractions[species_id, cell_id] = (change_in_species_molar_concentration * species_molecular_weights[species_id]) / rho
            end
        end
    end
end

grid_dimensions = (1, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "internal_cells", x -> x != "chicken") # all

# acetic acid, ethanol, water, ethyl acetate
initial_mass_fractions = [0.5, 0.5, 0.0, 0.0]

acetic_acid_ethanol_esterification_rxn = ChemicalReaction(
    -4.0, #Delta H
    -4.0, #Delta Gibbs free at ref temp
    298, #ref temp
    [1, 2], #reactant_ids
    [1, 1], #reactant_stoich_coeffs
    [3, 4], #product_ids
    [1, 1], #product_stoich_coeffs
    [-1, -1, 1, 1] #stoich coefficients -1 = reactant, 1 = product
)

reaction_physics = ChemPhysics(0.6e-3, 1000, 4.184, [acetic_acid_ethanol_esterification_rxn])

chem_phys_vec = [reaction_physics]

internal_cell_set_idxs = Set(getcellset(grid, "internal_cells"))

struct CellSet
    props_id::Int
    cell_set_idxs::Set{Int}
end

internal_cell_set = CellSet(1, internal_cell_set_idxs)

cell_sets = [internal_cell_set]

free_idxs = Int[]
dirichlet_idxs = Int[]
#should probably create separate vectors for situations where we need to fix temperature but not mass fractions

function my_bc_mapper(cell_id)
    if cell_id in cell_sets[1].cell_set_idxs
        chem_bc = ChemBC(initial_mass_fractions)
        heat_bc = HeatBC(500.0)
        push!(free_idxs, cell_id)
        return [chem_bc, heat_bc]
    end
end

chem_bcs = ChemBC[]
heat_bcs = HeatBC[]

n_cells = length(grid.cells)

for cell_id in 1:n_cells
    bcs = my_bc_mapper(cell_id) #returns vector [BC1, BC2]
    push!(chem_bcs, bcs[1])
    push!(heat_bcs, bcs[2])
end

boundary_map = MultiPhysicsBCs(chem_bcs, heat_bcs)

bc_sys = BoundarySystem(boundary_map, free_idxs, dirichlet_idxs)

n_species = length(initial_mass_fractions)
alloc_mass_fraction_vec = zeros(n_species, n_cells)

u_proto = ComponentArray(mass_fractions=alloc_mass_fraction_vec)

for cell_id in 1:n_cells
    u_proto.mass_fractions[:, cell_id] = chem_bcs[cell_id].initial_mass_fractions
end

cell_props_id_map = Int[]

for cell in CellIterator(grid)
    cell_id = cellid(cell)
    for i in eachindex(cell_sets)
        if cell_id in cell_sets[i].cell_set_idxs
            push!(cell_props_id_map, cell_sets[i].props_id)
        end
    end
end

initial_node_coordinates = get_node_coordinates(grid)

cell_neighbor_map = nothing
neighbor_map_respective_node_ids = nothing
unconnected_cell_face_map = nothing
unconnected_map_respective_node_ids = nothing

nodes_of_cells = get_nodes_of_cells(grid)

cell_volumes = [50u"ml" |> u"m^3"]
cell_centroids = nothing
connection_areas = nothing
connection_normals = nothing
connection_distances = nothing
unconnected_areas = nothing
unconnected_normals = nothing

#= 
NOTE: We might have to declare these as constant to prevent dynamic dispatch in the future but we're leaving this out for now
const cell_neighbor_map = cell_neighbor_map 
const const_neighbor_map_respective_node_ids = neighbor_map_respective_node_ids
=#

n_reactions = length(reaction_physics.chemical_reactions)
n_species = length(initial_mass_fractions)

u0 = Vector(u_proto)

du0 = u0 .* 0.0

#Experimental data retrieval and processing

using XLSX

experimental_data_path = joinpath(@__DIR__, "esterification_data_processing_for_julia.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

struct TrialPreprocessedData
    temperature_timestamps::Vector{Float64}
    temperatures::Vector{Float64}
    moles_timestamps::Vector{Float64}
    moles_respective_temp_idxs::Vector{Int64}
    acetic_acid_moles::Vector{Float64}
    ethanol_moles::Vector{Float64}
    ethyl_acetate_moles::Vector{Float64}
    water_moles::Vector{Float64}
end

T1_temp_end_idx = xf["30C Temps"]["I3"]
T1_temp_timestamps = float.(vec(xf["30C Temps"]["K2:K$T1_temp_end_idx"]))
T1_moles_timestamps = float.(vec(xf["30C"]["A3:A11"]))
T1_moles_respective_temp_idxs = findall(x -> x in T1_moles_timestamps, T1_temp_timestamps)
trial_1_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["30C Temps"]["K2:K$T1_temp_end_idx"])),
    #temperatures
    float.(vec(xf["30C Temps"]["L2:L$T1_temp_end_idx"])),

    #moles_timestamps 
    float.(vec(xf["30C"]["A3:A11"])),
    #moles_respective_temp_idxs
    T1_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["30C"]["G3:G11"])),
    #ethanol_moles
    float.(vec(xf["30C"]["H3:H11"])),
    #ethyl_acetate_moles
    float.(vec(xf["30C"]["J3:J11"])),
    #water_moles
    float.(vec(xf["30C"]["I3:I11"])),
)


T2_temp_end_idx = xf["40C Temps"]["I3"]
T2_temp_timestamps = float.(vec(xf["40C Temps"]["K2:K$T2_temp_end_idx"]))
T2_moles_timestamps = float.(vec(xf["40C"]["A3:A11"]))
T2_moles_respective_temp_idxs = findall(x -> x in T2_moles_timestamps, T2_temp_timestamps)
trial_2_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["40C Temps"]["K2:K$T2_temp_end_idx"])),
    #temperatures
    float.(vec(xf["40C Temps"]["L2:L$T2_temp_end_idx"])),

    #moles_timestamps 
    float.(vec(xf["40C"]["A3:A11"])),
    #moles_respective_temp_idxs
    T2_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["40C"]["E3:E11"])),
    #ethanol_moles
    float.(vec(xf["40C"]["F3:F11"])),
    #ethyl_acetate_moles
    float.(vec(xf["40C"]["H3:H11"])),
    #water_moles
    float.(vec(xf["40C"]["G3:G11"])),
)

T3_temp_end_idx = xf["50C Temps"]["I3"]
T3_temp_timestamps = float.(vec(xf["50C Temps"]["K2:K$T3_temp_end_idx"]))
T3_moles_timestamps = float.(vec(xf["50C"]["A3:A14"]))
T3_moles_respective_temp_idxs = findall(x -> x in T3_moles_timestamps, T3_temp_timestamps)
trial_3_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["50C Temps"]["K2:K$T3_temp_end_idx"])),
    #temperatures
    float.(vec(xf["50C Temps"]["L2:L$T3_temp_end_idx"])),

    #moles_timestamps
    float.(vec(xf["50C"]["A3:A14"])),
    #moles_respective_temp_idxs
    T3_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["50C"]["G3:G14"])),
    #ethanol_moles
    float.(vec(xf["50C"]["H3:H14"])),
    #ethyl_acetate_moles
    float.(vec(xf["50C"]["J3:J14"])),
    #water_moles
    float.(vec(xf["50C"]["I3:I14"])),
)

all_trials_data = [trial_1_data, trial_2_data, trial_3_data]

struct Trial
    temperature_timestamps::Vector{Float64}
    temperatures_vec::Matrix{Float64}

    moles_timestamps::Vector{Float64}
    moles_respective_temp_idxs::Vector{Int64}
    mass_fractions_matrix::Array{Float64, 3} #[timestamp, reactant_idx]
end

function get_mass_fraction(species_moles, species_molecular_weight, rho, mixture_volume)
    (((species_moles / mixture_volume) * species_molecular_weight) / rho)
end

trials = Trial[]

for trial in all_trials_data
    temp_timestamps = trial.temperature_timestamps .* u"s"

    temperatures_deg_C = trial.temperatures .* u"°C"

    trial_temperatures = temperatures_deg_C .|> u"K"

    moles_timestamps = trial.moles_timestamps
    moles_respective_temp_idxs = trial.moles_respective_temp_idxs

    acetic_acid_moles = trial.acetic_acid_moles .* u"mol"
    ethanol_moles = trial.ethanol_moles .* u"mol"
    ethyl_acetate_moles = trial.ethyl_acetate_moles .* u"mol"
    water_moles = trial.water_moles .* u"mol"

    mixture_volume = 50u"ml"
    mixture_rho = 1000.0u"kg/m^3"

    acetic_acid_mw = 60.05u"g/mol"
    ethanol_mw = 46.07u"g/mol"
    ethyl_acetate_mw = 88.11u"g/mol"
    water_mw = 18.02u"g/mol"

    acetic_acid_mass_fractions = get_mass_fraction.(acetic_acid_moles, acetic_acid_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
    ethanol_mass_fractions = get_mass_fraction.(ethanol_moles, ethanol_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
    ethyl_acetate_mass_fractions = get_mass_fraction.(ethyl_acetate_moles, ethyl_acetate_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
    water_mass_fractions = get_mass_fraction.(water_moles, water_mw, mixture_rho, mixture_volume) .|> u"kg/kg"

    mass_fractions_matrix = zeros(length(moles_timestamps), n_species, n_cells)
    temperature_vec = zeros(length(temp_timestamps), n_cells)

    for i in eachindex(moles_timestamps)
        mass_fractions_matrix[i, :, :] .= [acetic_acid_mass_fractions[i], ethanol_mass_fractions[i], ethyl_acetate_mass_fractions[i], water_mass_fractions[i]]
    end
    #mass_fractions_matrix is formatted as [time_idx, species_idx, cell_idx]

    for i in eachindex(temp_timestamps)
        temperature_vec[i, :] .= ustrip(trial_temperatures[i])
    end
    #temperature_vec is formatted as [time_idx, cell_idx]

    trial = Trial(ustrip.(temp_timestamps), ustrip.(temperature_vec), ustrip.(moles_timestamps), moles_respective_temp_idxs, mass_fractions_matrix)

    push!(trials, trial)
end

#trials[1]
#trials[2]
#trials[3]

errors = zeros(length(trials))

#=
function update_temp!(integrator)
    integrator.u[4] = (integrator.u[4] * m_initial + m_ice) / m_total
end

update_temp_callback = PresetTimeCallback(300.0, update_temp!)
=#

u_axes = getaxes(u_proto)[1]
#we use u_proto and u_axes in the function because KrylovJL_GMRES complains when a component array from ComponentArrays.jl is passed in

u0 = Vector(u_proto)

reaction_1_kf_A_guess = 38394.91787122687 / 1000
reaction_1_kf_Ea_guess = 88.80051948962193 / 1000

proto_p_guess = [reaction_1_kf_A_guess, reaction_1_kf_Ea_guess] #27 element Vector{Vector{}

p_guess = reduce(vcat, proto_p_guess) #81 element Vector{} # we do this because Optimization.solve(...) complains about vectors of vectors

mc_proto_length = n_cells * n_species

f_closure_implicit_pre = (du, u, p, t, temperature_timestamps, temperature_vec) -> FVM_iter_f!(
    du, u, p, t, temperature_timestamps, temperature_vec,
    cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances,
    unconnected_areas,
    [0.06005, 0.04607, 0.08811, 0.018015],
    cell_props_id_map, bc_sys, chem_phys_vec, u_axes, n_reactions, n_species
)

using SparseConnectivityTracer

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end

trials_mass_fractions_results = []
trials_timestamps = []

u0_for_jac = trials[1].mass_fractions_matrix[1, :, 1]
du0_for_jac = u0_for_jac .* 0.0

f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trials[1].temperature_timestamps, trials[1].temperatures_vec)

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_for_jac, u0_for_jac, detector)

for trial in trials
    trial_u0 = trial.mass_fractions_matrix[1, :, 1]
    trial_du0 = trial_u0 .* 0.0

    f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.temperature_timestamps, trial.temperatures_vec)

    ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

    trial_t0, trial_tMax = trial.temperature_timestamps[1], trial.temperature_timestamps[end]
    tspan = (trial_t0, trial_tMax)
    implicit_prob = ODEProblem(ode_func, trial_u0, tspan, p_guess)

    @time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), saveat=trial.moles_timestamps)

    push!(trials_mass_fractions_results, sol.u)
    push!(trials_timestamps, sol.t)
end

println(trials_mass_fractions_results[1])
println(trials_timestamps[1])
println(trials_mass_fractions_results[2])
println(trials_timestamps[2])
println(trials_mass_fractions_results[3])
println(trials_timestamps[3])

using DataFrames
ethanol_mass_fractions = []
acetic_acid_mass_fractions = []
ethyl_acetate_mass_fractions = []
water_mass_fractions = []

for mass_fractions in trials_mass_fractions_results[1]
    push!(ethanol_mass_fractions, mass_fractions[1])
    push!(acetic_acid_mass_fractions, mass_fractions[2])
    push!(ethyl_acetate_mass_fractions, mass_fractions[3])
    push!(water_mass_fractions, mass_fractions[4])
end

df1 = DataFrame(AA=trials_timestamps[1], AB=ethanol_mass_fractions, AC=acetic_acid_mass_fractions, AD=ethyl_acetate_mass_fractions, AE=water_mass_fractions)

ethanol_mass_fractions = []
acetic_acid_mass_fractions = []
ethyl_acetate_mass_fractions = []
water_mass_fractions = []
for mass_fractions in trials_mass_fractions_results[2]
    push!(ethanol_mass_fractions, mass_fractions[1])
    push!(acetic_acid_mass_fractions, mass_fractions[2])
    push!(ethyl_acetate_mass_fractions, mass_fractions[3])
    push!(water_mass_fractions, mass_fractions[4])
end

df2 = DataFrame(AA=trials_timestamps[2], AB=ethanol_mass_fractions, AC=acetic_acid_mass_fractions, AD=ethyl_acetate_mass_fractions, AE=water_mass_fractions)

ethanol_mass_fractions = []
acetic_acid_mass_fractions = []
ethyl_acetate_mass_fractions = []
water_mass_fractions = []
for mass_fractions in trials_mass_fractions_results[3]
    push!(ethanol_mass_fractions, mass_fractions[1])
    push!(acetic_acid_mass_fractions, mass_fractions[2])
    push!(ethyl_acetate_mass_fractions, mass_fractions[3])
    push!(water_mass_fractions, mass_fractions[4])
end
df3 = DataFrame(AA=trials_timestamps[3], AB=ethanol_mass_fractions, AC=acetic_acid_mass_fractions, AD=ethyl_acetate_mass_fractions, AE=water_mass_fractions)


XLSX.writetable("arrhenius_test_results_2.xlsx", 
    "T1" => df1, 
    "T2" => df2,
    "T3" => df3
)

println("done")