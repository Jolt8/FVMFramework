module FVMFramework

using Ferrite
using OrdinaryDiffEq
using LinearAlgebra
using SparseArrays
#using SciMLSensitivity
#using Optimization
#using OptimizationPolyalgorithms
#using Zygote
#using Enzyme
using ForwardDiff
using RecursiveArrayTools
#using OptimizationOptimJL
using ILUZero
using NonlinearSolve
using ComponentArrays
using StaticArrays
using ProfileView
using Dates
using WriteVTK
using PreallocationTools
using Polyester

import AlgebraicMultigrid
import SparseConnectivityTracer
import ADTypes
import Logging

#for future me, if you use this project on a new computer, in the FVMFramework environment just do:
#] activate .
#] instantiate

# ----- Geometry -----
include("geometry/geometry_helper_functions.jl")
export get_node_coordinates, get_cell_topology, get_nodes_of_cells, get_face_nodes
export get_cell_neighbors, get_unconnected_map, get_cell_face_map, cross_product

include("geometry/geometry_rebuilding_tetra.jl")
export rebuild_fvm_geometry_tetra!

include("geometry/geometry_rebuilding_hexa.jl")
export rebuild_fvm_geometry_hexa!

include("geometry/geometry_building.jl")
export build_fvm_geo_into_struct, FVMGeometry, FVMGeometryTetra, FVMGeometryHexa

include("geometry/topology_checking.jl")
export check_cellset_connectivity, check_grid_connectivity

# ----- Misc. Tools
include("named_tuple_interface/axes_generator.jl")
export create_axes

include("named_tuple_interface/inline_views.jl")
export create_views_inline

# ----- Physics ----
#   ---- Physics Types ----
include("physics_types/abstract_physics_types.jl")
export AbstractPhysics, AbstractFluidPhysics, AbstractSolidPhysics, AbstractReaction

#   ---- Flux Methods ----
#       --- Flux Physics ---
include("physics/flux_methods/flux_physics/advection.jl")
export species_advection!, all_species_advection!, enthalpy_advection!

include("physics/flux_methods/flux_physics/darcy_flow.jl")
export get_darcy_mass_flux, continuity_and_momentum_darcy

include("physics/flux_methods/flux_physics/diffusion.jl")
export species_numerical_flux, mass_fraction_diffusion!

include("physics/flux_methods/flux_physics/heat_transfer.jl")
export get_k_effective, numerical_flux, temp_diffusion!

#   ---- Internal Methods ----
#       --- Internal Capacities ---
include("physics/internal_methods/capacity_helper_functions.jl")
export cap_heat_flux_to_temp_change!, cap_mass_flux_to_pressure_change!, cap_species_mass_flux_to_mass_fraction_change!

#       --- Internal Physics ---
include("physics/internal_methods/internal_physics/chemistry.jl")
export power_law_react_cell!

include("physics/internal_methods/internal_physics/PAM_reforming.jl")
export PAM_reforming_react_cell!

#   ---- Helper Functions ----
include("physics/physics_helper_functions.jl")
export upwind, harmonic_mean #for fluxes
export R_gas #constant referenced almost everywhere
export van_t_hoff, arrenhius_k, K_gibbs_free #for chemical reactions
export mw_avg!, rho_ideal!, molar_concentrations!, get_cell_cp #other misc props

#   ---- Summation Functions ----
include("physics/physics_face_summation_functions.jl")
export sum_mass_flux_face_to_cell!

# ----- Setup and Recording Methods -----
#   ---- Sim Config ----
include("setup_and_recording/sim_config_helper_functions.jl")
export get_cell_set_total_volume, get_facet_set_total_area, get_facet_set_cells_respective_areas, get_cell_ids_in_facet_set

include("setup_and_recording/sim_config.jl")
export create_fvm_config, add_region!, add_patch!, add_controller!, finish_fvm_config
export SimulationConfigInfo, RegionSetupInfo, FVMSystem, AbstractController, ControllerSetupInfo


#   ---- Sim Recording ----
include("setup_and_recording/sim_recording.jl")
export sol_to_vtk

#   ---- Workflow Helper Functions ----
include("setup_and_recording/workflow_helper_functions.jl")
export rebuild_u_named, rebuild_u_named_vel, rebuild_u_named_mass_fractions

# ----- Solvers -----
#   ---- Preconditioners ----
include("solvers/preconditioners.jl")
export iluzero, algebraicmultigrid

#   ---- Callbacks ----
#       --- Progress Callbacks
include("solvers/callbacks/progress_callbacks.jl")
export show_t_progress, approximate_time_to_finish_cb

#   ---- FVM Operators ----
include("solvers/fvm_operators/methanol_reformer_op_different_connections.jl")
export methanol_reformer_f_test!

include("solvers/fvm_operators/methanol_reformer_op_steady_state.jl")
export methanol_reformer_f_steady_state!

include("solvers/fvm_operators/parameter_fitting_operator.jl")
export methylene_blue_diffuion_parameter_fitting_f!

include("solvers/fvm_operators/heat_transfer_minimal_allocs.jl")
export heat_transfer_f_test!

include("solvers/fvm_operators/pipe_f!.jl")
export pipe_f!

include("solvers/shared_group_functions/shared_group_functions.jl")
export solve_connection_group!, solve_controller_group!, solve_region_group!, solve_patch_group!

#   ---- Solver Debugging ----
include("solvers/solver_debugging.jl")
export debug_region!
end
