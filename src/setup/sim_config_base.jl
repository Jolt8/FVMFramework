mutable struct RegionSetupInfo{P <: AbstractPhysics} #this must be defined before SimulationConfigInfo
    name::String
    type::P
    initial_conditions::ComponentVector
    properties::ComponentVector
    cache_syms_and_units::NamedTuple
    region_function::Function
    region_cells::Vector{Int}
end

mutable struct PatchSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    properties::ComponentVector
    patch_function::Function
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}}
end

mutable struct ControllerSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    controller::ComponentVector
    monitored_cellset::String
    affected_cellset::String
    controller_function::Function
    monitored_cells::Vector{Int}
    affected_cells::Vector{Int}
end

struct SimulationConfigInfo
    grid::Ferrite.Grid
    geo::FVMGeometry
    top::ExclusiveTopology
    regions::Vector{RegionSetupInfo}
    patches::Vector{PatchSetupInfo}
    controllers::Vector{ControllerSetupInfo}
    optimized_parameters::Dict{Symbol, Any}
    u_proto::ComponentVector
end

function create_fvm_config(grid, u_proto)
    top = ExclusiveTopology(grid)

    geo = build_fvm_geo_into_struct(grid, top)

    return SimulationConfigInfo(
        grid,
        geo,
        top, 
        RegionSetupInfo[],
        PatchSetupInfo[],
        ControllerSetupInfo[],
        Dict{Symbol, Any}(), 
        u_proto
    )
end