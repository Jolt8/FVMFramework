struct MethanolReformerPhysics <: AbstractFluidPhysics
    k::Float64
    cp::Float64
    mu::Float64
    permeability::Float64
    species_diffusion_coeffs::Vector{Float64}
    species_molecular_weights::Vector{Float64}
    chemical_reactions::Vector{AbstractReaction}
    cell_kg_cat_per_m3_for_each_reaction::Vector{Float64}
end

struct WallPhysics <: AbstractSolidPhysics
    k::Float64
    rho::Float64
    cp::Float64
end

function solve_connection_group!(
        du, u, phys_a::TA, phys_b::TB, flux!, neighbors,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    ) where {TA, TB}
    #look into storing this data into a CSR in the future
    @batch for (idx_a, neighbor_list) in neighbors
        for (idx_b, face_idx) in neighbor_list
            flux!(
                du, u, phys_a, phys_b,
                idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                rho_cache, mw_avg_cache,
                change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
            )
        end
    end
end

function solve_controller_group!(
        du, u, controller::T, controller_id, control!, monitored_cells, affected_cells, 
        cell_volumes
    ) where T
    #no loop here because control! already loops over monitored_cells and affected_cells
    control!(
        du, u, controller, controller_id, monitored_cells, affected_cells,
        cell_volumes
    )
end

function solve_region_group!(
        du, u, phys::T, internal_physics!, region_cells,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    ) where T
    @batch for cell_id in region_cells
        internal_physics!(
            du, u, phys, cell_id,
            cell_volumes[cell_id],
            rho_cache, mw_avg_cache,
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
        )
    end
end


function methanol_reformer_f_test!(
        du, u, p, t,
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        unconnected_cell_face_map, cell_face_areas, cell_face_normals,
        connection_groups, controller_groups, region_groups,
        #phys::Vector{AbstractPhysics}, cell_phys_id_map::Vector{Int}, #we probably won't use these again, but they might be useful in the future
        u_merged_axes,
        rho_cache, mw_avg_cache,
        change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
    )

    u = ComponentVector(u, u_merged_axes)
    du = ComponentVector(du, u_merged_axes)
    du .= 0.0

    rho_cache = get_tmp(rho_cache, u)
    rho_cache .= 0.0

    mw_avg_cache = get_tmp(mw_avg_cache, u)
    mw_avg_cache .= 0.0

    change_in_molar_concentrations_cache = get_tmp(change_in_molar_concentrations_cache, u)
    change_in_molar_concentrations_cache .= 0.0

    molar_concentrations_cache = get_tmp(molar_concentrations_cache, u) #just using mass fractions for cell 1, this may cause some issues later!
    molar_concentrations_cache .= 0.0

    net_rates_cache = get_tmp(net_rates_cache, u)
    net_rates_cache .= 0.0
    #even though react_cell!() already sets change_in_molar_concentrations_cache to .= 0.0, we're just doing .*= 0.0 to all to make sure

    #Connection Loops
    for conn in connection_groups
        solve_connection_group!(
            du, u, conn.phys_a, conn.phys_b, conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            rho_cache, mw_avg_cache,
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
        )
    end
    
    #Controller Loops
    for cont in controller_groups
        solve_controller_group!(
            du, u, cont.id, cont.controller, cont.controller_function!, cont.monitored_cells, cont.affected_cells,
            cell_volumes
        )
    end

    #Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
    for reg in region_groups
        solve_region_group!(
            du, u, reg.phys, reg.region_function!, reg.region_cells,
            rho_cache, mw_avg_cache,
            change_in_molar_concentrations_cache, molar_concentrations_cache, net_rates_cache
        )
    end
end

#VERY IMPORTANT!!!!
#= For future reference when getting properties using u[field]:
    u_cv.temp[cell_id] = val
        works!
    u_cv[field][cell_id] = val
        fails :(
    view(u_cv, field)[cell_id] = val
        works!
    getproperty(u_cv, field)[cell_id] = val
        works!
=#

#instead of updating rho in each flux and internal physics, we could do an initial property retrieval loop
#the only disadvantage is that once again, solids can't do rho_ideal!() and just have phys.rho
#=
@batch for cell_id in eachindex(cell_volumes)
    mw_avg!(u, cell_id, molecular_weights, mw_avg_cache)
    rho_ideal!(u, cell_id, rho_cache, mw_avg_cache)
end
=#