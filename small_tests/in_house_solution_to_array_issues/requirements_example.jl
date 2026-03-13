using ComponentArrays 
using PreallocationTools
using Polyester

function ode_for_testing_f!(
    du_vec, u_vec, p, t, 
    p_axes,

    properties,
    
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    du_diff_cache_vec, u_diff_cache_vec,
    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes
)
    #ABSOLUTELY REQUIRED


    #START OF SETUP
    #get_tmp cached variables (ABSOLUTELY REQUIRED)
    du_cache_vec = get_tmp(du_diff_cache_vec, first(u_vec) + first(p))
    u_cache_vec = get_tmp(u_diff_cache_vec, first(u_vec) + first(p))
    
    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    #NOT ABSOLUTELY REQUIRED
    #Note that the specific syntax of create_views_inline is not required and would be better if it was not required. 
    #Whether or not a flat vector is passed in with axes or a data structure that already behaves like a vector and a set of axes is up to you

    #NOT ABSOLUTELY REQUIRED
    du_cache_in = create_views_inline(du_cache_vec, du_cache_axes)
    u_cache_in = create_views_inline(u_cache_vec, u_cache_axes)

    #NOT ABSOLUTELY REQUIRED
    du_vec .= 0.0
    du_in = create_views_inline(du_vec, du_proto_axes)
    u_in = create_views_inline(u_vec, u_proto_axes)

    #merging properties into something that can be referenced by all physics functions (doesn't have to be a named tuple)
    #ABSOLUTELY REQUIRED
    du = (; du_in..., du_cache_in...)
    u = (; properties..., u_in..., u_cache_in...) 

    #ideally the step below would not be required and properties would automatically the values (BUT NOT THE DUALS) of du_cache/u_cache
    #ABSOLUTELY REQUIRED
    u.rho .= (u_cache_in.rho .+ properties.rho)

    p_named = create_views_inline(p, p_axes)

    u.diffusion_pre_exponential_factor .= p_named.diffusion_pre_exponential_factor[1]
    u.diffusion_activation_energy .= p_named.diffusion_activation_energy[1]

    #it would be ideal if these were merged without allocations
    
    #END OF SETUP

    #CAPABILITIES THAT REQUIRE 0 ALLOCATIONS
    @batch for cell_id in 1:length(cell_volumes) #supports @batch
        du.mass_fractions[:methylene_blue][cell_id] += 1.0 #symbolic indexing
        du.mass_fractions.methylene_blu[cell_id] += 1.0 #dot indexing

        du[1] #looks or is an array under the hood to @batch can use it

        u.net_rates.reforming_reactions.WGS_rxn += 1.0 #supports non cell_id indexed fields

        for species_name in keys(u.mass_fractions)
            du.mass_fractions[species_name][cell_id] += 1.0 #supports fast looping through fields without allocations
            for reaction in keys(u.net_rates.reforming_reactions) #supports nested field looping
                du.net_rates.reforming_reactions[reaction] += 1.0 #supports fast looping through single value fields without allocations
            end
        end

        for face_idx in eachindex(cell_neighbor_areas)
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


n_cells = 100


ode_for_testing_f!(
    du_vec, u_vec, p, t,

    properties,
    
    cell_volumes, cell_centroids,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    unconnected_cell_face_map, cell_face_areas, cell_face_normals,

    p_axes,
    du_diff_cache_vec, u_diff_cache_vec,
    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes
)

n_cells = 1000

du0_vec = ComponentArray(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

u0_vec = ComponentArray(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells)
    ),
    pressure = zeros(n_cells)
)

du_caches = ComponentArray(
        u_caches = (
        mass_face = fill(zeros(n_faces), n_cells),
        net_rates = (
            reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
        ), #it would also be nice if commans were not required for fields with 1 field within them
    )
)

u_caches = ComponentArray(
        u_caches = (
        mass_face = fill(zeros(n_faces), n_cells),
        net_rates = (
            reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))),
        ), #it would also be nice if commans were not required for fields with 1 field within them
    )
)

du_vec = Vector(du_proto)
u_vec = Vector(u_proto)
du_cache_vec = Vector(du_caches)
u_cache_vec = Vector(u_caches)

du_proto_axes = getaxes(du_proto)
u_proto_axes = getaxes(u_proto)
du_cache_axes = getaxes(du_caches)
u_cache_axes = getaxes(u_caches)

du_diff_cache_vec = DiffCache(du_cache_vec)
u_diff_cache_vec = DiffCache(u_cache_vec)