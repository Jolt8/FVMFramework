using PreallocationTools

@inline create_views_inline(v, ax::NamedTuple) = map(a -> create_views_inline(v, a), ax)
@inline create_views_inline(v, ax::UnitRange) = view(v, ax)

function mass_fraction_update!(du, u, cell_id)
    map(keys(u.mass_fractions)) do species_name
        du.mass_fractions[species_name][cell_id] += 1.0
        u.mass_fractions[species_name][cell_id] += 1.0
    end
end

function cell_update!(du_vec, u_vec, u_diff_cache, state_vec, u_axes, cache_axes, state_axes, cell_ids)
    du_nt = create_views_inline(du_vec, u_axes)
    u_nt = create_views_inline(u_vec, u_axes)
    
    cache_vec = get_tmp(u_diff_cache, u_vec)
    cache_nt = create_views_inline(cache_vec, cache_axes)

    state_nt = create_views_inline(state_vec, state_axes)
    
    du = (; du_nt..., cache_nt...)
    u = (; u_nt..., cache_nt..., state_nt...)
    for cell_id in cell_ids
        du.temp[cell_id] += 1.0
        mass_fraction_update!(du, u, cell_id)
    end
end

u_vec = rand(300);
du_vec = rand(300);
u_axes = (temp = 1:100, mass_fractions = (methanol = 101:200, water = 201:300))

u_cache_vec = rand(200)
u_diff_cache = DiffCache(u_cache_vec)
cache_axes = (molar_concentrations = (methanol = 1:100, water = 101:200))

state_vec = rand(100)
state_axes = (cp = 1:100,)
#don't forget to add the , when there's only one field in the NamedTuple

cell_ids = collect(1:100)

@benchmark create_views_inline(u_vec, u_axes)

@benchmark cell_update!(
    du_vec, u_vec, 
    u_diff_cache, 
    state_vec, 
    u_axes, 
    cache_axes, 
    state_axes, 
    cell_ids
)


function test_f!(du_vec, u_vec, u_diff_cache, state_vec, u_axes, cache_axes, state_axes, cell_ids)
    du_nt = create_views_inline(du_vec, u_axes)
    u_nt = create_views_inline(u_vec, u_axes)
    
    cache_vec = get_tmp(u_diff_cache, u_vec)
    cache_nt = create_views_inline(cache_vec, cache_axes)

    state_nt = create_views_inline(state_vec, state_axes)
    
    du = (; du_nt..., cache_nt...)
    u = (; u_nt..., cache_nt..., state_nt...)
    for cell_id in cell_ids
        du.temp[cell_id] += 1.0
        mass_fraction_update!(du, u, cell_id)
    end
end


function test_f!(du_vec, u_vec, p, t, axes, cell_ids, face_idxs)
    du = create_views_inline(du_vec, axes)
    u = create_views_inline(u_vec, axes)

    @batch for cell_id in cell_ids
        du.temp[cell_id] += 1.0
        du.temp[cell_id] += 1.0 #oh no, just adding this increases the number of allocations from 3388 to 3548.
        #du.mass_fractions.methanol[cell_id] += 1.0 #just adding this (not even in addition to the one before) increases allocations from 3388 to 3842
    end

    @batch for cell_id in cell_ids
        for face_idx in face_idxs[cell_id]
            #update_mass_face(du, cell_id, face_idx)
        end
    end
end

test_f_closure! = (du, u, p, t) -> test_f!(du, u, p, t, u_axes, cell_ids, face_idxs)

prob = ODEProblem(test_f_closure!, u_vec, (0.0, 10000.0), [0.0])
@time sol = solve(prob, Tsit5())
