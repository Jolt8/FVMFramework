using ComponentArrays
using BenchmarkTools
using PreallocationTools
using Polyester

n_cells = 1000

du_proto = ComponentArray(
    mass_fractions = ComponentArray(methanol = zeros(n_cells), water = zeros(n_cells), carbon_monoxide = zeros(n_cells), hydrogen = zeros(n_cells), carbon_dioxide = zeros(n_cells)),
    temp = zeros(n_cells)
)

u_proto = ComponentArray(
    mass_fractions = ComponentArray(methanol = zeros(n_cells), water = zeros(n_cells), carbon_monoxide = zeros(n_cells), hydrogen = zeros(n_cells), carbon_dioxide = zeros(n_cells)),
    temp = zeros(n_cells)
)

du_cache_nt = ComponentArray(
    rho = zeros(n_cells),
)

u_cache_nt = ComponentArray(
    rho = zeros(n_cells),
)

properties = ComponentArray(rho = ones(n_cells) * 1100.0, k = 0.15)

du_vec = Vector(du_proto)
u_vec = Vector(u_proto)
du_cache_vec = Vector(du_cache_nt)
u_cache_vec = Vector(u_cache_nt)

du_proto_axes = getaxes(du_proto)
u_proto_axes = getaxes(u_proto)
du_cache_axes = getaxes(du_cache_nt)
u_cache_axes = getaxes(u_cache_nt)

du_diff_cache_vec = DiffCache(du_cache_vec)
u_diff_cache_vec = DiffCache(u_cache_vec)

function test_ode(du_vec, u_vec, p, t, du_diff_cache_vec, u_diff_cache_vec, properties, du_proto_axes, u_proto_axes, du_cache_axes, u_cache_axes, n_cells)
    du_cache_vec = get_tmp(du_diff_cache_vec, u_vec)
    u_cache_vec = get_tmp(u_diff_cache_vec, u_vec)

    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = ComponentArray(du_cache_vec, du_cache_axes)
    u_cache_nt = ComponentArray(u_cache_vec, u_cache_axes)

    du_vec .= 0.0
    du_nt = ComponentArray(du_vec, du_proto_axes)
    u_nt = ComponentArray(u_vec, u_proto_axes)

    du = ComponentArray(; du_nt..., du_cache_nt...)
    u = ComponentArray(; properties..., u_nt..., u_cache_nt...)

    u.rho .= (u_cache_nt.rho .+ properties.rho)

    @batch for cell_id in 1:n_cells
        du.mass_fractions.methanol[cell_id] += 1.0
        du.temp[cell_id] += 270.0
        du.rho[cell_id] += 1.0
        u.mass_fractions.methanol[cell_id] += 1.0
        u.temp[cell_id] += 270.0
        u.rho[cell_id] += 1.0
        u.k += 1.0
        map(keys(du.mass_fractions)) do species_name
            #ok, if we run this without anything inside it, we still get zero allocations
            #ahh, so using keys returned by map(keys) when doing @batch causes a lot of allocations
            #otherwise, it works fine
            #that really sucks
        end
        map(keys(du.mass_fractions)) do species_name
            du.mass_fractions[species_name][cell_id] += 1.0
        end
        #map(keys(u.mass_fractions)) do species_name #fuck, not using this map makes this way faster
            du.mass_fractions[:methanol][cell_id] += 1.0
            du.mass_fractions[:water][cell_id] += 1.0
            du.mass_fractions[:carbon_monoxide][cell_id] += 1.0
            du.mass_fractions[:hydrogen][cell_id] += 1.0
            du.mass_fractions[:carbon_dioxide][cell_id] += 1.0
        #end
    end
end

test_ode_closure = (du, u, p, t) -> test_ode(
    du, u, p, t,
    du_diff_cache_vec, u_diff_cache_vec, properties,

    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes,
    1000
)

p_guess = [0.0]
t = 0.0

#VSCodeServer.@profview test_ode_closure(du_vec, u_vec, p_guess, t)
test_ode_closure(du_vec, u_vec, p_guess, t)
@btime test_ode_closure($du_vec, $u_vec, $p_guess, $t) #18.093 ms (116944 allocations: 80.17 MiB)
VSCodeServer.@profview [test_ode_closure(du_vec, u_vec, p_guess, t) for _ in 1:100]
