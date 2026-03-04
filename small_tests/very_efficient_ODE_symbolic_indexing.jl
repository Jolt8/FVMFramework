using ComponentArrays
using Polyester
using BenchmarkTools
using OrdinaryDiffEq
using PreallocationTools

#NOTE: I've found out that accessing a ComponentArray within an ODE function is what's causing multithreading to fail

n_cells = 10000

du_proto = ComponentArray(
    mass_fractions = ComponentVector(methanol = zeros(n_cells), water = zeros(n_cells), carbon_monoxide = zeros(n_cells), hydrogen = zeros(n_cells), carbon_dioxide = zeros(n_cells)),
    temp = zeros(n_cells)
)

du_proto_nt = (; du_proto...)

u_proto = ComponentArray(
    mass_fractions=ComponentVector(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp = zeros(n_cells)
)

u_proto_nt = (; u_proto...)

du_cache = ComponentArray(
    rho = zeros(n_cells),
)

du_cache_nt = (; du_cache...)

u_cache = ComponentArray(
    rho = zeros(n_cells)
)

u_cache_nt = (; u_cache...)

properties = (rho = ones(n_cells) * 1100.0, cp = ones(n_cells) * 2500.0, k = ones(n_cells) * 0.15)

du_vec = Vector(ComponentArray(; du_proto_nt...))
u_vec = Vector(ComponentArray(; u_proto_nt...))
du_cache_vec = Vector(ComponentArray(; du_cache_nt...))
u_cache_vec = Vector(ComponentArray(; u_cache_nt...))

du_proto_axes = getaxes(ComponentArray(; du_proto_nt...))[1]
u_proto_axes = getaxes(ComponentArray(; u_proto_nt...))[1]
du_cache_axes = getaxes(ComponentArray(; du_cache_nt...))[1]
u_cache_axes = getaxes(ComponentArray(; u_cache_nt...))[1]

du_diff_cache_vec = DiffCache(du_cache_vec)
u_diff_cache_vec = DiffCache(u_cache_vec)

using ForwardDiff

function test_ode(du_vec, u_vec, p, t, du_diff_cache_vec, u_diff_cache_vec, properties, du_proto_axes, u_proto_axes, du_cache_axes, u_cache_axes, n_cells)
    du_cache_vec = get_tmp(du_diff_cache_vec, ForwardDiff.Dual(1.0))
    u_cache_vec = get_tmp(u_diff_cache_vec, ForwardDiff.Dual(1.0))

    du_cache_vec .= 0.0
    u_cache_vec .= 0.0

    du_cache_nt = (; NamedTuple(ComponentArray(du_cache_vec, du_cache_axes))...)
    u_cache_nt = (; NamedTuple(ComponentArray(u_cache_vec, u_cache_axes))...)

    du_nt = (; NamedTuple(ComponentArray(du_vec, du_proto_axes))...)
    u_nt = (; NamedTuple(ComponentArray(u_vec, u_proto_axes))...)

    du = (; du_nt..., du_cache_nt...)
    u = (; u_nt..., u_cache_nt..., properties...)

    @batch for cell_id in 1:n_cells
        du.mass_fractions.methanol[cell_id] += 1.0
        du.temp[cell_id] += 270.0
        du.rho[cell_id] += 1.0
        u.mass_fractions.methanol[cell_id] += 1.0
        u.temp[cell_id] += 270.0
        u.rho[cell_id] += 1.0
        u.cp[cell_id] += 1.0
        for (species_name, species_mass_fraction) in pairs(u.mass_fractions) #this is bad, do not do this
            du.mass_fractions[species_name][cell_id] += 1.0
        end
    end
end

test_ode_closure = (du, u, p, t) -> test_ode(
    du, u, p, t,
    du_diff_cache_vec, u_diff_cache_vec, properties,

    du_proto_axes, u_proto_axes,
    du_cache_axes, u_cache_axes,
    10000
)

p_guess = [0.0]
t = 0.0

#VSCodeServer.@profview test_ode_closure(du_vec, u_vec, p_guess, t)
test_ode_closure(du_vec, u_vec, p_guess, t)
@btime test_ode_closure($du_vec, $u_vec, $p_guess, $t)

test_prob = ODEProblem(test_ode_closure, Vector(ComponentArray(; u_proto_nt...)), (0.0, 0.00001), p_guess)
@time sol = solve(test_prob, Tsit5(), tspan = (0.0, 0.00001))
@btime sol = solve(test_prob, Tsit5(), tspan = (0.0, 0.00001))
VSCodeServer.@profview sol = solve(test_prob, Tsit5(), tspan = (0.0, 0.00001))
