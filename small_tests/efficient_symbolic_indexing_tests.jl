using ComponentArrays
using Polyester
using BenchmarkTools

n_cells = 10000

du_proto = ComponentArray(
    mass_fractions=ComponentVector(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp=zeros(n_cells)
)

du_proto_nt = (; du_proto...)

u_proto = ComponentArray(
    mass_fractions=ComponentVector(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp=zeros(n_cells)
)

u_proto_nt = (; u_proto...)

du_cache = ComponentArray(
    rho=zeros(n_cells),)

du_cache_nt = (; du_cache...)

u_cache = ComponentArray(
    rho=zeros(n_cells)
)

u_cache_nt = (; u_cache...)

properties = (cp=zeros(n_cells), k=zeros(n_cells))

du = ComponentVector(du_proto; NamedTuple(du_cache)..., properties...)
u = ComponentVector(u_proto; NamedTuple(u_cache)..., properties...)

du_axes = getaxes(du)[1]
u_axes = getaxes(u)[1]

du = Vector(du)
u = Vector(u)

#ahhh, so turning them all into a component vector is what's causing multithreading to fail

function test_accessors(du_proto_nt, du_cache_nt, u_proto_nt, u_cache_nt, properties, n_cells)
    du = (; du_proto_nt..., du_cache_nt...)
    u = (; u_proto_nt..., u_cache_nt..., properties...)
    for cell_id in 1:n_cells
        du.mass_fractions
        du.mass_fractions.methanol
        du.mass_fractions.methanol[cell_id]
        du.mass_fractions.methanol[cell_id] = 1.0
        u.mass_fractions
        u.mass_fractions.methanol
        u.mass_fractions.methanol[cell_id]
        u.mass_fractions.methanol[cell_id] = 1.0

        du.temp
        du.temp[cell_id]
        du.temp[cell_id] = 270.0
        u.temp
        u.temp[cell_id]
        u.temp[cell_id] = 270.0

        du.rho
        du.rho[cell_id]
        du.rho[cell_id] = 1.0
        u.rho
        u.rho[cell_id]
        u.rho[cell_id] = 1.0

        u.cp
        u.cp[cell_id]
        u.cp[cell_id] = 1.0
    end
end

function test_multiple_samples(du_proto, du_cache, u_proto, u_cache, properties, n_cells, n_samples)
    for _ in 1:n_samples
        test_accessors(du_proto, du_cache, u_proto, u_cache, properties, n_cells)
    end
end

#VSCodeServer.@profview test_multiple_samples(du_proto_nt, du_cache_nt, u_proto_nt, u_cache_nt, properties, n_cells, 1000)
#@btime test_accessors($du_proto_nt, $du_cache_nt, $u_proto_nt, $u_cache_nt, $properties, $n_cells)

function test_ode(du_vec, u_vec, p, t, du_cache_nt, u_cache_nt, properties, du_proto_axes, u_proto_axes, n_cells)
    du_nt = (; NamedTuple(ComponentArray(du_vec, du_proto_axes))...)
    u_nt = (; NamedTuple(ComponentArray(u_vec, u_proto_axes))...)

    du = (; du_nt..., du_cache_nt...)
    u = (; u_nt..., u_cache_nt..., properties...)
    @batch for cell_id in 1:n_cells
        du.mass_fractions.methanol[cell_id] += 1.0
        du.temp[cell_id] = 270.0
        du.rho[cell_id] = 1.0
        u.mass_fractions.methanol[cell_id] = 1.0
        u.temp[cell_id] = 270.0
        u.rho[cell_id] = 1.0
        u.cp[cell_id] = 1.0
    end
end

du_vec = Vector(ComponentArray(; du_proto_nt...))
u_vec = Vector(ComponentArray(; u_proto_nt...))

du_proto_axes = getaxes(ComponentArray(; du_proto_nt...))[1]
u_proto_axes = getaxes(ComponentArray(; u_proto_nt...))[1]

test_ode_closure = (du, u, p, t) -> test_ode(
    du, u, p, t,
    du_cache_nt, u_cache_nt, properties,
    du_proto_axes, u_proto_axes,
    1000
)

p_guess = [0.0]
t = 0.0

#VSCodeServer.@profview test_ode_closure(du_vec, u_vec, p_guess, t)
@btime test_ode_closure($du_vec, $u_vec, $p_guess, $t)

using OrdinaryDiffEq

test_prob = ODEProblem(test_ode_closure, Vector(ComponentArray(; u_proto_nt...)), (0.0, 0.00001), p_guess)
@time sol = solve(test_prob, Tsit5(), tspan = (0.0, 0.00001))
@btime sol = solve(test_prob, Tsit5(), tspan = (0.0, 0.00001))
