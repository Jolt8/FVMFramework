using PreallocationTools
using ComponentArrays
using ForwardDiff

N = 10

caches = DiffCache(
    ComponentVector(
        rho_cache = zeros(N),
        mw_avg_cache = zeros(N),
        change_in_molar_concentrations_cache = zeros(N),
        molar_concentrations_cache = zeros(N, 5),
        net_rates_cache = zeros(N, 3)
    ), 15
)

u_proto = ComponentArray(
    vel_x=zeros(N), vel_y=zeros(N), vel_z=zeros(N),
    pressure=zeros(N),
    mass_fractions=zeros(5, N),
    temp=zeros(N)
)

u_test = ForwardDiff.Dual(1)

test = get_tmp(caches, u_proto)

test.rho_cache

@time u = ComponentVector(u_proto; NamedTuple(test)...)

for i in eachindex(u.mass_fractions)
    println("hi")
end

@time u = [u_proto; test]

u.rho_cache

zoop = copy(u)

zoop.vel_x[1] == u

u.rho_cache

#I was really suprised that this works