using ComponentArrays
using PreallocationTools
using ForwardDiff
import SparseConnectivityTracer
import ADTypes

u = ComponentVector(test=zeros(10))
du = ComponentVector(test=zeros(10))

rho_cache = zeros(length(10))
mw_avg_cache = zeros(length(10))

n_species = length(zeros(10))
change_in_molar_concentrations_cache = zeros(10)
molar_concentrations_cache = zeros(10) #just using mass fractions for cell 1, this may cause some issues later!

net_rates_cache = zeros(10)

caches = ComponentVector(
    rho_cache=rho_cache,
    mw_avg_cache=mw_avg_cache,
    change_in_molar_concentrations_cache=change_in_molar_concentrations_cache,
    molar_concentrations_cache=molar_concentrations_cache,
    net_rates_cache=net_rates_cache
)

du_caches = copy(caches)

N::Int = ForwardDiff.pickchunksize(length(u))

caches_vec = DiffCache(Vector(caches), N)
du_caches_vec = DiffCache(Vector(du_caches), N)

caches_axes = getaxes(caches)[1]

function ode_test(du, u, p, t, du_caches_vec, caches_vec, caches_axes)
    #I've been trying to get SparseConnectivityTracer to get to work with this for about 3 hours now and I'm so fed up
    #this is really hacky, but I don't give a fuck because I'm so done with trying to get these two to work together
    #THE FUNNY THING IS: I remember these two packages working flawlessly together in the past, but they just don't now and I don't know why
    if eltype(u) == Float64
        caches_vec = get_tmp(caches_vec, Float64)
        du_caches_vec = get_tmp(du_caches_vec, Float64)
    else
        caches_vec = get_tmp(caches_vec, ForwardDiff.Dual(1.0))
        du_caches_vec = get_tmp(du_caches_vec, ForwardDiff.Dual(1.0))
    end

    println(caches_vec)

    caches = ComponentVector(caches_vec, caches_axes)
    du_caches = ComponentVector(du_caches_vec, caches_axes)

    u = ComponentVector(u; NamedTuple(caches)...)
    du = ComponentVector(du; NamedTuple(du_caches)...)
    du .= 0.0

    u.rho_cache .+= 1.0
end

ode_closure(du, u, p, t) = ode_test(du, u, p, t, du_caches_vec, caches_vec, caches_axes)

detector = SparseConnectivityTracer.TracerSparsityDetector()

p_guess = [0.0]

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> ode_closure(du, u, p_guess, 0.0), du, u, detector)
