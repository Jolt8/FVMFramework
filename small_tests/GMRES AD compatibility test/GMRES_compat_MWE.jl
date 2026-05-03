using OrdinaryDiffEq
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ForwardDiff
using Plots
using SparseConnectivityTracer
import ADTypes
using NonlinearSolve
using BenchmarkTools

function test_ode_function(du, u, p, t)
    du[1] += p[1] * u[1] + p[2] * u[2]
    du[2] += p[3] * u[1] + p[4] * u[2]

    du[1] -= u[1] * p[1]
    du[2] -= u[2] * p[3]

    return
end

du0 = [0.0, 0.0]
u0 = [1.0, 1.0]
tspan = (0.0, 1.0)
p = [1.0, 1.0, 1.0, 1.0]

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> test_ode_function(du, u, p, 0.0), du0, u0, detector
)

ode_func = ODEFunction(test_ode_function, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 1.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0, tspan, p)

@btime sol_krylov = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES())) #144.000 μs (501 allocations: 25.33 KiB)
@btime sol_klu = solve(implicit_prob, FBDF(linsolve = KLUFactorization())) #761.300 μs (295 allocations: 20.31 KiB)

function auto_forward_diff_loss(p)
    prob = ODEProblem(ode_func, u0, tspan, p)
    sol = solve(
        prob, 
        FBDF(linsolve = KLUFactorization()),
        sensealg = InterpolatingAdjoint(autodiff = AutoForwardDiff())
    )
    return sum(sol.u[1])
end

@time ForwardDiff.gradient(auto_forward_diff_loss, p) #1.711 ms (366 allocations: 28.83 KiB)

function auto_enzyme_loss(p)
    prob = ODEProblem(ode_func, u0, tspan, p)
    sol = solve(
        prob, 
        FBDF(linsolve = KrylovJL_GMRES()),
        sensealg = InterpolatingAdjoint(autodiff = AutoEnzyme())
    )
    return sum(sol.u[1])
end

@time ForwardDiff.gradient(auto_enzyme_loss, p) #1.731 ms (366 allocations: 28.83 KiB)

krylov_solve!()