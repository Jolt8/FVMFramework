using Symbolics 
using MethodOfLines
using Unitful
using Plots
using ModelingToolkit
using DifferentialEquations
using DomainSets

R_gas = ustrip(upreferred(8.314u"J/(K*mol)"))

mw_CH3OH = ustrip(upreferred(32.04u"g/mol"))
mw_H2O = ustrip(upreferred(18.02u"g/mol"))
mw_CO = ustrip(upreferred(28.01u"g/mol"))
mw_H2 = ustrip(upreferred(2.02u"g/mol"))
mw_CO2 = ustrip(upreferred(44.01u"g/mol"))
mw_air = ustrip(upreferred(28.97u"g/mol"))

cat_bulk_density = ustrip(upreferred(1250.0u"kg/m^3"))

pipe_inside_diameter = ustrip(upreferred(0.5u"inch"))
pipe_length = ustrip(upreferred(12.1u"inch"))
cross_sectional_area = pi * pipe_inside_diameter^2 / 4

bed_void_fraction = 0.4

cp = ustrip(upreferred(50.0u"J/(mol*K)"))

#reactions 
kf_A_MSR = ustrip(upreferred(1.59e10u"s^-1"))
kf_Ea_MSR = ustrip(upreferred(103000.0u"J/mol"))
k_eq_MSR_ref = ustrip(upreferred(1.5e9u"s^-1"))
ref_T_MSR = ustrip(upreferred(298.15u"K"))
Heat_rxn_MSR = ustrip(upreferred(49500.0u"J/mol"))

kf_A_WGS = ustrip(upreferred(4.63e10u"s^-1"))
kf_Ea_WGS = ustrip(upreferred(87500.0u"J/mol"))
k_eq_WGS_ref = ustrip(upreferred(4.6e9u"s^-1"))
ref_T_WGS = ustrip(upreferred(298.15u"K"))
Heat_rxn_WGS = ustrip(upreferred(-41100.0u"J/mol"))

kf_A_MD = ustrip(upreferred(1.46e13u"s^-1"))
kf_Ea_MD = ustrip(upreferred(170000.0u"J/mol"))
k_eq_MD_ref = ustrip(upreferred(1.46e12u"s^-1"))
ref_T_MD = ustrip(upreferred(298.15u"K"))
Heat_rxn_MD = ustrip(upreferred(90200.0u"J/mol"))


@independent_variables t x
@variables molar_flow_CH3OH(..) molar_flow_H2O(..) molar_flow_CO(..) molar_flow_H2(..) molar_flow_CO2(..) molar_flow_air(..) T(..) P(..)

Dt = Differential(t)
Dx = Differential(x)

molar_flow_total = molar_flow_CH3OH(t, x) + molar_flow_H2O(t, x) + molar_flow_CO(t, x) + molar_flow_H2(t, x) + molar_flow_CO2(t, x) + molar_flow_air(t, x)

y_CH3OH = molar_flow_CH3OH(t, x) / molar_flow_total
y_H2O = molar_flow_H2O(t, x) / molar_flow_total
y_CO = molar_flow_CO(t, x) / molar_flow_total
y_H2 = molar_flow_H2(t, x) / molar_flow_total
y_CO2 = molar_flow_CO2(t, x) / molar_flow_total
y_air = molar_flow_air(t, x) / molar_flow_total

p_CH3OH = y_CH3OH * P(t, x)
p_H2O = y_H2O * P(t, x)
p_CO = y_CO * P(t, x)
p_H2 = y_H2 * P(t, x)
p_CO2 = y_CO2 * P(t, x)
p_air = y_air * P(t, x)

mw_avg = y_CH3OH * mw_CH3OH + y_H2O * mw_H2O + y_CO * mw_CO + y_H2 * mw_H2 + y_CO2 * mw_CO2 + y_air * mw_air

density = (P(t, x) * mw_avg) / (R_gas * T(t, x))

kf_MSR = kf_A_MSR * exp(-kf_Ea_MSR/(R_gas*T(t, x)))
k_eq_MSR = k_eq_MSR_ref * exp(-(Heat_rxn_MSR)/R_gas * (1/T(t, x) - 1/ref_T_MSR))
r_net_MSR = kf_MSR * (p_CH3OH * p_H2O / (p_CO * p_H2^3) - k_eq_MSR)

kf_WGS = kf_A_WGS * exp(-kf_Ea_WGS/(R_gas*T(t, x)))
k_eq_WGS = k_eq_WGS_ref * exp(-(Heat_rxn_WGS)/R_gas * (1/T(t, x) - 1/ref_T_WGS))
r_net_WGS = kf_WGS * (p_CO * p_H2O / (p_CO2 * p_H2) - k_eq_WGS)

kf_MD = kf_A_MD * exp(-kf_Ea_MD/(R_gas*T(t, x)))
k_eq_MD = k_eq_MD_ref * exp(-(Heat_rxn_MD)/R_gas * (1/T(t, x) - 1/ref_T_MD))
r_net_MD = kf_MD * (p_CH3OH / (p_CO * p_H2^2) - k_eq_MD)

interstitial_velocity = (R_gas * T(t, x) * molar_flow_total) / (P(t, x) * (1-bed_void_fraction) * cross_sectional_area)

source_CH3OH = cat_bulk_density * cross_sectional_area * (-1 * r_net_MSR - 1 * r_net_MD)
source_H2O = cat_bulk_density * cross_sectional_area * (-1 * r_net_MSR - 1 * r_net_WGS)
source_CO = cat_bulk_density * cross_sectional_area * (r_net_MSR - r_net_WGS + r_net_MD)
source_H2 = cat_bulk_density * cross_sectional_area * (3 * r_net_MSR + r_net_WGS + 2 * r_net_MD)
source_CO2 = cat_bulk_density * cross_sectional_area * (r_net_MSR + r_net_WGS)
source_air = cat_bulk_density * cross_sectional_area * (0)

interstitial_velocity = (R_gas * T(t, x) * molar_flow_total) / (P(t, x) * (1-bed_void_fraction) * cross_sectional_area)

eq_F_CH3OH = Dt(molar_flow_CH3OH(t, x)) ~ interstitial_velocity * -Dx(molar_flow_CH3OH(t, x)) + source_CH3OH
eq_F_H2O = Dt(molar_flow_H2O(t, x)) ~ interstitial_velocity * -Dx(molar_flow_H2O(t, x)) + source_H2O
eq_F_CO = Dt(molar_flow_CO(t, x)) ~ interstitial_velocity * -Dx(molar_flow_CO(t, x)) + source_CO
eq_F_H2 = Dt(molar_flow_H2(t, x)) ~ interstitial_velocity * -Dx(molar_flow_H2(t, x)) + source_H2
eq_F_CO2 = Dt(molar_flow_CO2(t, x)) ~ interstitial_velocity * -Dx(molar_flow_CO2(t, x)) + source_CO2
eq_F_air = Dt(molar_flow_air(t, x)) ~ interstitial_velocity * -Dx(molar_flow_air(t, x)) + source_air
#eq_DT_Dt = Dt(T(t, x)) ~ (Dx(molar_flow_total * T(t, x)) * cp) / (cross_sectional_area * cat_bulk_density)
eq_DP_Dx = Dx(P(t, x)) ~ 0.0
eq_DT_Dt = Dt(T(t, x)) ~ 0.0

eqs = [eq_F_CH3OH, eq_F_H2O, eq_F_CO, eq_F_H2, eq_F_CO2, eq_F_air, eq_DT_Dt, eq_DP_Dx]

t_start = 0.0
t_end = 1000.0

x_min = 0.0
x_max = pipe_length

domains = [
    t ∈ Interval(t_start, t_end),
    x ∈ Interval(x_min, x_max)
]

bcs = [
    molar_flow_CH3OH(t_start, x) ~ 0.001,
    molar_flow_H2O(t_start, x) ~ 0.001,
    molar_flow_CO(t_start, x) ~ 0.001,
    molar_flow_H2(t_start, x) ~ 0.001,
    molar_flow_CO2(t_start, x) ~ 0.001,
    molar_flow_air(t_start, x) ~ 0.001,
    T(t_start, x) ~ 1.0,
    P(t_start, x) ~ 573.15,

    molar_flow_CH3OH(t, x_min) ~ 1.0,
    molar_flow_H2O(t, x_min) ~ 1.3,
    molar_flow_CO(t, x_min) ~ 0.001,
    molar_flow_H2(t, x_min) ~ 0.001,
    molar_flow_CO2(t, x_min) ~ 0.001,
    molar_flow_air(t, x_min) ~ 0.001,
    T(t, x_min) ~ 1.0,
    P(t, x_min) ~ 573.15,
]

@named pdesystem = PDESystem(
    eqs, 
    bcs, 
    domains, 
    [t, x], 
    [
        molar_flow_CH3OH(t, x), 
        molar_flow_H2O(t, x), 
        molar_flow_CO(t, x), 
        molar_flow_H2(t, x), 
        molar_flow_CO2(t, x), 
        molar_flow_air(t, x), 
        T(t, x), 
        P(t, x)
    ],
)

x_points = 10
discretization = MOLFiniteDifference([x => x_points], t) 
@time prob = discretize(pdesystem, discretization)
#holy moly this takes absolutely forever 

#prob = structural_simplify(prob)

sol = solve(prob, Tsit5())

sol = solve(prob, Rosenbrock23())
#sol = solve(prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero(), ))