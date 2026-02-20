#Structs
#=
struct VarPath
    root::Symbol #:du or :u
    path::Vector{Symbol} #[:reactions, :MSR_rxn, :kf_A]
end

struct VarAccess
    path::VarPath
    access_type::Symbol #:read or :write
end

struct TracerContext
    access_logs::Vector{VarAccess}
    pending_reads::Vector{VarPath}
    write_dependencies::Dict{VarPath, Set{VarPath}}
    loop_context::Dict{Symbol, Vector{Symbol}}
end

struct TracerVariable
    path::VarPath
    tracer_ref::Ref{TracerContext}
    value::Float64
end

#Functions and Overloads

function Base.getproperty(t::TracerVariable, sym::Symbol)
    parent_path = getfield(t, :path)
    new_path = VarPath(parent_path.root, vcat(parent_path.path, sym))

    return TracerVariable(new_path, getfield(t, :tracer_ref)[], rand())
end

function Base.setproperty!(t::TracerVariable, sym::Symbol, val)
    ctx = getfield(t, :tracer_ref)[]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, sym))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)

    empty!(ctx.pending_reads)
end

function Base.getindex(t::TracerVariable, idx::Int)
    parent_path = getfield(t, :path)
    return TracerVariable(
        VarPath(parent_path.root, vcat(parent_path.path, :_idx)),
        getfield(t, :tracer_ref),
        getfield(t, :value)
    )
end

function Base.getindex(t::TracerVariable, sym::Symbol)
    parent_path = getfield(t, :path)
    return TracerVariable(
        VarPath(parent_path.root, vcat(parent_path.path, sym)),
        getfield(t, :tracer_ref),
        getfield(t, :value)
    )
end

function Base.setindex!(t::TracerVariable, val, i::Int)
    ctx = getfield(t, :tracer_ref)[]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, :_idx))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)
    empty!(ctx.pending_reads)
end

function Base.setindex!(t::TracerVariable, val, sym::Symbol)
    ctx = getfield(t, :tracer_ref)[]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, sym))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)
    empty!(ctx.pending_reads)
end

function logread(t::TracerVariable)
    ctx = getfield(t, :tracer_ref)[]
    vp = getfield(t, :path)
    push!(ctx.access_logs, VarAccess(vp, :read))
    push!(ctx.pending_reads, vp)
end

function propertynames(t::TracerVariable) #this is just to make sure it doesn't yield path, tracer_ref, value
    ctx = getfield(t, :tracer_ref)[]
    context = ctx.loop_context[getfield(t, :path).path[end]]
    return context
end

#MATH
Base.:+(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) + b
)
Base.:-(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) - b
)
Base.:*(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) * b
)
Base.:/(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) / b
)
Base.:^(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) ^ b
)

Base.:+(a::Number, b::TracerVariable) = (logread(b); return a + getfield(b, :value))
Base.:-(a::Number, b::TracerVariable) = (logread(b); return a - getfield(b, :value))
Base.:*(a::Number, b::TracerVariable) = (logread(b); return a * getfield(b, :value))
Base.:/(a::Number, b::TracerVariable) = (logread(b); return a / getfield(b, :value))
Base.:^(a::Number, b::TracerVariable) = (logread(b); return a ^ getfield(b, :value))

Base.:+(a::TracerVariable, b::TracerVariable) = (logread(a); logread(b); 
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) + getfield(b, :value)
)
)
Base.:-(a::TracerVariable, b::TracerVariable) = (logread(a); logread(b); 
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) - getfield(b, :value)
)
)
Base.:*(a::TracerVariable, b::TracerVariable) = (logread(a); logread(b); 
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) * getfield(b, :value)
)
)
Base.:/(a::TracerVariable, b::TracerVariable) = (logread(a); logread(b); 
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) / getfield(b, :value)
)
)
Base.:^(a::TracerVariable, b::TracerVariable) = (logread(a); logread(b); 
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
        getfield(a, :value) ^ getfield(b, :value)
)
)
Base.:sqrt(a::TracerVariable) = (logread(a); return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    sqrt(getfield(a, :value))
)
)
Base.:exp(a::TracerVariable) = (logread(a); return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    exp(getfield(a, :value))
)
)
Base.:-(a::TracerVariable) = (logread(a); return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    -getfield(a, :value)
)
)


ctx_loop_context = Dict(
    :mass_fractions => [:methanol, :water, :carbon_monoxide, :hydrogen, :carbon_dioxide],
    :molar_concentrations => [:methanol, :water, :carbon_monoxide, :hydrogen, :carbon_dioxide],
    :reforming_reactions => [:MSR_rxn, :MD_rxn, :WGS_rxn]
)

ctx = TracerContext(VarAccess[], VarPath[], Dict(), ctx_loop_context)

u_tracer = TracerVariable(
    VarPath(:u, Symbol[]),
    ctx,
    1.0
)

du_tracer = TracerVariable(
    VarPath(:du, Symbol[]),
    ctx,
    1.0
)

function molar_concentrations!(u, cell_id)
    for species_name in propertynames(u.molar_concentrations)
        u.molar_concentrations[species_name][cell_id] = u.mass_fractions[species_name][cell_id] / u.molecular_weights[species_name]
    end
end

function mw_avg!(u, cell)
    for species_name in propertynames(u.mass_fractions)
        u.mw_avg[cell] += u.mass_fractions[species_name][cell] / u.molecular_weights[species_name]
    end

    u.mw_avg[cell] = u.mw_avg[cell]^-1.0
end

function rho_ideal!(u, cell)
    u.rho[cell] = (u.pressure[cell] * u.mw_avg[cell]) / (R_gas * u.temp[cell])
end

function van_t_hoff(A, dH, T)
    #K = A * exp(dH/RT)
    return A * exp(-dH / (R_gas * T))
end

function K_gibbs_free(u, cell, reaction)
    K_ref = exp(-reaction.ΔG_rxn_ref / (8.314e-3 * reaction.T_ref)) #R is in kJ

    ln_K_ratio = (-reaction.ΔH_rxn_ref / 8.314e-3) * (1 / u.temp[cell] - 1 / reaction.T_ref)

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

function PAM_reforming_react_cell!(du, u, cell_id, vol)
    #we divide by 1e-5 because PAM parameters are typically in bar
    # methanol, water, carbon_monoxide, hydrogen, carbon_dioxide
    conversion_factor = ((R_gas * u.temp[cell_id]) * 1e-5)
    P_CH3OH = u.molar_concentrations.methanol[cell_id] * conversion_factor
    P_H2O = u.molar_concentrations.water[cell_id] * conversion_factor
    P_CO = u.molar_concentrations.carbon_monoxide[cell_id] * conversion_factor
    P_H2 = u.molar_concentrations.hydrogen[cell_id] * conversion_factor
    P_CO2 = u.molar_concentrations.carbon_dioxide[cell_id] * conversion_factor

    K_CH3O = van_t_hoff(u.adsorption_A_vec.CH3O, u.adsorption_dH_vec.CH3O, u.temp[cell_id])
    K_HCOO = van_t_hoff(u.adsorption_A_vec.HCOO, u.adsorption_dH_vec.HCOO, u.temp[cell_id])
    K_OH = van_t_hoff(u.adsorption_A_vec.OH, u.adsorption_dH_vec.OH, u.temp[cell_id])

    term_CH3O = K_CH3O * (P_CH3OH / sqrt(P_H2))
    term_HCOO = K_HCOO * (P_CO2 * sqrt(P_H2))
    term_OH = K_OH * (P_H2O / sqrt(P_H2))

    DEN = 1.0 + term_CH3O + term_HCOO + term_OH

    #MSR
    k_MSR = u.reactions.reforming_reactions.MSR_rxn.kf_A * exp(-u.reactions.reforming_reactions.MSR_rxn.kf_Ea / (R_gas * u.temp[cell_id]))
    K_eq_MSR = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.MSR_rxn)

    driving_force = 1.0 - ((P_CO2 * P_H2^3) / (K_eq_MSR * P_CH3OH * P_H2O))

    u.net_rates.reforming_reactions.MSR_rxn = (k_MSR * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2) #max just to ensure we don't divide by zero

    #MD
    k_MD = u.reactions.reforming_reactions.MD_rxn.kf_A * exp(-u.reactions.reforming_reactions.MD_rxn.kf_Ea / (R_gas * u.temp[cell_id]))
    K_eq_MD = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.MD_rxn)
    driving_force = 1.0 - ((P_CO * P_H2^2) / (K_eq_MD * P_CH3OH))

    u.net_rates.reforming_reactions.MD_rxn = (k_MD * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2)

    #WGS
    k_WGS = u.reactions.reforming_reactions.WGS_rxn.kf_A * exp(-u.reactions.reforming_reactions.WGS_rxn.kf_Ea / (R_gas * u.temp[cell_id]))
    K_eq_WGS = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.WGS_rxn)
    driving_force = 1.0 - ((P_CO2 * P_H2) / (K_eq_WGS * P_CO * P_H2O))

    u.net_rates.reforming_reactions.WGS_rxn = (k_WGS * K_OH * ((P_CO * P_H2O) / sqrt(P_H2)) * driving_force) / (DEN^2)

    for species_name in propertynames(u.molar_concentrations)
        for reforming_reaction in propertynames(u.reactions.reforming_reactions)
            du.molar_concentrations[species_name][cell_id] += u.net_rates.reforming_reactions[reforming_reaction] * u.reactions.reforming_reactions[reforming_reaction].all_stoich_coeffs[species_name] #this causes GC
            #hold up, we need to decide if we should do u.net_rates.reactions.reforming_reactions[reaction_name] OR reforming_reaction.net_rate OR u.reactions.net_rates.reforming_reaction
            #the disadvantage of the second one is that we start mixing up cached values (net_rates) and fixed property values (u.reforming_reactions)
        end
        du.mass_fractions[species_name][cell_id] += (du.molar_concentrations[species_name][cell_id] * u.molecular_weights[species_name]) / u.rho[cell_id]
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
    end

    for reforming_reaction in propertynames(u.reactions.reforming_reactions)
        du.heat[cell_id] += u.net_rates.reforming_reactions[reforming_reaction] * (-u.reactions.reforming_reactions[reforming_reaction].heat_of_reaction) * vol #this causes GC
        # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
    end
end

#PAM_reforming_react_cell!(du_tracer, u_tracer, 1, 1)

function cap_heat_flux_to_temp_change!(
    du, u,
    cell_id,
    vol
)
    # J/s /= m^3 * kg*m^3 * J/(kg*K)
    # = K/s
    du.temp[cell_id] += du.heat[cell_id] / (vol * u.rho[cell_id] * u.cp[cell_id])
end

#cap_heat_flux_to_temp_change!(du_tracer, u_tracer, 1, 1)

#turns rho changes to pressure changes
function cap_mass_flux_to_pressure_change!(
    du, u,
    cell_id,
    vol
)
    # kg/s /= (m^3 / (J/(mol*K) * K))
    #remember, J = Pa*m^3
    # = Pa/s
    du.pressure[cell_id] += du.mass[cell_id] / (vol / (R_gas * u.temp[cell_id]))
end

#cap_mass_flux_to_pressure_change!(du_tracer, u_tracer, 1, 1)

function reforming_area!(du, u, cell_id, vol)
    #property updating/retrieval
    #I wonder if we could make these automatic by inferring which variables are caches that need to be updated, but this is fine for now
    molar_concentrations!(u, cell_id)
    mw_avg!(u, cell_id)
    rho_ideal!(u, cell_id)

    #internal physics

    #we have to figure out if we're going to pass in just a singular cell_volume or all cell_volumes and use cell_id
    #power_law_react_cell!(du, u, cell_id, u.reactions.example_reaction, vol)
    #example of how to do a power law reaction

    PAM_reforming_react_cell!(du, u, cell_id, vol)

    #sources

    #boundary conditions

    #capacities
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
end

R_gas = 8.314
reforming_area!(du_tracer, u_tracer, 1, 1.0)

function test_access_logs(du, u, cell_id, vol)
    for reforming_reaction in propertynames(u.reactions.reforming_reactions)
        du.heat[cell_id] += u.net_rates.reforming_reactions[reforming_reaction] * (-u.reactions.reforming_reactions[reforming_reaction].heat_of_reaction) * vol
    end
end

#test_access_logs(du_tracer, u_tracer, 1, 1.0)

du_tracer
ctx_du = getfield(du_tracer, :tracer_ref)[]
ctx_du.access_logs #this is the only thing we're after

ctx_du.access_logs
#=
for cell_id in 1:10
    reforming_area!(du_tracer, u_tracer, cell_id, 1.0)
end
=#

mutable struct VarAccessLog
    path::Vector{Symbol}
    accesses::Vector{Symbol}
    n_du_reads::Int
    n_du_writes::Int
    n_u_reads::Int
    n_u_writes::Int
end

ctx_du.access_logs

function merge_trace_results(region_contexts)
    var_access_logs = Dict{Vector{Symbol}, VarAccessLog}()

    encountered_paths = Vector{Symbol}[]

    for var_access in region_contexts
        path = var_access.path

        access_type = var_access.access_type

        if !(path.path in encountered_paths) && path.root != :_expr
            push!(encountered_paths, path.path)
            push!(var_access_logs,
                path.path => VarAccessLog(
                    path.path,
                    Symbol[],
                    0, 0,
                    0, 0
                )
            )
        end

        if path.root == :du
            if access_type == :read
                var_access_logs[path.path].n_du_reads += 1
                push!(var_access_logs[path.path].accesses, :du_read)
            elseif access_type == :write
                var_access_logs[path.path].n_du_writes += 1
                push!(var_access_logs[path.path].accesses, :du_write)
            end
        elseif path.root == :u
            if access_type == :read
                var_access_logs[path.path].n_u_reads += 1
                push!(var_access_logs[path.path].accesses, :u_read)
            elseif access_type == :write
                var_access_logs[path.path].n_u_writes += 1
                push!(var_access_logs[path.path].accesses, :u_write)
            end
        end
    end

    return var_access_logs, encountered_paths
end

var_access_logs, encountered_paths = merge_trace_results(ctx_du.access_logs)

var_access_logs



collect(values(var_access_logs))

function classify_variables(var_access_logs)
    state_vars = Dict{Vector{Symbol},VarAccessLog}()
    cache_vars = Dict{Vector{Symbol},VarAccessLog}()
    fixed_vars = Dict{Vector{Symbol},VarAccessLog}()

    for (path, log) in var_access_logs
        accesses = log.accesses

        # Find the last du-related access
        last_du_access = nothing
        for i in length(accesses):-1:1
            if accesses[i] in (:du_read, :du_write)
                last_du_access = accesses[i]
                break
            end
        end

        if last_du_access == :du_write
            # State variable: the last thing we did to du was write to it
            state_vars[path] = log
        elseif log.n_u_writes > 0 || log.n_du_writes > 0 || log.n_du_reads > 0
            #TODO: this is a little risky because we might have a state variable that is stored in du, but it's unlikely
            # Cache: we wrote to u, then read it later
            cache_vars[path] = log
        else
            # Fixed: only ever read from u
            fixed_vars[path] = log
        end
    end

    return state_vars, cache_vars, fixed_vars
end

state_vars, cache_vars, fixed_vars = classify_variables(var_access_logs)

println("=== State Variables (written to du) ===")
for (path, _) in state_vars
    println("  ", join(path, "."))
end

println("\n=== Cache Variables (written to u, then read) ===")
for (path, _) in cache_vars
    println("  ", join(path, "."))
end

println("\n=== Fixed Parameters (only read from u) ===")
for (path, _) in fixed_vars
    println("  ", join(path, "."))
end

# ── ComponentArray Builder ──────────────────────────────────────────────────
# Takes classified paths and builds a nested ComponentVector
# :_idx in path → per-cell vector (zeros(n_cells))
# no :_idx → scalar (0.0)

using ComponentArrays

function inflate_path!(tree::Dict, clean_path::Vector{Symbol}, value)
    curr = tree
    for i in 1:(length(clean_path)-1)
        if !haskey(curr, clean_path[i])
            curr[clean_path[i]] = Dict{Symbol,Any}()
        end
        curr = curr[clean_path[i]]
    end
    curr[clean_path[end]] = value
end

function dict_to_component_vector(d::Dict{Symbol,Any}, n_cells::Int)
    pairs = Pair{Symbol,Any}[]
    for (k, v) in d
        if v isa Dict
            push!(pairs, k => dict_to_component_vector(v, n_cells))
        elseif v == :vector
            push!(pairs, k => zeros(n_cells))
        elseif v == :scalar
            push!(pairs, k => 0.0)
        end
    end
    return ComponentVector(; pairs...)
end

var_access_logs[[:temp, :_idx]]
var_access_logs[[:heat, :_idx]]
var_access_logs[[:rho, :_idx]]
var_access_logs[[:cp, :_idx]]
var_access_logs[[:mass, :_idx]]

function build_component_array(classified_vars::Dict{Vector{Symbol},VarAccessLog}, n_cells::Int)
    tree = Dict{Symbol,Any}()

    for (path, _) in classified_vars
        # Remove :_idx to get the clean nested path
        clean_path = filter(s -> s != :_idx, path)

        # :_idx means per-cell vector, otherwise scalar
        if :_idx in path
            template = :vector
        else
            template = :scalar
        end

        inflate_path!(tree, clean_path, template)
    end

    return dict_to_component_vector(tree, n_cells)
end

n_cells = 100

state_cv = build_component_array(state_vars, n_cells)
cache_cv = build_component_array(cache_vars, n_cells)
fixed_cv = build_component_array(fixed_vars, n_cells)

println("\n=== State ComponentVector ===")
println(Base.propertynames(state_cv))

println("\n=== Cache ComponentVector ===")
println(Base.propertynames(cache_cv))

println("\n=== Fixed ComponentVector ===")
println(Base.propertynames(fixed_cv))
=#

