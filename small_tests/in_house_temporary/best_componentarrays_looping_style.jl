using ComponentArrays
using BenchmarkTools

@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{ComponentVector, N}) where {N}
    first_ax = groups[1].parameters[4].parameters[1]
    
    #the first parameter passed in is assumed to contain the same axes of the other groups passed in
    properties = keys(first_ax.parameters[1]) 
    
    exprs = []
    for n in properties
        view_exprs = [:(getproperty(groups[$i], $(QuoteNode(n)))) for i in 1:N]
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end
    
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

@generated function foreach_field_at!(f, groups::Vararg{ComponentVector, N}) where {N}
    first_ax = groups[1].parameters[4].parameters[1]
    
    #the first parameter passed in is assumed to contain the same axes of the other groups passed in
    layout = first_ax.parameters[1]
    
    exprs = []
    for n in layout
        push!(exprs, :(f($n, groups...)))
    end
    
    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

function create_component_array_of_size(n)
    return ComponentArray(
        mass_fractions = (
            methylene_blue = zeros(n),
            water = zeros(n),
        ),
        molar_concentrations = (
            methylene_blue = zeros(n), 
            water = zeros(n)
        ),
        net_rates = (
            MSR_rxn = 0.0,
            WGS_rxn = 0.0
        )
    )
end

function test_foreach_field_allocations(u, n_cells)
    for cell_id in 1:n_cells
        foreach_field_at!(cell_id, u.mass_fractions, u.molar_concentrations) do species, mass_fractions, molar_concentrations
            mass_fractions[species] += 1.0
            molar_concentrations[species] += 1.0
        end
    end
    foreach_field_at!(u.net_rates) do rxn, net_rates
        net_rates[rxn] += 1.0
    end
    return nothing
end

function test_raw_field_access(u, n_cells)
    for cell_id in 1:n_cells
        u.mass_fractions.methylene_blue[cell_id] += 1.0
        u.mass_fractions.water[cell_id] += 1.0
        u.molar_concentrations.methylene_blue[cell_id] += 1.0
        u.molar_concentrations.water[cell_id] += 1.0
    end
    u.net_rates.MSR_rxn += 1.0
    u.net_rates.WGS_rxn += 1.0
    return nothing
end

function compare_raw(u_vec, n_cells)
    for cell_id in 1:n_cells
        u_vec[cell_id] += 1.0
    end
end

n_cells = 10000000

u = create_component_array_of_size(n_cells)

@btime test_foreach_field_allocations($u, $n_cells) 
#165 ns for 100 cells
#15.4 μs for 10000 cells
#31.4 ms for 10000000 cells
#include("different_style_componentarrays_looping.jl") is actually faster than this somehow, 
#this is weird because it means that it's also faster than the raw field access

@btime test_raw_field_access($u, $n_cells) 
#same as above

u_vec = Vector(u)
relative_size = length(u_vec)
@btime compare_raw($u_vec, $relative_size) 
#34 ns for 100 cells
#6.9 μs for 10000 cells
#26.529 ms for 10000000 cells 