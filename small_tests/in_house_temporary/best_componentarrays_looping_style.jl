using ComponentArrays

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

u = ComponentArray(
    mass_fractions = (
        methylene_blue = zeros(100),
        water = zeros(100),
    ),
    molar_concentrations = (
        methylene_blue = zeros(100), 
        water = zeros(100)
    ),
    net_rates = (
        MSR_rxn = 0.0,
        WGS_rxn = 0.0
    )
)

function test_foreach_field_allocations(u)
    for cell_id in 1:100
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

@btime test_foreach_field_allocations($u) 