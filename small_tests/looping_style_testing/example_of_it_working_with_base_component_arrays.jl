using ComponentArrays
using BenchmarkTools

@generated function foreach_field_at!(f, groups::Vararg{ComponentVector, N}) where {N}
    # For a ComponentVector, the axes are in the 4th type parameter
    AxTuple = groups[1].parameters[4]
    
    # Extract the first Axis
    Ax = AxTuple.parameters[1]
    
    # The first parameter of Axis is the NamedTuple containing the layout
    layout = Ax.parameters[1]
    properties = keys(layout)
    
    exprs = []
    offset = 0
    for (i, n) in enumerate(properties)
        layout_val = layout[i]
        if layout_val isa ViewAxis
            # Nested arrays: pass cell_id as index and the sub-fields as containers
            args = [:(view(groups[$j], :)) for j in 1:N]
            arr_length = length(layout_val.ax)
            push!(exprs, :(f($offset+1:$(offset+arr_length), $(args...))))
            offset += arr_length
        else
            # Flat fields: pass property index/axis and the groups themselves as containers
            push!(exprs, :(f($layout_val, groups...)))
        end
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
    molecular_weights = (
        methylene_blue = 271.18,
        water = 18.015
    ),
    net_rates = (
        WGS_rxn = 1.0,
        MD_rxn = 2.0,
    )
)

function test_foreach_field(u)
    for cell_id in 1:100
        foreach_field_at!(u.mass_fractions, u.molar_concentrations, u.molecular_weights) do species, mass_fractions, molar_concentrations, molecular_weights
            mass_fractions[species[cell_id]] += 1.0
            molar_concentrations[species[cell_id]] += 1.0
            u.net_rates.WGS_rxn += 1.0
            #molecular_weights[species] += 1.0 #this doesn't work and throws and error, 
            #if this is absolutely impossible to implement, then it can be left out

            #I kinda have a feeling that it is impossible to implement because species is returned as a 
            #UnitRange, and returning an Int or UnitRange depending on the type of the field is probably not possible
        end
    end

    foreach_field_at!(u.net_rates) do rxn_idx, net_rates
        net_rates[rxn_idx] += 1.0
    end
    return nothing
end

test_foreach_field(u)

@btime test_foreach_field($u) #393.069 ns, 0 allocations! 
