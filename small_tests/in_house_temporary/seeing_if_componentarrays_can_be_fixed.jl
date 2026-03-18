using ComponentArrays

@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{ComponentVector, N}) where {N}
    first_ax = groups[1].parameters[4].parameters[1]
    
    #the first parameter passed in is assumed to be the NamedTuple containing the layout
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

n_cells = 100

u = ComponentArray(
    mass_fractions = (
        methylene_blue = zeros(n_cells),
        water = zeros(n_cells),
    ),
    molar_concentrations = (
        methylene_blue = zeros(n_cells), 
        water = zeros(n_cells)
    )
)

function test_foreach_field_allocations(u, n_cells)
    for cell_id in 1:n_cells
        foreach_field_at!(cell_id, u.mass_fractions, u.molar_concentrations) do species, mass_fractions, molar_concentrations
            mass_fractions[species] += 1.0
            molar_concentrations[species] += 1.0
        end
    end
    return nothing
end

@btime test_foreach_field_allocations($u, $n_cells) 