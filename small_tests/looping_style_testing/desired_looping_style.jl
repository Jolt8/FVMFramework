#This is the ideal syntax:

#pretend u.mass_fractions is indexed by u.mass_fractions.methanol[cell_id]
#pretend molecular_weights are indexed like u.molecular_weights.methanol #(notice no cell_id because its a scalar for every cell)

for cell_id in 1:100
    for_fields!(u.mass_fractions, u.molecular_weights) do species, mass_fractions, molecular_weights
        mass_fractions[species[cell_id]] += 1.0
        #molecular_weights[species] += 1.0 #this is not 100% necessary, if it's impossible to get it to work with this, then that's fine, most properties are stored per cell anyways for consistency
        u.rho[cell_id] += 1.0 #stuff like this can still be indexed inside this function and not incur any allocations
    end
end

foreach_field_at!(u.net_rates) do rxn_idx, net_rates
    net_rates[rxn_idx] += 1.0
end