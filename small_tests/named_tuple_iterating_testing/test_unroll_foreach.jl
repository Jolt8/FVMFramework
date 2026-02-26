using BenchmarkTools

function test_unroll_foreach(named_tuple, other_thing)
    # Using foreach with Val over keys
    # What if we just iterate over properties?
    # but the user needs the name!

    # foreach over a tuple of keys:
    # keys(nt) returns a tuple of symbols: (:c, :d)
    # we can do Val.(keys(nt)) to get (Val{:c}(), Val{:d}())

    foreach(Val.(keys(named_tuple))) do val_species
        species_name = typeof(val_species).parameters[1]
        species_arr = getproperty(named_tuple, species_name)
        other_arr = getproperty(other_thing, species_name)

        for i in eachindex(species_arr)
            species_arr[i] += other_arr[i]
        end
    end
end

v1 = rand(100)
v2 = rand(100)
axes = (c=1:50, d=51:100)

named_tuple = (c=view(v1, axes.c), d=view(v1, axes.d))
other_tuple = (c=view(v2, axes.c), d=view(v2, axes.d))


@btime test_unroll_foreach($named_tuple, $other_tuple)

