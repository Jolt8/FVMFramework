using BenchmarkTools

function test_pairs_foreach(named_tuple, other_thing)
    # Testing if pairs() unrolls
    foreach(pairs(named_tuple)) do (species_name, species_arr)
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

@btime test_pairs_foreach($named_tuple, $other_tuple)
