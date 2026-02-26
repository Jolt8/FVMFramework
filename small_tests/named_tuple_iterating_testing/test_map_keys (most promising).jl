using BenchmarkTools

function test_map_keys(named_tuple, other_thing, third_thing)
    # If the user wants to access `other_thing` and `third_thing` via `name_sym`
    map(keys(named_tuple), values(named_tuple)) do name_sym, species_arr
        other_arr = getproperty(other_thing, name_sym)
        third_arr = getproperty(third_thing, name_sym)

        for i in eachindex(species_arr)
            species_arr[i] += other_arr[i] + third_arr[i]
        end
    end
end

v1 = rand(100);
v2 = rand(100);
v3 = rand(100);
axes = (c=1:50, d=51:100)

named_tuple = (c=view(v1, axes.c), d=view(v1, axes.d))
other_tuple = (c=view(v2, axes.c), d=view(v2, axes.d))
third_tuple = (c=view(v3, axes.c), d=view(v3, axes.d))

@btime test_map_keys($named_tuple, $other_tuple, $third_tuple)
