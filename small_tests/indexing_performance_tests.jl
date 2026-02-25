
using BenchmarkTools
using ComponentArrays

function bench_named_tuple_pairs(mass_fractions, molar_concentrations)
    testing = 0.0
    for (key, value) in pairs(mass_fractions)
        testing += value + molar_concentrations[key]
    end
    return testing
end

function bench_named_tuple_keys(mass_fractions, molar_concentrations)
    testing = 0.0
    for key in keys(mass_fractions)
        testing += mass_fractions[key] + molar_concentrations[key]
    end
    return testing
end

function bench_named_tuple_map(nt1, nt2)
    return sum(map(+, nt1, nt2))
end

function bench_named_tuple_zip(nt1, nt2)
    testing = 0.0
    for (v1, v2) in zip(values(nt1), values(nt2))
        testing += v1 + v2
    end
    return testing
end

function bench_tuple_zip(t1, t2)
    testing = 0.0
    for (v1, v2) in zip(t1, t2)
        testing += v1 + v2
    end
    return testing
end

function bench_raw_indexing(t1, t2)
    testing = 0.0
    for i in eachindex(t1)
        testing += t1[i] + t2[i]
    end
    return testing
end

function bench_component_vector_propertynames(mass_fractions, molar_concentrations)
    testing = 0.0
    for property_name in propertynames(mass_fractions)
        testing += mass_fractions[property_name] + molar_concentrations[property_name]
    end
    return testing
end

function bench_component_vector_pairs(mass_fractions, molar_concentrations)
    testing = 0.0
    for (key, value) in pairs(mass_fractions)
        testing += mass_fractions[key] + molar_concentrations[key]
    end
    return testing
end

nt_mass = (methanol = 1.0, water = 1.3, hydrogen = 0.1)
nt_molar = (methanol = 0.4, water = 0.5, hydrogen = 0.05)

t_mass = (1.0, 1.3, 0.1)
t_molar = (0.4, 0.5, 0.05)
keys_to_idx = (methanol = 1, water = 2, hydrogen = 3)

cv_mass = ComponentVector(nt_mass)
cv_molar = ComponentVector(nt_molar)

println("NamedTuple pairs:")
@btime bench_named_tuple_pairs($nt_mass, $nt_molar)

println("NamedTuple keys:")
@btime bench_named_tuple_keys($nt_mass, $nt_molar)

println("NamedTuple map(+):")
@btime bench_named_tuple_map($nt_mass, $nt_molar)

println("NamedTuple zip(values):")
@btime bench_named_tuple_zip($nt_mass, $nt_molar)

println("Tuple zip:")
@btime bench_tuple_zip($t_mass, $t_molar)

println("Raw Tuple indexing:")
@btime bench_raw_indexing($t_mass, $t_molar)

println("ComponentVector propertynames:")
@btime bench_component_vector_propertynames($cv_mass, $cv_molar)

println("ComponentVector pairs:")
@btime bench_component_vector_pairs($cv_mass, $cv_molar)

