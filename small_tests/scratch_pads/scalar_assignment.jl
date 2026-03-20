using ComponentArrays

test_1 = ComponentVector(temp = [1.0])
test_2 = ComponentVector(temp = [2.0])

for property_name in propertynames(test_1)
    test_1_val = getproperty(test_1, property_name)
    test_2_val = getproperty(test_2, property_name)

    println(test_1_val)
    println(test_2_val)
    println(typeof(test_1_val))

    test_1[property_name] = test_2_val
end

test_1