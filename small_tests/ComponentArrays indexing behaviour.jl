using ComponentArrays

#views are no necessary when indexing with Int
test = ComponentVector(mass_face = [[1, 2, 3], [4, 5, 6]])
saved_test = ComponentVector(mass_face = [[1, 2, 3], [4, 5, 6]])

cell_id = 1
face_idx = 2

test.mass_face[cell_id][face_idx] = 10

test == saved_test

#views are necessary when indexing a symbol that contains a vector

test_2 = ComponentVector(mass_fractions = (methanol = [1.0, 2.0], water = [3.0, 4.0]))
saved_test_2 = ComponentVector(mass_fractions = (methanol = [1.0, 2.0], water = [3.0, 4.0]))

test_2 == saved_test_2

test_2.mass_fractions[:methanol][1] = 10

test_2 == saved_test_2

view(test_2.mass_fractions, :methanol)[1] = 10

test_2 == saved_test_2

#views are not necessary when indexing a symbol that contains a scalar

test_3 = ComponentVector(mass_fractions = (methanol = 1.0, water = 2.0))
saved_test_3 = ComponentVector(mass_fractions = (methanol = 1.0, water = 2.0))

test_3 == saved_test_3

test_3.mass_fractions[:methanol] = 10

test_3 == saved_test_3

#view(test_3.mass_fractions, :methanol) = 10 this errors