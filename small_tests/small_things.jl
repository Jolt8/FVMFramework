abstract type TestType end

struct TestType1 <: TestType end
struct TestType2 <: TestType end

vec = Vector{TestType}(undef, 5)

vec[1] = TestType1()
vec[2] = TestType2()
vec[3] = TestType1()
vec[4] = TestType2()
vec[5] = TestType1()

vec

struct TestGroup
    type::TestType
end

struct TestGroupParametric{T <: TestType}
    type::T
end

region_group = TestGroup(TestType1())

region_group.type

region_group = TestGroupParametric(TestType1())

region_group.type

