
function test_get_tmp(du, u)
    println(typeof(du))
    test = get_tmp(DiffCache([0.0, 0.0], 10), du)

    println(test)
end

jac_sparsity = ADTypes.jacobian_sparsity(
    test_get_tmp, du_vec, u_vec, TracerSparsityDetector()
)
