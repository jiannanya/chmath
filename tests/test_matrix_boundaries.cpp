#include "boundary_helpers.hpp"
#include "test.hpp"
using namespace chm;

namespace matrix_boundary {
template <class T, std::size_t N> void inverse_check() {
    const auto a = boundary::well_conditioned<T, N>();
    const auto inv = inverse(a);
    CHECK(inv);
    NEAR(a * (*inv), mat<T, N, N>::identity(), T(3e-5));
    const auto b = a * boundary::sequence<T, N>();
    const auto x = solve(a, b);
    CHECK(x);
    NEAR(*x, boundary::sequence<T, N>(), T(3e-5));
}
} // namespace matrix_boundary
CH_TEST(matrix_inverse_dimension_one) {
    CHECK(inverse(mat<double, 1, 1>{4})->elements[0] == 0.25);
}
CH_TEST(matrix_inverse_dimension_two) {
    matrix_boundary::inverse_check<double, 2>();
}
CH_TEST(matrix_inverse_dimension_three) {
    matrix_boundary::inverse_check<double, 3>();
}
CH_TEST(matrix_inverse_dimension_four) {
    matrix_boundary::inverse_check<double, 4>();
}
CH_TEST(matrix_inverse_dimension_five) {
    matrix_boundary::inverse_check<double, 5>();
}
CH_TEST(matrix_inverse_dimension_eight_float) {
    matrix_boundary::inverse_check<float, 8>();
}
CH_TEST(matrix_rectangular_product_reference) {
    mat<double, 2, 3> a;
    mat<double, 3, 4> b;
    for (std::size_t i = 0; i < a.elements.size(); ++i)
        a.elements[i] = double(i + 1);
    for (std::size_t i = 0; i < b.elements.size(); ++i)
        b.elements[i] = double(i + 1) / 7;
    const auto p = a * b;
    for (std::size_t r = 0; r < 2; ++r)
        for (std::size_t c = 0; c < 4; ++c) {
            long double expected = 0;
            for (std::size_t k = 0; k < 3; ++k)
                expected += static_cast<long double>(a(r, k)) * b(k, c);
            NEAR(p(r, c), static_cast<double>(expected), 1e-13);
        }
}
CH_TEST(matrix_rectangular_identity_diagonal) {
    const mat<int, 2, 4> a{3};
    CHECK(a(0, 0) == 3 && a(1, 1) == 3 && a(0, 3) == 0);
}
CH_TEST(matrix_determinant_balanced_extreme_diagonal) {
    auto a = mat4d::identity();
    a(0, 0) = a(1, 1) = 1e200;
    a(2, 2) = a(3, 3) = 1e-200;
    NEAR(determinant(a), 1.0, 1e-13);
}
CH_TEST(matrix_determinant_balanced_underflow_diagonal) {
    // Include gradual underflow, where the intermediate is nonzero but rounded.
    for (double small : {1e-200, 1e-160, 1e-155}) {
        auto a = mat4d::identity();
        a(0, 0) = a(1, 1) = small;
        a(2, 2) = a(3, 3) = 1 / small;
        NEAR(determinant(a), 1.0, 1e-13);
    }
}
CH_TEST(matrix_determinant_row_exchange_sign) {
    auto a = mat3d::identity();
    a.set_column(0, vec3d{0, 1, 0});
    a.set_column(1, vec3d{1, 0, 0});
    CHECK(determinant(a) == -1.0);
}
CH_TEST(matrix_determinant_cyclic_permutation) {
    const auto a = mat3d::from_rows({vec3d{0, 1, 0}, vec3d{0, 0, 1}, vec3d{1, 0, 0}});
    CHECK(determinant(a) == 1.0);
}
CH_TEST(matrix_determinant_overflow_not_false_zero) {
    const auto a = mat2d::from_rows({vec2d{boundary::maximum, boundary::maximum},
                                     vec2d{-boundary::maximum, boundary::maximum}});
    CHECK(std::isinf(determinant(a)));
    CHECK(determinant(a) > 0);
}
CH_TEST(matrix_singular_duplicate_rows) {
    const auto a = mat3d::from_rows({vec3d{1, 2, 3}, vec3d{4, 5, 6}, vec3d{1, 2, 3}});
    CHECK(!inverse(a));
    CHECK(determinant(a) == 0);
}
CH_TEST(matrix_singular_zero_column) {
    auto a = mat4d::identity();
    a.set_column(2, vec4d{});
    CHECK(!factor_lu(a));
}
CH_TEST(matrix_lu_negative_tolerance) {
    CHECK(!factor_lu(mat2d::identity(), -0.01));
}
CH_TEST(matrix_lu_nan_tolerance) {
    CHECK(!factor_lu(mat2d::identity(), boundary::nan));
}
CH_TEST(matrix_lu_infinite_tolerance) {
    CHECK(!factor_lu(mat2d::identity(), boundary::inf));
}
CH_TEST(matrix_lu_explicit_zero_tolerance) {
    auto a = mat2d::from_rows({vec2d{1, 1}, vec2d{1, 1 + 1e-15}});
    CHECK(!factor_lu(a));
    CHECK(factor_lu(a, 0.0));
}
CH_TEST(matrix_lu_many_rhs_reuse) {
    const auto a = boundary::well_conditioned<double, 4>();
    const auto lu = factor_lu(a);
    CHECK(lu);
    for (int i = 0; i < 32; ++i) {
        const vec4d x{double(i), 1, 2, 3};
        NEAR(lu->solve(a * x), x, 1e-12);
    }
}
CH_TEST(matrix_solve_nonfinite_rhs) {
    CHECK(!solve(mat3d::identity(), vec3d{0, boundary::inf, 0}));
    CHECK(!solve(mat3d::identity(), vec3d{0, 0, boundary::nan}));
}
CH_TEST(matrix_inverse_input_unchanged) {
    const auto original = boundary::well_conditioned<double, 4>();
    auto a = original;
    CHECK(inverse(a));
    CHECK(a == original);
}
CH_TEST(matrix_self_multiplication_alias) {
    auto a = mat2d::from_rows({vec2d{1, 2}, vec2d{3, 4}});
    a = a * a;
    CHECK(a == mat2d::from_rows({vec2d{7, 10}, vec2d{15, 22}}));
}
CH_TEST(matrix_transpose_involution_nonsquare) {
    mat<int, 3, 7> a;
    for (std::size_t i = 0; i < a.elements.size(); ++i)
        a.elements[i] = int(i);
    CHECK(transpose(transpose(a)) == a);
}
CH_TEST(matrix_outer_product_rank_one) {
    const auto a = outer_product(vec3d{1, 2, 3}, vec3d{4, 5, 6});
    CHECK(!inverse(a));
    CHECK(a * vec3d{5, -4, 0} == vec3d{});
}
CH_TEST(matrix_cholesky_asymmetric_rejection) {
    auto a = mat3d{2};
    a(1, 0) = 0.1;
    CHECK(!cholesky(a));
}
CH_TEST(matrix_cholesky_semidefinite_rejection) {
    auto a = mat3d::identity();
    a(1, 1) = 0;
    CHECK(!cholesky(a));
}
CH_TEST(matrix_cholesky_negative_pivot) {
    auto a = mat3d::identity();
    a(2, 2) = -1;
    CHECK(!cholesky(a));
}
CH_TEST(matrix_qr_tall_reconstruction) {
    mat<double, 7, 3> a;
    for (std::size_t i = 0; i < 7; ++i) {
        a(i, 0) = 1;
        a(i, 1) = double(i);
        a(i, 2) = double(i * i);
    }
    const auto qr = factor_qr(a);
    CHECK(qr);
    NEAR(qr->q * qr->r, a, 1e-12);
    NEAR(transpose(qr->q) * qr->q, mat<double, 7, 7>::identity(), 1e-12);
}
CH_TEST(matrix_qr_negative_leading_column) {
    const auto a = mat2d::from_rows({vec2d{-3, 1}, vec2d{-4, 2}});
    const auto qr = factor_qr(a);
    CHECK(qr);
    NEAR(qr->q * qr->r, a, 1e-12);
}
CH_TEST(matrix_least_squares_orthogonal_residual) {
    const auto a = mat<double, 3, 2>::from_rows({vec2d{1, 0}, vec2d{0, 1}, vec2d{1, 1}});
    const vec3d b{1, 2, 4};
    const auto x = least_squares(a, b);
    CHECK(x);
    NEAR(transpose(a) * (a * (*x) - b), vec2d{}, 1e-12);
}
CH_TEST(matrix_least_squares_rank_deficiency) {
    const auto a = mat<double, 3, 2>::from_rows({vec2d{1, 2}, vec2d{2, 4}, vec2d{3, 6}});
    CHECK(!least_squares(a, vec3d{1, 2, 3}));
}
CH_TEST(matrix_eigen_repeated_values) {
    const auto e = symmetric_eigen(mat3d{7});
    CHECK(e && e->converged);
    CHECK(e->values == vec3d{7});
    NEAR(transpose(e->vectors) * e->vectors, mat3d::identity(), 1e-14);
}
CH_TEST(matrix_eigen_indefinite_sorted) {
    auto a = mat3d::identity();
    a(0, 0) = -3;
    a(1, 1) = 2;
    const auto e = symmetric_eigen(a);
    CHECK(e && e->converged);
    CHECK(e->values == vec3d{-3, 1, 2});
}
CH_TEST(matrix_eigen_zero_budget_reports_failure) {
    auto a = mat3d{2};
    a(0, 1) = a(1, 0) = 0.5;
    const auto e = symmetric_eigen(a, 1e-12, 0);
    CHECK(e && !e->converged && e->sweeps == 0);
}
CH_TEST(matrix_vector_product_reference) {
    // Covers the 4x4 float SIMD kernel, the generic accumulator path, and a
    // rectangular shape against a long-double oracle.
    mat4f a;
    for (std::size_t i = 0; i < a.elements.size(); ++i)
        a.elements[i] = float(i % 7) - 3.0F;
    const vec4f v{-1.5F, 0.25F, 2.0F, -0.75F};
    const auto r = a * v;
    for (std::size_t i = 0; i < 4; ++i) {
        long double expected = 0;
        for (std::size_t j = 0; j < 4; ++j)
            expected += static_cast<long double>(a(i, j)) * static_cast<long double>(v[j]);
        NEAR(r[i], static_cast<float>(expected), 1e-5F);
    }
    const vec4d vd{-1.5, 0.25, 2.0, -0.75};
    mat4d b;
    for (std::size_t i = 0; i < b.elements.size(); ++i)
        b.elements[i] = double(i % 5) * 0.25 - 0.5;
    const auto rd = b * vd;
    for (std::size_t i = 0; i < 4; ++i) {
        long double expected = 0;
        for (std::size_t j = 0; j < 4; ++j)
            expected += static_cast<long double>(b(i, j)) * static_cast<long double>(vd[j]);
        NEAR(rd[i], static_cast<double>(expected), 1e-13);
    }
    mat<double, 3, 5> c;
    for (std::size_t i = 0; i < c.elements.size(); ++i)
        c.elements[i] = double(i) / 3.0;
    const vec<double, 5> v5{1, -2, 3, -4, 5};
    const auto r5 = c * v5;
    for (std::size_t i = 0; i < 3; ++i) {
        long double expected = 0;
        for (std::size_t j = 0; j < 5; ++j)
            expected += static_cast<long double>(c(i, j)) * static_cast<long double>(v5[j]);
        NEAR(r5[i], static_cast<double>(expected), 1e-13);
    }
    // Constant evaluation must not use the runtime kernel and must agree.
    constexpr mat4f constant_matrix{2.0F};
    constexpr vec4f constant_vector{1.0F, 2.0F, 3.0F, 4.0F};
    constexpr auto constant_result = constant_matrix * constant_vector;
    CHECK(constant_result == vec4f{2.0F, 4.0F, 6.0F, 8.0F});
}
