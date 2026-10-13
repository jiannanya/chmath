#include "test.hpp"
#include <chmath/numeric/decomposition.hpp>

using namespace chm;
namespace numeric_hotpaths {
template <floating T, std::size_t R, std::size_t C> void fitting() {
    const T tol = std::same_as<T, float> ? T(2e-4) : T(2e-12);
    for (int trial = 0; trial < 20; ++trial) {
        mat<T, R, C> a;
        vec<T, C> expected;
        vec<T, R> noise;
        for (T &v : a.elements)
            v = test::sample<T>(T(-1), T(1));
        for (T &v : expected.elements)
            v = test::sample<T>();
        for (T &v : noise.elements)
            v = test::sample<T>(T(-0.01), T(0.01));
        const auto b = a * expected;
        const auto x = least_squares(a, b);
        CHECK(x);
        NEAR(*x, expected, tol);
        const auto fit = least_squares(a, b + noise);
        CHECK(fit);
        // Orthogonality of the residual characterizes the least-squares solution.
        NEAR(transpose(a) * (a * *fit - b - noise), vec<T, C>{}, tol * T(R));
        const auto qr = factor_qr(a);
        CHECK(qr);
        NEAR(qr->q * qr->r, a, tol);
        NEAR(transpose(qr->q) * qr->q, (mat<T, R, R>::identity()), tol);
    }
}
template <floating T, std::size_t N> void lu_residuals() {
    const T tol = std::same_as<T, float> ? T(3e-4) : T(3e-12);
    for (int trial = 0; trial < 20; ++trial) {
        mat<T, N, N> a;
        for (T &v : a.elements)
            v = test::sample<T>(T(-1), T(1));
        for (std::size_t i = 0; i < N; ++i)
            a(i, i) += T(N);
        // Force row interchanges as well as nontrivial elimination.
        for (std::size_t j = 0; j < N; ++j)
            std::swap(a(0, j), a(N - 1, j));
        const auto lu = factor_lu(a);
        CHECK(lu);
        mat<T, N, N> l = mat<T, N, N>::identity(), u, pa;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                (i > j ? l(i, j) : u(i, j)) = lu->factors(i, j);
                pa(i, j) = a(lu->permutation[i], j);
            }
        NEAR(l * u, pa, tol);
        NEAR(a * lu->solve(mat<T, N, N>::identity()), (mat<T, N, N>::identity()), tol);
    }
}
} // namespace numeric_hotpaths

CH_TEST(numeric_hotpath_least_squares_float_tall) {
    numeric_hotpaths::fitting<float, 8, 3>();
    numeric_hotpaths::fitting<float, 32, 4>();
}
CH_TEST(numeric_hotpath_least_squares_double_tall) {
    numeric_hotpaths::fitting<double, 32, 4>();
    numeric_hotpaths::fitting<double, 64, 8>();
}
CH_TEST(numeric_hotpath_least_squares_square_and_singleton) {
    numeric_hotpaths::fitting<double, 1, 1>();
    numeric_hotpaths::fitting<double, 8, 8>();
}
CH_TEST(numeric_hotpath_least_squares_scaling_and_failures) {
    const auto a =
        mat<double, 4, 2>::from_rows({vec2d{1, 0}, vec2d{1, 1}, vec2d{1, 2}, vec2d{1, 3}});
    for (double scale : {1e-200, 1.0, 1e200}) {
        const auto x = least_squares(a * scale, vec4d{1, 3, 5, 7} * scale);
        CHECK(x);
        NEAR(*x, vec2d{1, 2}, 1e-12);
    }
    const double huge = std::numeric_limits<double>::max() * 0.8;
    const auto x = least_squares(mat2d::identity(), vec2d{huge, -huge});
    CHECK(x);
    NEAR(*x / huge, vec2d{1, -1}, 1e-12);
    const auto projected = least_squares(
        mat<double, 2, 1>::from_rows({vec<double, 1>{0}, vec<double, 1>{1}}), vec2d{huge, huge});
    CHECK(projected);
    NEAR((*projected)[0] / huge, 1.0, 1e-12);
    auto deficient = a;
    deficient.set_column(1, a.column(0));
    CHECK(!least_squares(deficient, vec4d{1}));
    deficient.set_column(0, vec4d{});
    CHECK(!least_squares(deficient, vec4d{1}));
    CHECK(!least_squares(a, vec4d{1}, -1.0));
    CHECK(!least_squares(a, vec4d{1}, std::numeric_limits<double>::infinity()));
    CHECK(!least_squares(a, vec4d{std::numeric_limits<double>::quiet_NaN()}));
    deficient(0, 0) = std::numeric_limits<double>::infinity();
    CHECK(!least_squares(deficient, vec4d{1}));
    auto threshold = mat<double, 8, 2>{};
    threshold(0, 0) = 1;
    threshold(1, 1) = 1e-8;
    CHECK(!least_squares(threshold, vec<double, 8>{1}, 1e-7));
    CHECK(least_squares(threshold, vec<double, 8>{1}, 1e-9));
}
CH_TEST(numeric_hotpath_lu_large_float) {
    numeric_hotpaths::lu_residuals<float, 5>();
    numeric_hotpaths::lu_residuals<float, 16>();
}
CH_TEST(numeric_hotpath_lu_large_double) {
    numeric_hotpaths::lu_residuals<double, 8>();
    numeric_hotpaths::lu_residuals<double, 32>();
}
CH_TEST(numeric_hotpath_lu_large_tiny_and_constexpr) {
    constexpr auto lu = factor_lu(mat<double, 16, 16>::identity());
    static_assert(lu && lu->solve(vec<double, 16>{1}) == vec<double, 16>{1});
    const auto a = mat<double, 16, 16>::identity() * 1e-310;
    const auto tiny = factor_lu(a);
    CHECK(tiny);
    NEAR(tiny->solve(vec<double, 16>{1e-310}), (vec<double, 16>{1}), 1e-12);
}
