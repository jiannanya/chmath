#include "test.hpp"
#include <chmath/numeric/decomposition.hpp>
#include <chmath/numeric/numeric.hpp>

using namespace chm;

CH_TEST(quadratic_stability_and_degeneracy) {
    const auto roots = solve_quadratic(1.0, -3.0, 2.0);
    CHECK(roots);
    CHECK(roots->count == 2);
    NEAR(roots->values[0], 1.0, 1e-14);
    NEAR(roots->values[1], 2.0, 1e-14);
    const auto stable = solve_quadratic(1.0, 1e8, 1.0);
    CHECK(stable);
    CHECK(stable->count == 2);
    NEAR(stable->values[1], -1e-8, 1e-20);
    NEAR(stable->values[0], -1e8, 1e-14);
    CHECK(solve_quadratic(1.0, 0.0, 1.0)->count == 0);
    CHECK(solve_quadratic(1.0, 2.0, 1.0)->count == 1);
    CHECK(solve_quadratic(0.0, 2.0, 4.0)->values[0] == -2.0);
    CHECK(solve_quadratic(0.0, 0.0, 1.0)->count == 0);
    CHECK(solve_quadratic(0.0, 0.0, 0.0)->all_reals);
    CHECK(!solve_quadratic(std::numeric_limits<double>::infinity(), 1.0, 2.0));
    CHECK(solve_quadratic(1.0, 0.0, 0.0)->values[0] == 0.0);
    for (double scale : {1e-200, 1.0, 1e200}) {
        const auto r = solve_quadratic(scale, -3.0 * scale, 2.0 * scale);
        CHECK(r);
        CHECK(r->count == 2);
        NEAR(r->values[0], 1.0, 1e-14);
        NEAR(r->values[1], 2.0, 1e-14);
    }
    for (int i = 0; i < 1000; ++i) {
        const double a = test::sample<double>(-10, -1), b = test::sample<double>(1, 10);
        const auto r = solve_quadratic(1.0, -(a + b), a * b);
        CHECK(r);
        CHECK(r->count == 2);
        NEAR(r->values[0], a, 1e-12);
        NEAR(r->values[1], b, 1e-12);
    }
}
CH_TEST(polynomials_roots_and_sum) {
    const std::array<double, 3> coefficients{1, 2, 3};
    CHECK(evaluate_polynomial<double>(coefficients, 2.0) == 17.0);
    CHECK(evaluate_polynomial<double>({}, 2.0) == 0.0);
    const auto root = bisect([](double x) { return x * x - 2; }, 0.0, 2.0, 1e-13, 1e-13);
    CHECK(root && root->converged);
    NEAR(root->value, std::sqrt(2.0), 1e-12);
    const auto endpoint = bisect([](double x) { return x; }, 0.0, 2.0);
    CHECK(endpoint && endpoint->converged && endpoint->iterations == 0);
    CHECK(!bisect([](double x) { return x * x + 1; }, -1.0, 1.0));
    CHECK(!bisect([](double x) { return x; }, 2.0, 1.0));
    CHECK(!bisect([](double) { return std::numeric_limits<double>::quiet_NaN(); }, 0.0, 1.0));
    const auto limited = bisect([](double x) { return x * x - 2; }, 0.0, 2.0, 0.0, 0.0, 1);
    CHECK(limited && !limited->converged);
    compensated_sum<double> sum;
    sum.add(1e16);
    sum.add(1.0);
    sum.add(-1e16);
    CHECK(sum.value() == 1.0);
    sum.reset();
    CHECK(sum.value() == 0.0);
    sum.add(1);
    sum.add(1e16);
    sum.add(-1e16);
    CHECK(sum.value() == 1.0);
}
CH_TEST(adaptive_integration) {
    const auto cubic = integrate([](double x) { return x * x * x; }, 0.0, 1.0, 1e-12);
    CHECK(cubic && cubic->converged);
    NEAR(cubic->value, 0.25, 1e-13);
    const auto sine = integrate([](double x) { return std::sin(x); }, 0.0, pi<double>, 1e-10);
    CHECK(sine && sine->converged);
    NEAR(sine->value, 2.0, 1e-10);
    CHECK(sine->evaluations > 5);
    const auto backwards = integrate([](double x) { return x * x; }, 1.0, 0.0, 1e-12);
    CHECK(backwards);
    NEAR(backwards->value, -1.0 / 3.0, 1e-12);
    const auto zero = integrate([](double) { return 9.0; }, 1.0, 1.0);
    CHECK(zero && zero->value == 0 && zero->evaluations == 0);
    const auto limited = integrate([](double x) { return std::exp(x); }, 0.0, 1.0, 1e-14, 20, 5);
    CHECK(limited && !limited->converged && limited->evaluations <= 5);
    const auto depth = integrate([](double x) { return std::exp(x); }, 0.0, 1.0, 1e-14, 0);
    CHECK(depth && !depth->converged);
    CHECK(!integrate([](double) { return std::numeric_limits<double>::infinity(); }, 0.0, 1.0));
    CHECK(!integrate([](double x) { return x; }, 0.0, 1.0, -1.0));
    CHECK(!integrate([](double x) { return x; }, 0.0, 1.0, 1e-6, 33));
}
template <class T, std::size_t N> void decomposition_properties() {
    const T tol = T(3e-4);
    for (int i = 0; i < 150; ++i) {
        mat<T, N, N> a;
        for (T &x : a.elements)
            x = test::sample<T>(T(-1), T(1));
        const auto positive = transpose(a) * a + mat<T, N, N>::identity();
        const auto l = cholesky(positive);
        CHECK(l);
        NEAR((*l) * transpose(*l), positive, tol);
        const auto qr = factor_qr(a);
        CHECK(qr);
        NEAR(qr->q * qr->r, a, tol);
        NEAR(transpose(qr->q) * qr->q, mat<T, N, N>::identity(), tol);
        const auto eigen = symmetric_eigen(positive);
        CHECK(eigen && eigen->converged);
        NEAR(transpose(eigen->vectors) * eigen->vectors, mat<T, N, N>::identity(), tol);
        for (std::size_t j = 0; j < N; ++j) {
            NEAR(positive * eigen->vectors.column(j), eigen->vectors.column(j) * eigen->values[j],
                 tol);
            if (j > 0)
                CHECK(eigen->values[j] >= eigen->values[j - 1]);
        }
        NEAR(dot(eigen->values, vec<T, N>(T(1))), trace(positive), tol);
    }
}
CH_TEST(decomposition_random_properties) {
    decomposition_properties<float, 3>();
    decomposition_properties<double, 3>();
    decomposition_properties<double, 6>();
}
CH_TEST(qr_least_squares_and_failures) {
    const auto a =
        mat<double, 4, 2>::from_rows({vec2d{1, 0}, vec2d{1, 1}, vec2d{1, 2}, vec2d{1, 3}});
    const vec<double, 4> observations{1, 3, 5, 7};
    const auto x = least_squares(a, observations);
    CHECK(x);
    NEAR(*x, vec2d{1, 2}, 1e-12);
    const auto noisy = observations + vec4d{0.1, -0.2, 0.1, 0.2};
    const auto fit = least_squares(a, noisy);
    CHECK(fit);
    NEAR(transpose(a) * (a * (*fit) - noisy), vec2d{}, 1e-12);
    CHECK(!least_squares(mat<double, 4, 2>{}, observations));
    const auto deficient = mat<double, 3, 2>::from_rows({vec2d{1, 2}, vec2d{2, 4}, vec2d{3, 6}});
    CHECK(!least_squares(deficient, vec3d{1, 2, 3}));
    auto asym = mat3d::identity();
    asym(1, 0) = 1;
    CHECK(!cholesky(asym));
    CHECK(!symmetric_eigen(asym));
    CHECK(!cholesky(mat3d{-1}));
    CHECK(!cholesky(mat3d{}));
    CHECK(symmetric_eigen(mat3d{})->converged);
    CHECK(factor_qr(mat<double, 3, 2>{}));
    CHECK(symmetric_eigen(mat<double, 1, 1>{2})->values[0] == 2.0);
    auto coupled = mat3d::identity();
    coupled(0, 1) = coupled(1, 0) = 0.5;
    const auto limited = symmetric_eigen(coupled, 1e-12, 0);
    CHECK(limited && !limited->converged);
    auto invalid = mat3d::identity();
    invalid(0, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK(!factor_qr(invalid));
    CHECK(!cholesky(invalid));
    CHECK(!symmetric_eigen(invalid));
    for (double scale : {1e-200, 1e200}) {
        const auto positive = mat3d{2} * scale;
        const auto l = cholesky(positive);
        CHECK(l);
        NEAR((*l) * transpose(*l) / scale, mat3d{2}, 1e-12);
        const auto e = symmetric_eigen(positive);
        CHECK(e && e->converged);
        NEAR(e->values / scale, vec3d{2}, 1e-12);
    }
}
CH_TEST(symmetric_eigen_extreme_off_diagonal_ratios) {
    // A tiny off-diagonal against a large diagonal difference drives the Jacobi
    // tangent far outside the ordinary range; the hypot-based rotation must
    // still converge and produce an orthonormal eigenbasis.
    for (double off : {1e-14, 1e-10, 1e-3, 1e3}) {
        auto a = mat3d{1};
        a(0, 1) = a(1, 0) = off;
        a(1, 1) = a(2, 2) = 2.0;
        const auto e = symmetric_eigen(a, 1e-14, 64);
        CHECK(e && e->converged);
        for (std::size_t i = 0; i < 3; ++i)
            NEAR(a * e->vectors.column(i), e->vectors.column(i) * e->values[i], 1e-10);
        NEAR(transpose(e->vectors) * e->vectors, mat3d::identity(), 1e-12);
        NEAR(dot(e->values, vec3d{1, 1, 1}), trace(a), 1e-12);
    }
}
CH_TEST(decomposition_is_symmetric_within_tolerance_only) {
    // The element-wise symmetry test must accept a mirrored pair that differs
    // by less than the tolerance and reject one that differs by more.
    auto near_symmetric = mat3d{2};
    near_symmetric(0, 1) = 1e-9;
    near_symmetric(1, 0) = -1e-9;
    CHECK(cholesky(near_symmetric, 1e-6));
    CHECK(symmetric_eigen(near_symmetric, 1e-6));
    auto asymmetric = mat3d{2};
    asymmetric(0, 1) = 1e-3;
    asymmetric(1, 0) = -1e-3;
    CHECK(!cholesky(asymmetric, 1e-6));
    CHECK(!symmetric_eigen(asymmetric, 1e-6));
}
