#include "boundary_helpers.hpp"
#include "test.hpp"
#include <stdexcept>

using namespace chm;
namespace curve_boundary {
constexpr bezier<double, 2, 3> cubic{{vec2d{-1, 2}, vec2d{3, 5}, vec2d{4, -2}, vec2d{7, 1}}};
constexpr std::array control{vec2d{1, 0}, vec2d{1, 1}, vec2d{0, 1}};
constexpr std::array<double, 6> knots{0, 0, 0, 1, 1, 1};
} // namespace curve_boundary

CH_TEST(curve_split_at_zero) {
    const auto [left, right] = curve_boundary::cubic.split(0.0);
    for (double t : {0.0, 0.25, 0.9, 1.0}) {
        CHECK(left.evaluate(t) == curve_boundary::cubic.control[0]);
        CHECK(right.evaluate(t) == curve_boundary::cubic.evaluate(t));
    }
}
CH_TEST(curve_split_at_one) {
    const auto [left, right] = curve_boundary::cubic.split(1.0);
    for (double t : {0.0, 0.25, 0.9, 1.0}) {
        CHECK(left.evaluate(t) == curve_boundary::cubic.evaluate(t));
        CHECK(right.evaluate(t) == curve_boundary::cubic.control[3]);
    }
}
CH_TEST(curve_constant_split_extrapolation) {
    const bezier<double, 2, 0> c{{vec2d{2, 9}}};
    const auto [a, b] = c.split(-3.0);
    CHECK(a.evaluate(2.0) == c.control[0]);
    CHECK(b.evaluate(-1.0) == c.control[0]);
}
CH_TEST(curve_float_quadratic_analytic) {
    const bezier<float, 2, 2> c{{vec2f{0, 0}, vec2f{1, 0}, vec2f{2, 1}}};
    for (int i = 0; i <= 128; ++i) {
        const float t = float(i) / 128;
        NEAR(c.evaluate(t), vec2f{2 * t, t * t}, 1e-6f);
    }
}
CH_TEST(curve_degree_eight_bernstein_reference) {
    bezier<double, 3, 8> c;
    for (auto &p : c.control)
        p = {test::sample<double>(), test::sample<double>(), test::sample<double>()};
    constexpr std::array<double, 9> choose{1, 8, 28, 56, 70, 56, 28, 8, 1};
    for (int j = 1; j < 100; ++j) {
        const double t = double(j) / 100;
        vec3d reference;
        for (std::size_t i = 0; i <= 8; ++i)
            reference += c.control[i] *
                         (choose[i] * std::pow(t, double(i)) * std::pow(1 - t, double(8 - i)));
        NEAR(c.evaluate(t), reference, 1e-12);
    }
}
CH_TEST(curve_second_derivative_reference) {
    const auto d = curve_boundary::cubic.derivative();
    const auto dd = d.derivative();
    for (double t : {0.0, 0.2, 0.8, 1.0})
        NEAR(dd.evaluate(t), (d.evaluate(t + 1e-5) - d.evaluate(t - 1e-5)) / 2e-5, 1e-8);
}
CH_TEST(curve_patch_corners) {
    bezier_patch<double, 3, 2, 3> p;
    for (auto &row : p.control)
        for (auto &v : row)
            v = {test::sample<double>(), test::sample<double>(), test::sample<double>()};
    CHECK(p.evaluate(0, 0) == p.control[0][0]);
    CHECK(p.evaluate(1, 0) == p.control[0][2]);
    CHECK(p.evaluate(0, 1) == p.control[3][0]);
    CHECK(p.evaluate(1, 1) == p.control[3][2]);
}
CH_TEST(curve_patch_mixed_partials_commute) {
    bezier_patch<double, 3, 2, 2> p;
    for (auto &row : p.control)
        for (auto &v : row)
            v = {test::sample<double>(), test::sample<double>(), test::sample<double>()};
    constexpr double u = 0.3, v = 0.7, h = 1e-5;
    NEAR((p.derivative_u(u, v + h) - p.derivative_u(u, v - h)) / (2 * h),
         (p.derivative_v(u + h, v) - p.derivative_v(u - h, v)) / (2 * h), 1e-8);
}
CH_TEST(curve_bspline_maximum_degree) {
    std::array<vec2d, 33> control{};
    std::array<double, 66> knots{};
    for (std::size_t i = 0; i < 33; ++i) {
        control[i] = {double(i), 2 * double(i)};
        knots[i + 33] = 1;
    }
    const auto c = bspline_view<double, 2, 32>::create(control, knots);
    CHECK(c);
    for (double t : {0.0, 0.01, 0.5, 0.99, 1.0}) {
        const auto p = c->evaluate(t), d = c->derivative(t);
        CHECK(p && d);
        NEAR(*p, vec2d{32 * t, 64 * t}, 1e-11);
        NEAR(*d, vec2d{32, 64}, 1e-11);
    }
}
CH_TEST(curve_bspline_insufficient_controls) {
    const std::array<vec2d, 2> c{};
    const std::array<double, 5> k{0, 0, 0, 1, 1};
    CHECK(!bspline_view<double, 2, 2>::create(c, k));
}
CH_TEST(curve_bspline_nan_knot) {
    auto k = curve_boundary::knots;
    k[4] = boundary::nan;
    CHECK(!bspline_view<double, 2, 2>::create(curve_boundary::control, k));
}
CH_TEST(curve_bspline_infinite_knot) {
    auto k = curve_boundary::knots;
    k[5] = boundary::inf;
    CHECK(!bspline_view<double, 2, 2>::create(curve_boundary::control, k));
}
CH_TEST(curve_bspline_zero_domain) {
    const std::array<double, 6> k{};
    CHECK(!bspline_view<double, 2, 2>::create(curve_boundary::control, k));
}
CH_TEST(curve_bspline_nonfinite_control) {
    auto c = curve_boundary::control;
    c[1][1] = boundary::nan;
    CHECK(!bspline_view<double, 2, 2>::create(c, curve_boundary::knots));
}
CH_TEST(curve_bspline_discontinuity_right_hand_value) {
    const std::array c{vec2d{0, 0}, vec2d{1, 1}, vec2d{5, 5}, vec2d{6, 6}};
    const std::array<double, 6> k{0, 0, 0.5, 0.5, 1, 1};
    const auto s = bspline_view<double, 2, 1>::create(c, k);
    CHECK(s);
    CHECK(s->evaluate(0.5) == c[2]);
    NEAR(*s->evaluate(std::nextafter(0.5, 0.0)), c[1], 1e-14);
    CHECK(s->derivative(0.5) == vec2d{2, 2});
}
CH_TEST(curve_bspline_derivative_overflow_cancellation) {
    const std::array c{vec2d{-1e308, 0}, vec2d{1e308, 0}};
    const std::array<double, 4> k{-1e308, -1e308, 1e308, 1e308};
    const auto s = bspline_view<double, 2, 1>::create(c, k);
    CHECK(s);
    const auto d = s->derivative(0.0);
    CHECK(d);
    NEAR(*d, vec2d{1, 0}, 1e-14);
    const std::array<double, 4> finite_width{0, 0, 1e308, 1e308};
    const auto finite = bspline_view<double, 2, 1>::create(c, finite_width);
    CHECK(finite);
    const auto twice = finite->derivative(5e307);
    CHECK(twice);
    NEAR(*twice, vec2d{2, 0}, 1e-14);
}
CH_TEST(curve_bspline_derivative_subnormal_domain) {
    const std::array c{vec2d{0, 0}, vec2d{1e-310, 2e-310}};
    const std::array<double, 4> k{0, 0, 1e-310, 1e-310};
    const auto s = bspline_view<double, 2, 1>::create(c, k);
    CHECK(s);
    const auto d = s->derivative(5e-311);
    CHECK(d);
    NEAR(*d, vec2d{1, 2}, 1e-12);
}
CH_TEST(curve_nurbs_nan_weight) {
    const std::array w{1.0, boundary::nan, 1.0};
    CHECK(!nurbs_view<double, 2, 2>::create(curve_boundary::control, w, curve_boundary::knots));
}
CH_TEST(curve_nurbs_infinite_weight) {
    const std::array w{1.0, boundary::inf, 1.0};
    CHECK(!nurbs_view<double, 2, 2>::create(curve_boundary::control, w, curve_boundary::knots));
}
CH_TEST(curve_nurbs_subnormal_weight_scale_invariance) {
    const std::array w{1e-310, 2e-310, 1e-310}, reference_weights{1.0, 2.0, 1.0};
    const auto c =
        nurbs_view<double, 2, 2>::create(curve_boundary::control, w, curve_boundary::knots);
    const auto r = nurbs_view<double, 2, 2>::create(curve_boundary::control, reference_weights,
                                                    curve_boundary::knots);
    CHECK(c && r);
    for (double t : {0.0, 0.1, 0.5, 0.9, 1.0}) {
        const auto p = c->evaluate(t), q = r->evaluate(t);
        CHECK(p && q);
        NEAR(*p, *q, 1e-12);
    }
}
CH_TEST(numeric_quadratic_zero_constant_distinct_roots) {
    const auto r = solve_quadratic(1.0, -3.0, 0.0);
    CHECK(r && r->count == 2);
    CHECK(r->values[0] == 0 && r->values[1] == 3);
}
CH_TEST(numeric_quadratic_negative_leading_sorted_roots) {
    const auto r = solve_quadratic(-1.0, 1.0, 6.0);
    CHECK(r && r->count == 2);
    NEAR(r->values[0], -2.0, 1e-14);
    NEAR(r->values[1], 3.0, 1e-14);
}
CH_TEST(numeric_quadratic_nan_linear_coefficient) {
    CHECK(!solve_quadratic(1.0, boundary::nan, 1.0));
}
CH_TEST(numeric_quadratic_subnormal_coefficients) {
    const auto r = solve_quadratic(boundary::tiny, -3 * boundary::tiny, 2 * boundary::tiny);
    CHECK(r && r->count == 2);
    NEAR(r->values[0], 1.0, 1e-14);
    NEAR(r->values[1], 2.0, 1e-14);
}
CH_TEST(numeric_linear_root_unrepresentable) {
    CHECK(!solve_quadratic(0.0, boundary::tiny, 1.0));
}
CH_TEST(numeric_bisection_zero_iteration_budget) {
    std::size_t calls = 0;
    const auto r = bisect(
        [&](double x) {
            ++calls;
            return x * x - 2;
        },
        0.0, 2.0, 0.0, 0.0, 0);
    CHECK(r && !r->converged && r->iterations == 0);
    CHECK(calls == 2);
    CHECK(r->residual == -2);
}
CH_TEST(numeric_bisection_single_point_root) {
    const auto r = bisect([](double x) { return x - 2; }, 2.0, 2.0);
    CHECK(r && r->converged && r->value == 2 && r->iterations == 0);
}
CH_TEST(numeric_bisection_extreme_midpoint) {
    const auto r =
        bisect([](double x) { return x; }, -boundary::maximum, boundary::maximum, 0.0, 0.0);
    CHECK(r && r->converged && r->value == 0 && r->iterations == 1);
}
CH_TEST(numeric_bisection_nonfinite_interior) {
    CHECK(!bisect([](double x) { return x == 0.5 ? boundary::nan : x - 0.5; }, 0.0, 1.0));
}
CH_TEST(numeric_bisection_callback_exception_propagates) {
    bool caught = false;
    try {
        (void)bisect([](double) -> double { throw std::runtime_error("callback"); }, 0.0, 1.0);
    } catch (const std::runtime_error &) {
        caught = true;
    }
    CHECK(caught);
}
CH_TEST(numeric_integration_minimum_budget_three) {
    std::size_t calls = 0;
    const auto r = integrate(
        [&](double x) {
            ++calls;
            return std::exp(x);
        },
        0.0, 1.0, 1e-12, 20, 3);
    CHECK(r && !r->converged && r->evaluations == 3 && calls == 3);
    CHECK(std::isinf(r->estimated_error));
}
CH_TEST(numeric_integration_odd_remaining_budget) {
    std::size_t calls = 0;
    const auto r = integrate(
        [&](double x) {
            ++calls;
            return std::exp(x);
        },
        0.0, 1.0, 1e-12, 20, 4);
    CHECK(r && !r->converged && r->evaluations == 3 && calls == 3);
}
CH_TEST(numeric_integration_zero_depth_budget) {
    std::size_t calls = 0;
    const auto r = integrate(
        [&](double x) {
            ++calls;
            return x * x;
        },
        0.0, 1.0, 1e-12, 0, 100);
    CHECK(r && !r->converged && calls == 3);
    NEAR(r->value, 1.0 / 3, 1e-14);
}
CH_TEST(numeric_integration_zero_interval_does_not_call) {
    const auto r =
        integrate([](double) -> double { throw std::runtime_error("must not call"); }, 2.0, 2.0);
    CHECK(r && r->converged && r->evaluations == 0 && r->value == 0);
}
CH_TEST(numeric_integration_nan_tolerance) {
    CHECK(!integrate([](double x) { return x; }, 0.0, 1.0, boundary::nan));
}
CH_TEST(numeric_integration_nonfinite_interval) {
    CHECK(!integrate([](double x) { return x; }, 0.0, boundary::inf));
}
CH_TEST(numeric_integration_large_samples_small_interval) {
    const auto r = integrate([](double) { return 1e308; }, 0.0, 1e-308, 1e-12);
    CHECK(r && r->converged);
    NEAR(r->value, 1.0, 1e-12);
}
CH_TEST(numeric_integration_overflow_width_finite_value) {
    const auto r = integrate([](double) { return 1e-308; }, -1e308, 1e308, 1e-12);
    CHECK(r && r->converged);
    NEAR(r->value, 2.0, 1e-12);
}
CH_TEST(numeric_compensated_negative_cancellation) {
    compensated_sum<double> sum;
    for (int i = 0; i < 1000; ++i) {
        sum.add(-1e20);
        sum.add(-0.25);
        sum.add(1e20);
    }
    CHECK(sum.value() == -250);
    sum.reset();
    sum.add(boundary::tiny);
    CHECK(sum.value() == boundary::tiny);
}
CH_TEST(numeric_polynomial_high_degree_zero_argument) {
    std::array<double, 65> coefficients{};
    coefficients[0] = 7;
    coefficients[64] = 1e200;
    CHECK(evaluate_polynomial<double>(coefficients, 0.0) == 7);
    NEAR(evaluate_polynomial<double>(coefficients, 0.5), 7 + std::ldexp(1e200, -64), 1e-14);
}
