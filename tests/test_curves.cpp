#include "test.hpp"
#include <chmath/curves/curves.hpp>

using namespace chm;
constexpr bezier<double, 2, 1> constexpr_line{{vec2d{1, 2}, vec2d{3, 4}}};
static_assert(constexpr_line.evaluate(0.5) == vec2d{2, 3});

CH_TEST(bezier_evaluation_split_and_derivative) {
    const bezier<double, 3, 3> curve{
        {vec3d{0, 0, 0}, vec3d{1, 3, 0}, vec3d{3, 2, 1}, vec3d{4, 0, 2}}};
    CHECK(curve.evaluate(0) == curve.control[0]);
    CHECK(curve.evaluate(1) == curve.control[3]);
    const auto d = curve.derivative();
    CHECK(d.evaluate(0) == (curve.control[1] - curve.control[0]) * 3.0);
    for (int i = 0; i <= 100; ++i) {
        const double t = double(i) / 100;
        const auto v = curve.evaluate(t);
        const double s = 1 - t;
        const auto reference = curve.control[0] * (s * s * s) + curve.control[1] * (3 * s * s * t) +
                               curve.control[2] * (3 * s * t * t) + curve.control[3] * (t * t * t);
        NEAR(v, reference, 1e-12);
        const double h = 1e-5;
        NEAR(d.evaluate(t), (curve.evaluate(t + h) - curve.evaluate(t - h)) / (2 * h), 1e-8);
        const auto [left, right] = curve.split(0.3);
        NEAR(left.evaluate(t), curve.evaluate(0.3 * t), 1e-12);
        NEAR(right.evaluate(t), curve.evaluate(0.3 + 0.7 * t), 1e-12);
    }
    const bezier<double, 3, 0> constant{{vec3d{1, 2, 3}}};
    CHECK(constant.evaluate(0.9) == vec3d{1, 2, 3});
    CHECK(hermite(vec2d{1, 2}, vec2d{3, 4}, vec2d{5, 6}, vec2d{7, 8}, 0.0) == vec2d{1, 2});
    CHECK(hermite(vec2d{1, 2}, vec2d{3, 4}, vec2d{5, 6}, vec2d{7, 8}, 1.0) == vec2d{5, 6});
    CHECK(catmull_rom(vec2d{0, 0}, vec2d{1, 1}, vec2d{2, 0}, vec2d{3, 1}, 0.0) == vec2d{1, 1});
    CHECK(catmull_rom(vec2d{0, 0}, vec2d{1, 1}, vec2d{2, 0}, vec2d{3, 1}, 1.0) == vec2d{2, 0});
}
CH_TEST(bezier_patch_and_partials) {
    bezier_patch<double, 3, 1, 1> patch;
    patch.control = {std::array{vec3d{0, 0, 0}, vec3d{2, 0, 0}},
                     std::array{vec3d{0, 3, 0}, vec3d{2, 3, 2}}};
    NEAR(patch.evaluate(0.5, 0.5), vec3d{1, 1.5, 0.5}, 1e-14);
    NEAR(patch.derivative_u(0.5, 0.5), vec3d{2, 0, 1}, 1e-14);
    NEAR(patch.derivative_v(0.5, 0.5), vec3d{0, 3, 1}, 1e-14);
}
CH_TEST(bspline_matches_bezier_and_derivative) {
    const std::array control{vec3d{0, 0, 0}, vec3d{1, 3, 0}, vec3d{3, 2, 1}, vec3d{4, 0, 2}};
    const std::array<double, 8> knots{0, 0, 0, 0, 1, 1, 1, 1};
    const auto spline = bspline_view<double, 3, 3>::create(control, knots);
    CHECK(spline);
    CHECK(spline->domain_min() == 0 && spline->domain_max() == 1);
    const bezier<double, 3, 3> b{control};
    for (int i = 0; i <= 100; ++i) {
        const double t = double(i) / 100;
        CHECK(spline->evaluate(t));
        NEAR(*spline->evaluate(t), b.evaluate(t), 1e-12);
        CHECK(spline->derivative(t));
        NEAR(*spline->derivative(t), b.derivative().evaluate(t), 1e-12);
    }
    CHECK(!spline->evaluate(-0.1));
    CHECK(!spline->evaluate(1.1));
    CHECK(!spline->evaluate(std::numeric_limits<double>::quiet_NaN()));
    CHECK(!spline->derivative(2.0));
    const std::array<double, 8> bad{0, 0, 0, 0, 0, 1, 1, 1};
    CHECK(!bspline_view<double, 3, 3>::create(control, bad));
    auto unordered = knots;
    unordered[5] = -1;
    CHECK(!bspline_view<double, 3, 3>::create(control, unordered));
    CHECK(!bspline_view<double, 3, 3>::create(control, std::span<const double>{knots}.first(7)));
}
CH_TEST(bspline_repeated_knots_and_nonclamped) {
    const std::array control{vec2d{0, 0}, vec2d{1, 1}, vec2d{2, 0}, vec2d{3, 1}, vec2d{4, 0}};
    const std::array<double, 8> knots{0, 0, 0, 0.5, 0.5, 1, 1, 1};
    const auto spline = bspline_view<double, 2, 2>::create(control, knots);
    CHECK(spline);
    NEAR(*spline->evaluate(0.5), control[2], 1e-12);
    for (double t : {0.1, 0.3, 0.7, 0.9})
        NEAR(*spline->derivative(t),
             (*spline->evaluate(t + 1e-6) - *spline->evaluate(t - 1e-6)) / (2e-6), 1e-8);
    const std::array<double, 8> uniform{0, 1, 2, 3, 4, 5, 6, 7};
    const auto open = bspline_view<double, 2, 2>::create(control, uniform);
    CHECK(open);
    NEAR(*open->evaluate(2.0), (control[0] + control[1]) / 2.0, 1e-12);
    NEAR(*open->evaluate(5.0), (control[3] + control[4]) / 2.0, 1e-12);
    const std::array constant{vec2d{1, 2}, vec2d{3, 4}};
    const std::array<double, 3> steps{0, 1, 2};
    const auto degree0 = bspline_view<double, 2, 0>::create(constant, steps);
    CHECK(degree0);
    CHECK(degree0->evaluate(0.5) == constant[0]);
    CHECK(degree0->evaluate(1.0) == constant[1]);
    CHECK(degree0->evaluate(2.0) == constant[1]);
    const std::array linear{vec2d{0, 0}, vec2d{2, 4}};
    const std::array<double, 4> linear_knots{0, 0, 1, 1};
    const auto line = bspline_view<double, 2, 1>::create(linear, linear_knots);
    CHECK(line);
    CHECK(line->derivative(0.0) == vec2d{2, 4});
    CHECK(line->derivative(1.0) == vec2d{2, 4});
}
CH_TEST(nurbs_exact_quarter_circle) {
    const std::array control{vec2d{1, 0}, vec2d{1, 1}, vec2d{0, 1}};
    const std::array<double, 3> weights{1, std::sqrt(0.5), 1};
    const std::array<double, 6> knots{0, 0, 0, 1, 1, 1};
    const auto curve = nurbs_view<double, 2, 2>::create(control, weights, knots);
    CHECK(curve);
    CHECK(curve->domain_min() == 0 && curve->domain_max() == 1);
    CHECK(curve->evaluate(0) == control[0]);
    CHECK(curve->evaluate(1) == control[2]);
    for (int i = 0; i <= 100; ++i) {
        const auto p = curve->evaluate(double(i) / 100);
        CHECK(p);
        NEAR(length(*p), 1.0, 1e-12);
    }
    NEAR(*curve->evaluate(0.5), vec2d{std::sqrt(0.5), std::sqrt(0.5)}, 1e-12);
    const std::array<double, 3> huge{1e300, std::sqrt(0.5) * 1e300, 1e300};
    const auto same = nurbs_view<double, 2, 2>::create(control, huge, knots);
    CHECK(same);
    NEAR(*same->evaluate(0.3), *curve->evaluate(0.3), 1e-12);
    auto bad = weights;
    bad[0] = 0;
    CHECK(!nurbs_view<double, 2, 2>::create(control, bad, knots));
    bad[0] = -1;
    CHECK(!nurbs_view<double, 2, 2>::create(control, bad, knots));
    CHECK(!nurbs_view<double, 2, 2>::create(control, std::span<const double>{weights}.first(2),
                                            knots));
    CHECK(!curve->evaluate(-1));
    const std::array<double, 3> ones{1, 1, 1};
    const auto polynomial = nurbs_view<double, 2, 2>::create(control, ones, knots);
    CHECK(polynomial);
    NEAR(*polynomial->evaluate(0.7), (bezier<double, 2, 2>{control}.evaluate(0.7)), 1e-12);
}
