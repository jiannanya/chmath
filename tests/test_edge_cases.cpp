#include "test.hpp"
#include <chmath/chmath.hpp>

using namespace chm;

CH_TEST(extreme_finite_geometry_regressions) {
    NEAR(multiply_divide(2.0, 1e308, 1e308), 2.0, 1e-15);
    NEAR(multiply_divide(1e-200, 1e-200, 1e-300), 1e-100, 1e-114);
    CHECK(multiply_divide(0.0, 1e308, 1e308) == 0.0);
    CHECK(multiply_divide(2.0, 3.0, 4.0) == 1.5);
    const ray3d ray{{1e308, 0, 0}, {-1e308, 0, 0}};
    const auto box = intersect(ray, aabb3d{{-1e308, -1, -1}, {-9e307, 1, 1}});
    CHECK(box);
    NEAR(box->enter, 1.9, 1e-14);
    NEAR(box->exit, 2.0, 1e-14);
    const auto ball = intersect(ray, sphered{{}, 1e308});
    CHECK(ball);
    NEAR(ball->enter, 0.0, 1e-14);
    NEAR(ball->exit, 2.0, 1e-14);
    const auto far = intersect(ray, sphered{{-1e308, 0, 0}, 1e307});
    CHECK(far);
    NEAR(far->enter, 1.9, 1e-12);
    NEAR(far->exit, 2.1, 1e-12);
    const triangle<double> thin{{0, 0, 0}, {1, 0, 0}, {1, 1e-9, 0}};
    const vec3d query{0.75, 0.25e-9, 1};
    const auto weights = barycentric(thin, query, 0.0);
    CHECK(weights);
    NEAR(*weights, vec3d{0.25, 0.5, 0.25}, 1e-14);
    const auto projected = closest_point(thin, query);
    NEAR(projected[0], 0.75, 1e-14);
    NEAR(projected[1], 0.25e-9, 1e-23);
    CHECK(projected[2] == 0.0);
    const std::array controls{vec2d{0, 0}, vec2d{2, 4}};
    const std::array<double, 4> knots{-1e308, -1e308, 1e308, 1e308};
    const auto spline = bspline_view<double, 2, 1>::create(controls, knots);
    CHECK(spline);
    NEAR(*spline->evaluate(0.0), vec2d{1, 2}, 1e-14);
    NEAR((*spline->derivative(0.0)) * 1e308, vec2d{1, 2}, 1e-14);
}
CH_TEST(nonfinite_input_rejection) {
    const double nan = std::numeric_limits<double>::quiet_NaN(),
                 inf = std::numeric_limits<double>::infinity();
    const vec3d invalid{nan, 0, 0};
    const quatd qnan{nan, 0, 0, 1};
    CHECK(!try_normalize(qnan));
    CHECK(!inverse(qnan));
    CHECK(std::isinf(length(vec3d{inf, 0, 0})));
    CHECK(!make_plane(invalid, 1.0));
    CHECK(!make_plane(vec3d{0, 1, 0}, nan));
    CHECK(!from_matrix(mat3d::identity(), nan));
    CHECK(!factor_lu(mat3d::identity(), nan));
    CHECK(!symmetric_eigen(mat3d::identity(), -1.0));
    CHECK(!cholesky(mat3d::identity(), -1.0));
    CHECK(!least_squares(mat3d::identity(), invalid));
    affine3d bad;
    bad.matrix(0, 0) = nan;
    CHECK(!inverse(bad));
    CHECK(!decompose(bad));
    CHECK(!look_at(invalid, vec3d{}, vec3d{0, 1, 0}));
    auto m = mat4d::identity();
    m(0, 0) = nan;
    CHECK(!extract_frustum(m));
    CHECK(!to_affine(m));
    const auto fr = extract_frustum(mat4d::identity());
    CHECK(fr);
    CHECK(classify(*fr, invalid) == containment::outside);
    CHECK(classify(*fr, sphered{{}, -1}) == containment::outside);
    CHECK(classify(*fr, aabb3d{{-0.1, -0.1, -0.1}, {0.1, 0.1, 0.1}}) == containment::intersecting);
    CHECK(classify(*fr, aabb3d{{10, 10, 10}, {11, 11, 11}}) == containment::outside);
    CHECK(!project(vec3d{}, mat4d::identity(), viewport<double>{nan, 0, 1, 1}));
    CHECK(!unproject(vec3d{}, mat4d::identity(), viewport<double>{nan, 0, 1, 1}));
    const ray3d r{{0, 0, 1}, {0, 0, -1}};
    CHECK(!intersect(r, sphered{{}, -1}));
    CHECK(!intersect(r, plane<double>{{}, 0}));
    CHECK(!intersect(r, triangle<double>{invalid, {}, {}}));
    const std::array cp{invalid, vec3d{1, 2, 3}};
    const std::array<double, 4> knots{0, 0, 1, 1};
    const std::array<double, 2> w{1, 1};
    CHECK(!bspline_view<double, 3, 1>::create(cp, knots));
    CHECK(!nurbs_view<double, 3, 1>::create(cp, w, knots));
}
CH_TEST(nonrepresentable_results_and_numeric_failure_paths) {
    const double tiny = std::numeric_limits<double>::denorm_min(),
                 huge = std::numeric_limits<double>::max(),
                 nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(!inverse(mat3d{tiny}));
    CHECK(!solve(mat3d{tiny}, vec3d{1, 1, 1}));
    CHECK(!inverse(quatd{tiny, 0, 0, 0}));
    CHECK(!solve_quadratic(0.0, tiny, 1.0));
    CHECK(!solve_quadratic(tiny, huge, 1.0));
    CHECK(!solve_quadratic(tiny, 1.0, 1.0));
    CHECK(!make_plane(vec3d{tiny, 0, 0}, 1.0));
    const auto identity_rotation = rotation_between(vec3d{0, 0, 1}, vec3d{0, 0, 2});
    CHECK(identity_rotation);
    CHECK(*identity_rotation == quatd{});
    CHECK(to_axis_angle(quatd{}).angle == 0.0);
    const auto root_at_hi = bisect([](double x) { return x - 1; }, 0.0, 1.0);
    CHECK(root_at_hi && root_at_hi->value == 1.0);
    CHECK(!bisect([nan](double x) { return x == 0.5 ? nan : x - 0.4; }, 0.0, 1.0));
    CHECK(!integrate([nan](double x) { return x == 0.25 ? nan : x * x; }, 0.0, 1.0));
    const auto representable = integrate([huge](double) { return huge; }, 0.0, 1.0);
    CHECK(representable && representable->converged);
    CHECK(representable->value == huge);
    CHECK(!integrate([huge](double) { return huge; }, 0.0, 2.0));
    CHECK(!integrate([](double x) { return x; }, 0.0, 1.0, 1e-6, 10, 2));
    CHECK(!intersect(ray3d{{0, 0, 0}, {1, 0, 0}}, sphered{{}, 0}, 1.0, 2.0));
}
CH_TEST(remaining_value_operations) {
    const mat2d a{2}, b{3};
    CHECK(a + b == mat2d{5});
    CHECK(a - b == mat2d{-1});
    CHECK(-a == mat2d{-2});
    CHECK(2.0 * a == mat2d{4});
    CHECK(a / 2.0 == mat2d::identity());
    const mat2f f(a);
    CHECK(f(0, 0) == 2.0f);
    const std::array<double, 4> data{1, 2, 3, 4};
    CHECK(mat2d(data)(1, 0) == 2.0);
    NEAR(to_matrix(rotation(quatd{})), mat4d::identity(), 1e-15);
    const quatd q{1, 2, 3, 4};
    CHECK(q - q == quatd{0, 0, 0, 0});
    CHECK(2.0 * q == q * 2.0);
    aabb3d box;
    box.expand(aabb3d{{1, 2, 3}, {4, 5, 6}});
    CHECK(box.valid());
    const aabb<double, 2> planar{{0, 0}, {1, 1}};
    CHECK(intersect(ray<double, 2>{{2, 0.5}, {-1, 0}}, planar));
    const auto singleton = bspline_view<double, 2, 0>::create({}, {});
    CHECK(!singleton);
}
