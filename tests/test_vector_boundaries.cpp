#include "boundary_helpers.hpp"
#include "test.hpp"
using namespace chm;

CH_TEST(vector_normalize_one_dimension_negative) {
    CHECK(normalize(vec<double, 1>{-7.0}) == vec<double, 1>{-1.0});
}
CH_TEST(vector_normalize_two_dimensions) {
    NEAR(normalize(vec2d{5, 12}), vec2d{5.0 / 13, 12.0 / 13}, 1e-14);
}
CH_TEST(vector_normalize_four_dimensions) {
    NEAR(normalize(vec4d{1, 1, 1, 1}), vec4d{0.5}, 1e-14);
}
CH_TEST(vector_normalize_seven_dimensions) {
    NEAR(length(normalize(boundary::sequence<double, 7>())), 1.0, 1e-14);
}
CH_TEST(vector_normalize_maximum_float) {
    const auto v = normalize(vec3f{std::numeric_limits<float>::max()});
    NEAR(length(v), 1.0f, 1e-6f);
}
CH_TEST(vector_normalize_minimum_double) {
    const auto v = normalize(vec3d{std::numeric_limits<double>::min()});
    NEAR(length(v), 1.0, 1e-14);
}
CH_TEST(vector_normalize_denormal_axis) {
    CHECK(normalize(vec3d{0, boundary::tiny, 0}) == vec3d{0, 1, 0});
}
CH_TEST(vector_normalize_all_denormals) {
    NEAR(length(normalize(vec3d{boundary::tiny})), 1.0, 1e-14);
}
CH_TEST(vector_normalize_mixed_magnitude) {
    NEAR(normalize(vec3d{1e300, 1e-300, -1e300}), vec3d{std::sqrt(0.5), 0, -std::sqrt(0.5)}, 1e-14);
}
CH_TEST(vector_normalize_negative_zero) {
    CHECK(!try_normalize(vec3d{-0.0, 0.0, -0.0}));
}
CH_TEST(vector_nan_dominates_infinity_norm) {
    CHECK(std::isnan(length(vec3d{boundary::inf, boundary::nan, 1})));
    CHECK(std::isnan(length(vec3d{1, boundary::nan, boundary::inf})));
}
CH_TEST(vector_infinite_component_norm) {
    CHECK(std::isinf(length(vec3d{1, -boundary::inf, 3})));
}
CH_TEST(vector_norm_representable_max_axis) {
    CHECK(length(vec3d{boundary::maximum, 0, 0}) == boundary::maximum);
}
CH_TEST(vector_distance_underflow_resistance) {
    NEAR(distance(vec3d{1e-300, 0, 0}, vec3d{}), 1e-300, 1e-313);
}
CH_TEST(vector_angle_nearly_parallel) {
    const auto a = angle_between(vec3d{1, 0, 0}, vec3d{1, 1e-12, 0});
    CHECK(a);
    NEAR(*a, 1e-12, 1e-24);
}
CH_TEST(vector_angle_nearly_antiparallel) {
    const auto a = angle_between(vec3d{1, 0, 0}, vec3d{-1, 1e-8, 0});
    CHECK(a);
    NEAR(*a, pi<double> - 1e-8, 1e-14);
}
CH_TEST(vector_angle_nonfinite_rejected) {
    CHECK(!angle_between(vec3d{boundary::inf, 0, 0}, vec3d{1, 0, 0}));
}
CH_TEST(vector_projection_rejects_nonfinite_source) {
    CHECK(!project(vec3d{boundary::nan, 0, 0}, vec3d{1, 0, 0}));
}
CH_TEST(vector_projection_extreme_axis) {
    const auto p = project(vec3d{3, 4, 5}, vec3d{1e300, 0, 0});
    CHECK(p);
    CHECK(*p == vec3d{3, 0, 0});
}
CH_TEST(vector_projection_orthogonal_zero) {
    const auto p = project(vec3d{0, 2, 0}, vec3d{1, 0, 0});
    CHECK(p);
    CHECK(*p == vec3d{});
}
CH_TEST(vector_reflection_involution) {
    const auto n = normalize(vec3d{1, 2, 3});
    const vec3d v{4, -1, 2};
    NEAR(reflect(reflect(v, n), n), v, 1e-13);
}
CH_TEST(vector_reflection_preserves_length) {
    const vec3d v{4, -1, 2};
    NEAR(length(reflect(v, normalize(vec3d{1, 2, 3}))), length(v), 1e-13);
}
CH_TEST(vector_refraction_nan_ratio_rejected) {
    CHECK(!refract(vec3d{0, -1, 0}, vec3d{0, 1, 0}, boundary::nan));
}
CH_TEST(vector_refraction_nonpositive_ratio_rejected) {
    CHECK(!refract(vec3d{0, -1, 0}, vec3d{0, 1, 0}, 0.0));
    CHECK(!refract(vec3d{0, -1, 0}, vec3d{0, 1, 0}, -1.0));
}
CH_TEST(vector_refraction_critical_angle) {
    const auto r = refract(vec3d{0.5, -std::sqrt(0.75), 0}, vec3d{0, 1, 0}, 2.0);
    CHECK(r);
    NEAR(*r, vec3d{1, 0, 0}, 1e-7);
}
CH_TEST(vector_refraction_normal_extreme_ratio) {
    const auto r = refract(vec3d{0, -1, 0}, vec3d{0, 1, 0}, 1e300);
    CHECK(r);
    CHECK(*r == vec3d{0, -1, 0});
}
CH_TEST(vector_cross_lagrange_identity) {
    const vec3d a{2, -3, 4}, b{5, 6, -7};
    NEAR(length_squared(cross(a, b)), length_squared(a) * length_squared(b) - square(dot(a, b)),
         1e-12);
}
CH_TEST(vector_cross_orientation_2d) {
    CHECK(cross(vec2d{1, 2}, vec2d{3, 4}) == -2.0);
    CHECK(cross(vec2d{3, 4}, vec2d{1, 2}) == 2.0);
}
CH_TEST(vector_unsigned_values_and_hadamard) {
    const vec<unsigned, 3> v{0U, 1U, 65535U};
    CHECK(hadamard(v, vec<unsigned, 3>{2U}) == vec<unsigned, 3>{0U, 2U, 131070U});
}
CH_TEST(vector_explicit_narrowing_representable) {
    const vec3d d{1.5, -2.5, 3.25};
    const vec3f f(d);
    CHECK(vec3d(f) == d);
    CHECK(vec3i(d) == vec3i{1, -2, 3});
}
