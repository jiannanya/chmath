#include "boundary_helpers.hpp"
#include "test.hpp"
using namespace chm;

CH_TEST(scalar_clamp_lowest_integer) {
    CHECK(clamp(std::numeric_limits<int>::lowest(), -7, 9) == -7);
}
CH_TEST(scalar_clamp_highest_integer) {
    CHECK(clamp(std::numeric_limits<int>::max(), -7, 9) == 9);
}
CH_TEST(scalar_clamp_equal_endpoints) {
    CHECK(clamp(-100.0, 3.0, 3.0) == 3.0);
    CHECK(clamp(100.0, 3.0, 3.0) == 3.0);
}
CH_TEST(scalar_saturate_infinities) {
    CHECK(saturate(boundary::inf) == 1);
    CHECK(saturate(-boundary::inf) == 0);
}
CH_TEST(scalar_clamp_nan_propagation) {
    CHECK(std::isnan(clamp(boundary::nan, 0.0, 1.0)));
}
CH_TEST(scalar_smoothstep_huge_interval) {
    NEAR(smoothstep(-boundary::maximum, boundary::maximum, 0.0), 0.5, 1e-15);
}
CH_TEST(scalar_smoothstep_zero_width) {
    CHECK(smoothstep(4.0, 4.0, 3.0) == 0);
    CHECK(smoothstep(4.0, 4.0, 4.0) == 1);
}
CH_TEST(scalar_smoothstep_clipped_domain) {
    CHECK(smoothstep(-2.0, 2.0, -3.0) == 0);
    CHECK(smoothstep(-2.0, 2.0, 3.0) == 1);
}
CH_TEST(scalar_smoothstep_monotonicity) {
    double prior = 0;
    for (int i = 0; i <= 100; ++i) {
        const auto v = smoothstep(0.0, 1.0, double(i) / 100);
        CHECK(v >= prior);
        prior = v;
    }
}
CH_TEST(scalar_lerp_extreme_opposite_signs) {
    CHECK(lerp(-boundary::maximum, boundary::maximum, 0.5) == 0);
}
CH_TEST(scalar_lerp_same_sign_large) {
    NEAR(lerp(1e308, 1.2e308, 0.5), 1.1e308, 1e-15);
}
CH_TEST(scalar_lerp_extrapolation) {
    CHECK(lerp(2.0, 4.0, -1.0) == 0);
    CHECK(lerp(2.0, 4.0, 2.0) == 6);
}
CH_TEST(scalar_lerp_exact_endpoints) {
    CHECK(lerp(-17.25, 91.5, 0.0) == -17.25);
    CHECK(lerp(-17.25, 91.5, 1.0) == 91.5);
}
CH_TEST(scalar_relative_tolerance_overflow) {
    CHECK(!almost_equal(boundary::maximum, -boundary::maximum, 1.5, 0.0));
    CHECK(almost_equal(boundary::maximum, -boundary::maximum, 2.0, 0.0));
}
CH_TEST(scalar_rejects_negative_tolerance) {
    CHECK(!almost_equal(1.0, 1.0, -1.0, 0.0));
    CHECK(!almost_equal(1.0, 1.0, 0.0, -1.0));
}
CH_TEST(scalar_rejects_nan_tolerance) {
    CHECK(!almost_equal(1.0, 1.0, boundary::nan, 0.0));
    CHECK(!almost_equal(1.0, 1.0, 0.0, boundary::nan));
}
CH_TEST(scalar_rejects_infinite_tolerance) {
    CHECK(!almost_equal(1.0, 2.0, boundary::inf, 0.0));
    CHECK(!almost_equal(1.0, 2.0, 0.0, boundary::inf));
}
CH_TEST(scalar_comparison_subnormal_gap) {
    CHECK(!almost_equal(0.0, boundary::tiny, 0.0, 0.0));
    CHECK(almost_equal(0.0, boundary::tiny, 0.0, boundary::tiny));
}
CH_TEST(scalar_signed_zero_equality) {
    CHECK(almost_equal(-0.0, 0.0, 0.0, 0.0));
    CHECK(is_finite(-0.0));
}
CH_TEST(scalar_wrap_negative_pi) {
    CHECK(wrap_angle(-pi<double>) == -pi<double>);
}
CH_TEST(scalar_wrap_many_revolutions) {
    NEAR(wrap_angle(64 * tau<double> + 0.25), 0.25, 1e-12);
}
CH_TEST(scalar_wrap_nonfinite) {
    CHECK(std::isnan(wrap_angle(boundary::inf)));
    CHECK(std::isnan(wrap_angle(boundary::nan)));
}
CH_TEST(scalar_angle_unit_round_trip) {
    for (double a : {-1e6, -180.0, 0.0, 90.0, 1e6})
        NEAR(degrees(radians(a)), a, 1e-14);
}
CH_TEST(scalar_constants_float_precision) {
    CHECK(epsilon<float> == 32 * std::numeric_limits<float>::epsilon());
    NEAR(tau<float>, 2 * pi<float>, 1e-6f);
}
CH_TEST(scalar_multiply_divide_overflow_cancel) {
    NEAR(multiply_divide(1e300, 1e200, 1e300), 1e200, 1e-14);
}
CH_TEST(scalar_multiply_divide_underflow_cancel) {
    NEAR(multiply_divide(1e-300, 1e-200, 1e-300), 1e-200, 1e-213);
}
CH_TEST(scalar_multiply_divide_negative_sign) {
    NEAR(multiply_divide(-1e300, 1e200, -1e300), 1e200, 1e-14);
}
CH_TEST(scalar_multiply_divide_subnormal_exact) {
    CHECK(multiply_divide(boundary::tiny, 2.0, 2.0) == boundary::tiny);
}
CH_TEST(scalar_multiply_divide_ieee_specials) {
    CHECK(std::isnan(multiply_divide(0.0, boundary::inf, 1.0)));
    CHECK(std::isinf(multiply_divide(1.0, 1.0, 0.0)));
}
CH_TEST(scalar_unsigned_modular_arithmetic) {
    using U = unsigned;
    CHECK(square(U(65535)) == U(4294836225U));
    CHECK(clamp(U(0), U(1), U(5)) == 1U);
}
