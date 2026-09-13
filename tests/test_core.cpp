#include "test.hpp"
#include <chmath/matrix/matrix.hpp>
#include <type_traits>

using namespace chm;
constexpr vec3i constexpr_a{1, 2, 3}, constexpr_b{4, 5, 6};
static_assert(dot(constexpr_a, constexpr_b) == 32);
static_assert(cross(constexpr_a, constexpr_b) == vec3i{-3, 6, -3});
static_assert((mat<int, 3, 3>::identity() * constexpr_a) == constexpr_a);
static_assert(sizeof(vec3f) == 12 && sizeof(vec4f) == 16 && sizeof(mat3f) == 36 &&
              sizeof(mat4f) == 64);
static_assert(alignof(vec3f) == alignof(float));
static_assert(std::is_trivially_copyable_v<vec3f> && std::is_standard_layout_v<vec3f>);
static_assert(std::is_trivially_copyable_v<mat4d> && std::is_standard_layout_v<mat4d>);
static_assert(!std::is_convertible_v<vec3d, vec3f>);

CH_TEST(scalar_contracts) {
    CHECK(clamp(5, 0, 3) == 3);
    CHECK(saturate(-2.0) == 0.0);
    NEAR(radians(180.0), pi<double>, 1e-15);
    NEAR(degrees(pi<double>), 180.0, 1e-15);
    CHECK(lerp(-1e308, 1e308, 0.5) == 0.0);
    CHECK(smoothstep(0.0, 1.0, 0.5) == 0.5);
    CHECK(smoothstep(1.0, 1.0, 0.0) == 0.0);
    CHECK(smoothstep(1.0, 1.0, 1.0) == 1.0);
    CHECK(wrap_angle(pi<double>) == -pi<double>);
    CHECK(almost_equal(1e6, 1e6 + 1.0, 1e-5, 0.0));
    CHECK(!almost_equal(0.0, 1e-4));
    const double inf = std::numeric_limits<double>::infinity(),
                 nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(almost_equal(inf, inf));
    CHECK(!almost_equal(inf, -inf));
    CHECK(!almost_equal(inf, 1.0));
    CHECK(!almost_equal(nan, nan));
}
template <class T> void vector_checks() {
    const T tol = T(1e-5);
    const vec<T, 3> x{T(3), T(4), T(0)};
    NEAR(length(x), T(5), tol);
    NEAR(normalize(x), vec<T, 3>{T(0.6), T(0.8), T(0)}, tol);
    CHECK(!try_normalize(vec<T, 3>{}));
    CHECK(normalize(vec<T, 3>{}) == vec<T, 3>{});
    CHECK(!try_normalize(vec<T, 3>{std::numeric_limits<T>::infinity(), T(0), T(0)}));
    const T big = std::numeric_limits<T>::max() / T(2), small = std::numeric_limits<T>::min();
    NEAR(length(normalize(vec<T, 3>{big, big, big})), T(1), tol);
    NEAR(length(normalize(vec<T, 3>{small, small, small})), T(1), tol);
    NEAR(*angle_between(vec<T, 3>{T(1), T(0), T(0)}, vec<T, 3>{T(-1), T(0), T(0)}), pi<T>, tol);
    CHECK(!angle_between(x, vec<T, 3>{}));
    CHECK(reflect(vec<T, 3>{T(1), T(-1), T(0)}, vec<T, 3>{T(0), T(1), T(0)}) ==
          vec<T, 3>{T(1), T(1), T(0)});
    CHECK(refract(vec<T, 3>{T(0), T(-1), T(0)}, vec<T, 3>{T(0), T(1), T(0)}, T(1)) ==
          vec<T, 3>{T(0), T(-1), T(0)});
    CHECK(!refract(normalize(vec<T, 3>{T(1), T(-0.1), T(0)}), vec<T, 3>{T(0), T(1), T(0)}, T(2)));
    CHECK(!project(x, vec<T, 3>{}));
    NEAR(*project(x, vec<T, 3>{T(2), T(0), T(0)}), vec<T, 3>{T(3), T(0), T(0)}, tol);
    for (int i = 0; i < 1000; ++i) {
        const vec<T, 3> u{test::sample<T>(), test::sample<T>(), test::sample<T>()},
            v{test::sample<T>(), test::sample<T>(), test::sample<T>()};
        const auto c = cross(u, v);
        NEAR(dot(c, u), T(0), T(0.001));
        NEAR(dot(c, v), T(0), T(0.001));
        NEAR(length(normalize(u)), T(1), tol);
        NEAR(lerp(u, v, T(0)), u, tol);
        NEAR(lerp(u, v, T(1)), v, tol);
    }
}
CH_TEST(vectors_float_and_double) {
    vector_checks<float>();
    vector_checks<double>();
}
CH_TEST(vector_storage_and_arithmetic) {
    const auto a = constexpr_a, b = constexpr_b;
    vec3f v;
    CHECK(v == vec3f{});
    v.x() = 1;
    v.y() = 2;
    v.z() = 3;
    CHECK(v.data()[2] == 3);
    CHECK(hadamard(a, b) == vec3i{4, 10, 18});
    CHECK(a + b == vec3i{5, 7, 9});
    CHECK(a - b == vec3i{-3, -3, -3});
    CHECK(a * 2 == vec3i{2, 4, 6});
    CHECK((a * 2) / 2 == a);
    CHECK(-a == vec3i{-1, -2, -3});
    CHECK(min(a, b) == a);
    CHECK(max(a, b) == b);
    CHECK(clamp(b, a, vec3i{5}) == vec3i{4, 5, 5});
    CHECK(cross(vec2i{1, 0}, vec2i{0, 1}) == 1);
    CHECK(vec3f(a) == vec3f{1, 2, 3});
    CHECK(vec<double, 1>{2.0}[0] == 2.0);
    CHECK(vec<double, 8>(2.0).data()[7] == 2.0);
}
template <class T, std::size_t N> void matrix_checks() {
    const T tol = T(2e-5);
    for (int trial = 0; trial < 200; ++trial) {
        mat<T, N, N> m;
        vec<T, N> x;
        for (std::size_t i = 0; i < N; ++i) {
            x[i] = test::sample<T>();
            for (std::size_t j = 0; j < N; ++j)
                m(i, j) = test::sample<T>(T(-1), T(1));
            m(i, i) += T(N + 1);
        }
        const auto inv = inverse(m);
        CHECK(inv);
        NEAR(m * (*inv), mat<T, N, N>::identity(), tol);
        NEAR((*inv) * m, mat<T, N, N>::identity(), tol);
        const auto solution = solve(m, m * x);
        CHECK(solution);
        NEAR(*solution, x, tol);
        const auto lu = factor_lu(m);
        CHECK(lu);
        NEAR(lu->solve(m * x), x, tol);
        NEAR(determinant(m), determinant(transpose(m)), T(0.003));
    }
}
CH_TEST(matrix_properties) {
    matrix_checks<float, 2>();
    matrix_checks<float, 3>();
    matrix_checks<float, 4>();
    matrix_checks<float, 6>();
    matrix_checks<double, 1>();
    matrix_checks<double, 2>();
    matrix_checks<double, 3>();
    matrix_checks<double, 4>();
    matrix_checks<double, 8>();
}
CH_TEST(matrix_singular_pivot_and_scale) {
    const auto m = mat3d::from_rows({vec3d{0, 2, 1}, vec3d{1, 0, 3}, vec3d{4, 1, 0}});
    NEAR(determinant(m), 25.0, 1e-12);
    CHECK(inverse(m));
    NEAR(m * (*inverse(m)), mat3d::identity(), 1e-12);
    CHECK(!inverse(mat3d{}));
    CHECK(determinant(mat3d{}) == 0.0);
    auto singular = m;
    singular.set_column(2, singular.column(1));
    CHECK(!inverse(singular));
    for (double scale : {1e-200, 1e-20, 1e20, 1e200}) {
        const auto s = m * scale;
        const auto inv = inverse(s);
        CHECK(inv);
        NEAR(s * (*inv), mat3d::identity(), 1e-12);
    }
    auto invalid = m;
    invalid(0, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK(!inverse(invalid));
    CHECK(std::isnan(determinant(invalid)));
    CHECK(!factor_lu(m, -1.0));
    CHECK(!solve(m, vec3d{std::numeric_limits<double>::infinity(), 0, 0}));
    const mat<double, 2, 3> r = mat<double, 2, 3>::from_rows({vec3d{1, 2, 3}, vec3d{4, 5, 6}});
    CHECK(r.data()[1] == 4);
    CHECK(r.row(1) == vec3d{4, 5, 6});
    CHECK(r.column(2) == vec2d{3, 6});
    CHECK(r * vec3d{1, 0, 0} == vec2d{1, 4});
    CHECK((r * transpose(r))(1, 1) == 77);
    CHECK(outer_product(vec2d{2, 3}, vec3d{4, 5, 6})(1, 2) == 18);
}
