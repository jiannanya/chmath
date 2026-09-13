#include "test.hpp"
#include <chmath/transform/transform.hpp>
#include <type_traits>

using namespace chm;
static_assert(sizeof(quatf) == 16 && sizeof(affine3f) == 48 && sizeof(affine3d) == 96);
static_assert(std::is_trivially_copyable_v<quatf> && std::is_standard_layout_v<affine3f>);
static_assert(transform_point(affine3f{}, vec3f{1, 2, 3}) == vec3f{1, 2, 3});
static_assert(rotate(quatf{}, vec3f{1, 2, 3}) == vec3f{1, 2, 3});

template <class T> void quaternion_properties() {
    const T tol = T(4e-5);
    for (int i = 0; i < 1000; ++i) {
        const auto q =
            from_euler_xyz(vec<T, 3>{test::sample<T>(), test::sample<T>(), test::sample<T>()});
        const auto r =
            from_euler_xyz(vec<T, 3>{test::sample<T>(), test::sample<T>(), test::sample<T>()});
        const vec<T, 3> v{test::sample<T>(), test::sample<T>(), test::sample<T>()};
        NEAR(dot(q, q), T(1), tol);
        NEAR(rotate(q, v), to_matrix(q) * v, tol);
        NEAR(rotate(q * r, v), rotate(q, rotate(r, v)), tol);
        NEAR(rotate(conjugate(q), rotate(q, v)), v, tol);
        const auto inv = inverse(q * T(3));
        CHECK(inv);
        NEAR(((q * T(3)) * (*inv)).coefficients, quat<T>{}.coefficients, tol);
        const auto recovered = from_matrix(to_matrix(q));
        CHECK(recovered);
        NEAR(to_matrix(*recovered), to_matrix(q), tol);
        const auto aa = to_axis_angle(q);
        const auto from_aa = from_axis_angle(aa.axis, aa.angle);
        CHECK(from_aa);
        NEAR(to_matrix(*from_aa), to_matrix(q), tol);
        NEAR(to_matrix(slerp(q, r, T(0))), to_matrix(q), tol);
        NEAR(to_matrix(slerp(q, r, T(1))), to_matrix(r), tol);
        NEAR(to_matrix(slerp(q, -q, T(0.5))), to_matrix(q), tol);
        NEAR(dot(slerp(q, r, T(0.37)), slerp(q, r, T(0.37))), T(1), tol);
        const auto between = rotation_between(v, rotate(q, v));
        CHECK(between);
        NEAR(rotate(*between, normalize(v)), normalize(rotate(q, v)), tol);
    }
}
CH_TEST(quaternion_random_properties) {
    quaternion_properties<float>();
    quaternion_properties<double>();
}
CH_TEST(quaternion_known_and_degenerate) {
    const auto z = from_axis_angle(vec3d{0, 0, 2}, pi<double> / 2);
    CHECK(z);
    NEAR(rotate(*z, vec3d{1, 0, 0}), vec3d{0, 1, 0}, 1e-14);
    CHECK(!from_axis_angle(vec3d{}, 1.0));
    CHECK(!from_axis_angle(vec3d{1, 0, 0}, std::numeric_limits<double>::infinity()));
    CHECK(!inverse(quatd{0, 0, 0, 0}));
    CHECK(!try_normalize(quatd{0, 0, 0, 0}));
    CHECK(normalize(quatd{0, 0, 0, 0}) == quatd{});
    CHECK(!from_matrix(mat3d{2}));
    CHECK(!from_matrix(mat3d{}));
    CHECK(!rotation_between(vec3d{}, vec3d{1, 0, 0}));
    const auto antipodal = rotation_between(vec3d{1, 2, 3}, vec3d{-1, -2, -3});
    CHECK(antipodal);
    NEAR(rotate(*antipodal, vec3d{1, 2, 3}), vec3d{-1, -2, -3}, 1e-12);
    const auto close = rotation_between(vec3d{1, 0, 0}, vec3d{-1, 1e-7, 0});
    CHECK(close);
    NEAR(rotate(*close, vec3d{1, 0, 0}), normalize(vec3d{-1, 1e-7, 0}), 1e-12);
    const auto half = slerp(quatd{}, *z, 0.5);
    NEAR(rotate(half, vec3d{1, 0, 0}), vec3d{std::sqrt(0.5), std::sqrt(0.5), 0}, 1e-12);
    for (const vec3d axis : {vec3d{1, 0, 0}, vec3d{0, 1, 0}, vec3d{0, 0, 1}}) {
        const auto q = from_axis_angle(axis, pi<double>);
        CHECK(q);
        const auto back = from_matrix(to_matrix(*q));
        CHECK(back);
        NEAR(to_matrix(*back), to_matrix(*q), 1e-12);
    }
    quatd q;
    q.x() = 1;
    q.y() = 2;
    q.z() = 3;
    q.w() = 4;
    CHECK(q.data()[3] == 4);
    CHECK(q.imaginary() == vec3d{1, 2, 3});
    const double large = std::numeric_limits<double>::max() / 4;
    CHECK(inverse(quatd{large, large, large, large}));
    NEAR(dot(normalize(q * 1e-300), normalize(q * 1e-300)), 1.0, 1e-14);
}
template <class T> void transform_properties() {
    const T tol = T(1e-4);
    for (int i = 0; i < 500; ++i) {
        trs<T> description;
        description.position = {test::sample<T>(), test::sample<T>(), test::sample<T>()};
        description.orientation =
            from_euler_xyz(vec<T, 3>{test::sample<T>(), test::sample<T>(), test::sample<T>()});
        description.scale = {test::sample<T>(T(0.2), T(5)), test::sample<T>(T(0.2), T(5)),
                             test::sample<T>(T(0.2), T(5))};
        if (i % 2 == 0)
            description.scale[0] *= T(-1);
        const auto a = compose(description);
        const auto inv = inverse(a);
        CHECK(inv);
        const vec<T, 3> v{test::sample<T>(), test::sample<T>(), test::sample<T>()};
        NEAR(transform_point(*inv, transform_point(a, v)), v, tol);
        NEAR(to_matrix(a * (*inv)), mat<T, 4, 4>::identity(), tol);
        const auto dec = decompose(a);
        CHECK(dec);
        NEAR(compose(*dec).matrix, a.matrix, tol);
        const auto m = to_matrix(a);
        CHECK(to_affine(m));
        NEAR(to_affine(m)->matrix, a.matrix, tol);
        CHECK(transform_point(m, v));
        NEAR(*transform_point(m, v), transform_point(a, v), tol);
        NEAR(transform_vector(m, v), transform_vector(a, v), tol);
        const auto n = transform_normal(a, vec<T, 3>{T(0), T(1), T(0)});
        CHECK(n);
        NEAR(dot(*n, transform_vector(a, vec<T, 3>{T(1), T(0), T(0)})), T(0), tol);
    }
}
CH_TEST(affine_properties) {
    transform_properties<float>();
    transform_properties<double>();
}
CH_TEST(affine_failure_and_composition) {
    const auto a = translation(vec3d{1, 2, 3}) * scaling(vec3d{2, 3, 4});
    CHECK(transform_point(a, vec3d{1, 1, 1}) == vec3d{3, 5, 7});
    CHECK(transform_vector(a, vec3d{1, 1, 1}) == vec3d{2, 3, 4});
    auto shear = a;
    shear.matrix(0, 1) = 1;
    CHECK(!decompose(shear));
    CHECK(inverse(shear));
    const auto zero = scaling(vec3d{0, 1, 1});
    CHECK(!inverse(zero));
    CHECK(!normal_matrix(zero));
    CHECK(!decompose(zero));
    CHECK(!transform_normal(a, vec3d{}));
    CHECK(!transform_point(mat4d{}, vec3d{1, 2, 3}));
    auto m = mat4d::identity();
    m(3, 0) = 1;
    CHECK(!to_affine(m));
    CHECK(!look_at(vec3d{}, vec3d{}, vec3d{0, 1, 0}));
    CHECK(!look_at(vec3d{}, vec3d{0, 1, 0}, vec3d{0, 1, 0}));
    for (auto hand : {handedness::right, handedness::left}) {
        const auto view = look_at(vec3d{1, 2, 3}, vec3d{}, vec3d{0, 1, 0}, hand);
        CHECK(view);
        NEAR(transform_point(*view, vec3d{1, 2, 3}), vec3d{}, 1e-14);
        const auto target = transform_point(*view, vec3d{});
        CHECK(hand == handedness::right ? target[2] < 0 : target[2] > 0);
    }
}
template <class T> void projection_properties() {
    const T near = T(0.5), far = T(100), tol = T(3e-4);
    const viewport<T> vp{T(10), T(20), T(1280), T(720), T(0.2), T(0.9)};
    for (auto hand : {handedness::right, handedness::left})
        for (auto range : {depth_range::zero_to_one, depth_range::minus_one_to_one})
            for (auto depth : {depth_direction::forward, depth_direction::reverse}) {
                const T sign = hand == handedness::right ? T(-1) : T(1);
                T dn = range == depth_range::zero_to_one ? T(0) : T(-1), df = T(1);
                if (depth == depth_direction::reverse)
                    std::swap(dn, df);
                const auto p = perspective(pi<T> / T(2), T(2), near, far, hand, range, depth);
                CHECK(p);
                NEAR((*transform_point(*p, vec<T, 3>{T(0), T(0), sign * near}))[2], dn, tol);
                NEAR((*transform_point(*p, vec<T, 3>{T(0), T(0), sign * far}))[2], df, tol);
                NEAR((*transform_point(*p, vec<T, 3>{near * T(2), near, sign * near})),
                     vec<T, 3>{T(1), T(1), dn}, tol);
                const auto inf =
                    perspective(pi<T> / T(2), T(2), near, std::numeric_limits<T>::infinity(), hand,
                                range, depth);
                CHECK(inf);
                NEAR((*transform_point(*inf, vec<T, 3>{T(0), T(0), sign * near}))[2], dn, tol);
                NEAR((*transform_point(*inf, vec<T, 3>{T(0), T(0), sign * T(1e7)}))[2], df, tol);
                const auto ortho =
                    orthographic(T(-4), T(4), T(-2), T(2), near, far, hand, range, depth);
                CHECK(ortho);
                NEAR(*transform_point(*ortho, vec<T, 3>{T(-4), T(-2), sign * near}),
                     vec<T, 3>{T(-1), T(-1), dn}, tol);
                NEAR(*transform_point(*ortho, vec<T, 3>{T(4), T(2), sign * far}),
                     vec<T, 3>{T(1), T(1), df}, tol);
                const auto inv = inverse(*p);
                CHECK(inv);
                for (int i = 0; i < 100; ++i) {
                    const vec<T, 3> point{test::sample<T>(T(-2), T(2)),
                                          test::sample<T>(T(-2), T(2)),
                                          sign * test::sample<T>(T(1), T(20))};
                    const auto screen = project(point, *p, vp, range);
                    CHECK(screen);
                    const auto recovered = unproject(*screen, *inv, vp, range);
                    CHECK(recovered);
                    NEAR(*recovered, point, tol);
                }
            }
}
CH_TEST(projection_all_conventions) {
    projection_properties<float>();
    projection_properties<double>();
}
CH_TEST(projection_rejects_invalid) {
    CHECK(!perspective(0.0, 1.0, 0.1, 100.0));
    CHECK(!perspective(pi<double>, 1.0, 0.1, 100.0));
    CHECK(!perspective(1.0, 0.0, 0.1, 100.0));
    CHECK(!perspective(1.0, 1.0, 0.0, 100.0));
    CHECK(!perspective(1.0, 1.0, 1.0, 0.1));
    CHECK(!perspective(1.0, 1.0, 0.1, std::numeric_limits<double>::quiet_NaN()));
    CHECK(!orthographic(0.0, 0.0, -1.0, 1.0, 0.1, 100.0));
    CHECK(!orthographic(-1.0, 1.0, 0.0, 0.0, 0.1, 100.0));
    CHECK(!project(vec3d{}, mat4d::identity(), viewport<double>{}));
    CHECK(!unproject(vec3d{}, mat4d::identity(), viewport<double>{}));
}
