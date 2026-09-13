#include "boundary_helpers.hpp"
#include "test.hpp"
using namespace chm;

CH_TEST(rotation_axis_angle_negative_axis) {
    const auto q = from_axis_angle(vec3d{0, 0, -1}, pi<double> / 2);
    CHECK(q);
    NEAR(rotate(*q, vec3d{1, 0, 0}), vec3d{0, -1, 0}, 1e-14);
}
CH_TEST(rotation_axis_angle_denormal_axis) {
    const auto q = from_axis_angle(vec3d{0, boundary::tiny, 0}, pi<double>);
    CHECK(q);
    NEAR(rotate(*q, vec3d{1, 0, 0}), vec3d{-1, 0, 0}, 1e-14);
}
CH_TEST(rotation_axis_angle_huge_axis) {
    const auto q = from_axis_angle(vec3d{boundary::maximum, 0, 0}, pi<double> / 2);
    CHECK(q);
    NEAR(rotate(*q, vec3d{0, 1, 0}), vec3d{0, 0, 1}, 1e-14);
}
CH_TEST(rotation_identity_axis_angle_canonical) {
    const auto a = to_axis_angle(quatd{});
    CHECK(a.axis == vec3d{1, 0, 0});
    CHECK(a.angle == 0);
}
CH_TEST(rotation_negative_identity_same_rotation) {
    NEAR(to_matrix(quatd{0, 0, 0, -1}), mat3d::identity(), 1e-14);
    CHECK(to_axis_angle(quatd{0, 0, 0, -1}).angle == 0);
}
CH_TEST(rotation_inverse_large_nonunit) {
    const quatd q{1e200, -2e200, 3e200, 4e200};
    const auto inv = inverse(q);
    CHECK(inv);
    NEAR((q * (*inv)).coefficients, quatd{}.coefficients, 1e-13);
}
CH_TEST(rotation_inverse_small_nonunit) {
    const quatd q{1e-200, -2e-200, 3e-200, 4e-200};
    const auto inv = inverse(q);
    CHECK(inv);
    NEAR((q * (*inv)).coefficients, quatd{}.coefficients, 1e-13);
}
CH_TEST(rotation_nlerp_antipodes) {
    const auto q = from_euler_xyz(vec3d{1, 2, 3});
    NEAR(to_matrix(nlerp(q, -q, 0.5)), to_matrix(q), 1e-13);
}
CH_TEST(rotation_slerp_nearly_equal) {
    const auto q = from_axis_angle(vec3d{1, 0, 0}, 1e-8);
    CHECK(q);
    NEAR(to_matrix(slerp(quatd{}, *q, 0.5)), to_matrix(*from_axis_angle(vec3d{1, 0, 0}, 0.5e-8)),
         1e-14);
}
CH_TEST(rotation_slerp_half_turn_midpoint) {
    const auto q = from_axis_angle(vec3d{0, 0, 1}, pi<double>);
    CHECK(q);
    NEAR(rotate(slerp(quatd{}, *q, 0.5), vec3d{1, 0, 0}), vec3d{0, 1, 0}, 1e-14);
}
CH_TEST(rotation_between_y_antipodes) {
    const auto q = rotation_between(vec3d{0, 1, 0}, vec3d{0, -1, 0});
    CHECK(q);
    NEAR(rotate(*q, vec3d{0, 1, 0}), vec3d{0, -1, 0}, 1e-14);
}
CH_TEST(rotation_between_z_antipodes) {
    const auto q = rotation_between(vec3d{0, 0, 1}, vec3d{0, 0, -1});
    CHECK(q);
    NEAR(rotate(*q, vec3d{0, 0, 1}), vec3d{0, 0, -1}, 1e-14);
}
CH_TEST(rotation_matrix_reflection_rejected) {
    auto m = mat3d::identity();
    m(0, 0) = -1;
    CHECK(!from_matrix(m));
}
CH_TEST(rotation_matrix_shear_rejected) {
    auto m = mat3d::identity();
    m(0, 1) = 0.2;
    CHECK(!from_matrix(m));
}
CH_TEST(rotation_euler_axis_order) {
    const vec3d angles{0.2, 0.4, 0.7};
    const auto qx = from_axis_angle(vec3d{1, 0, 0}, angles[0]),
               qy = from_axis_angle(vec3d{0, 1, 0}, angles[1]),
               qz = from_axis_angle(vec3d{0, 0, 1}, angles[2]);
    CHECK(qx && qy && qz);
    NEAR(to_matrix(from_euler_xyz(angles)), to_matrix(*qz * *qy * *qx), 1e-13);
}
CH_TEST(rotation_constexpr_unit_transform) {
    constexpr auto v = rotate(quatf{}, vec3f{2, 3, 4});
    static_assert(v == vec3f{2, 3, 4});
    CHECK(v == vec3f{2, 3, 4});
}
CH_TEST(affine_translation_composition_order) {
    const auto m = scaling(vec3d{2}) * translation(vec3d{1, 2, 3});
    CHECK(transform_point(m, vec3d{}) == vec3d{2, 4, 6});
}
CH_TEST(affine_direction_ignores_translation) {
    CHECK(transform_vector(translation(vec3d{1e100}), vec3d{1, 2, 3}) == vec3d{1, 2, 3});
}
CH_TEST(affine_negative_y_scale_decomposition) {
    const auto a =
        compose(trs<double>{{1, 2, 3}, from_euler_xyz(vec3d{0.4, 0.7, 0.1}), {2, -3, 4}});
    const auto d = decompose(a);
    CHECK(d);
    CHECK(d->scale[0] < 0);
    NEAR(compose(*d).matrix, a.matrix, 1e-13);
}
CH_TEST(affine_two_negative_scales_decomposition) {
    const auto a = scaling(vec3d{-2, -3, 4});
    const auto d = decompose(a);
    CHECK(d);
    NEAR(compose(*d).matrix, a.matrix, 1e-13);
}
CH_TEST(affine_three_negative_scales_decomposition) {
    const auto a = scaling(vec3d{-2, -3, -4});
    const auto d = decompose(a);
    CHECK(d);
    NEAR(compose(*d).matrix, a.matrix, 1e-13);
}
CH_TEST(affine_tiny_uniform_scale_inverse) {
    const auto a = scaling(vec3d{1e-200});
    const auto i = inverse(a);
    CHECK(i);
    NEAR(transform_point(*i, transform_point(a, vec3d{1, 2, 3})), vec3d{1, 2, 3}, 1e-13);
}
CH_TEST(affine_normal_uses_inverse_transpose) {
    const auto a = scaling(vec3d{2, 1, 1});
    const auto n = transform_normal(a, normalize(vec3d{1, 1, 0}));
    CHECK(n);
    NEAR(*n, normalize(vec3d{0.5, 1, 0}), 1e-14);
}
CH_TEST(affine_to_matrix_last_row_exact) {
    CHECK(to_matrix(translation(vec3d{1, 2, 3})).row(3) == vec4d{0, 0, 0, 1});
}
CH_TEST(affine_homogeneous_zero_w_rejected) {
    auto a = mat4d::identity();
    a(3, 3) = 0;
    CHECK(!transform_point(a, vec3d{}));
}
CH_TEST(affine_homogeneous_negative_w_division) {
    auto a = mat4d::identity();
    a(3, 3) = -2;
    const auto p = transform_point(a, vec3d{2, 4, 6});
    CHECK(p);
    CHECK(*p == vec3d{-1, -2, -3});
}
CH_TEST(camera_look_at_extreme_eye_target) {
    const auto v = look_at(vec3d{1e308, 0, 0}, vec3d{-1e308, 0, 0}, vec3d{0, 1, 0});
    CHECK(v);
    NEAR(transform_point(*v, vec3d{1e308, 0, 0}), vec3d{}, 1e-12);
}
CH_TEST(camera_look_at_scaled_up) {
    const auto a = look_at(vec3d{0, 0, 4}, vec3d{}, vec3d{0, 1, 0}),
               b = look_at(vec3d{0, 0, 4}, vec3d{}, vec3d{0, 1e300, 0});
    CHECK(a && b);
    NEAR(a->matrix, b->matrix, 1e-14);
}
CH_TEST(camera_parallel_up_rejected) {
    CHECK(!look_at(vec3d{0, 0, 1}, vec3d{}, vec3d{0, 0, 2}));
}
CH_TEST(camera_nan_up_rejected) {
    CHECK(!look_at(vec3d{0, 0, 1}, vec3d{}, vec3d{0, boundary::nan, 0}));
}
CH_TEST(projection_equal_near_far_rejected) {
    CHECK(!perspective(1.0, 1.0, 1.0, 1.0));
    CHECK(!orthographic(-1.0, 1.0, -1.0, 1.0, 1.0, 1.0));
}
CH_TEST(projection_negative_aspect_rejected) {
    CHECK(!perspective(1.0, -1.0, 0.1, 100.0));
}
CH_TEST(projection_infinite_fovy_rejected) {
    CHECK(!perspective(boundary::inf, 1.0, 0.1, 100.0));
}
CH_TEST(projection_reverse_infinite_far_coefficients) {
    const auto p = perspective(1.0, 1.0, 0.25, boundary::inf, handedness::right,
                               depth_range::zero_to_one, depth_direction::reverse);
    CHECK(p);
    CHECK((*p)(2, 2) == 0);
    CHECK((*p)(2, 3) == 0.25);
}
CH_TEST(projection_orthographic_negative_near) {
    const auto p = orthographic(-1.0, 1.0, -1.0, 1.0, -2.0, 2.0);
    CHECK(p);
    NEAR(*transform_point(*p, vec3d{0, 0, 2}), vec3d{0, 0, 0}, 1e-14);
}
CH_TEST(viewport_infinite_width_rejected) {
    CHECK(!unproject(vec3d{}, mat4d::identity(), viewport<double>{0, 0, boundary::inf, 1}));
}
CH_TEST(viewport_infinite_depth_rejected) {
    CHECK(!unproject(vec3d{}, mat4d::identity(), viewport<double>{0, 0, 1, 1, 0, boundary::inf}));
}
CH_TEST(viewport_nan_height_rejected) {
    CHECK(!project(vec3d{}, mat4d::identity(), viewport<double>{0, 0, 1, boundary::nan}));
}
CH_TEST(viewport_offset_depth_roundtrip) {
    const viewport<double> v{100, 200, 300, 400, -1, 2};
    const vec3d p{0.25, -0.5, 0.75};
    const auto w = project(p, mat4d::identity(), v);
    CHECK(w);
    const auto r = unproject(*w, mat4d::identity(), v);
    CHECK(r);
    NEAR(*r, p, 1e-14);
}
CH_TEST(viewport_zero_depth_interval_rejected) {
    CHECK(!project(vec3d{}, mat4d::identity(), viewport<double>{0, 0, 1, 1, 0.5, 0.5}));
}
