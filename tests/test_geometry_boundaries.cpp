#include "boundary_helpers.hpp"
#include "test.hpp"
using namespace chm;

CH_TEST(geometry_empty_box_never_overlaps) {
    CHECK(!overlaps(aabb3d{}, aabb3d{}));
}
CH_TEST(geometry_inverted_box_invalid) {
    const aabb3d b{{1, 0, 0}, {0, 1, 1}};
    CHECK(b.empty());
    CHECK(!b.valid());
    CHECK(!b.contains(vec3d{}));
}
CH_TEST(geometry_nan_box_invalid) {
    const aabb3d b{{boundary::nan, 0, 0}, {1, 1, 1}};
    CHECK(b.empty());
    CHECK(!b.valid());
}
CH_TEST(geometry_box_huge_midpoint) {
    const aabb3d b{{-boundary::maximum, 0, 0}, {boundary::maximum, 0, 0}};
    CHECK(b.center() == vec3d{});
}
CH_TEST(geometry_box_point_degenerate) {
    const aabb3d b{{1, 2, 3}, {1, 2, 3}};
    CHECK(b.valid());
    CHECK(b.contains(vec3d{1, 2, 3}));
    CHECK(volume(b) == 0);
}
CH_TEST(geometry_ray_negative_interval) {
    const auto h =
        intersect(ray3d{{2, 0, 0}, {1, 0, 0}}, aabb3d{{-1, -1, -1}, {1, 1, 1}}, -10.0, 0.0);
    CHECK(h);
    CHECK(h->enter == -3 && h->exit == -1);
}
CH_TEST(geometry_ray_box_tangent_single_parameter) {
    const auto h = intersect(ray3d{{-2, 0, 0}, {1, 0, 0}}, aabb3d{{0, 0, 0}, {0, 0, 0}});
    CHECK(h);
    CHECK(h->enter == 2 && h->exit == 2);
}
CH_TEST(geometry_ray_box_parallel_on_face) {
    CHECK(intersect(ray3d{{1, 0, 2}, {-0.0, 0, -1}}, aabb3d{{-1, -1, -1}, {1, 1, 1}}));
}
CH_TEST(geometry_ray_box_parallel_outside_face) {
    CHECK(!intersect(ray3d{{1.01, 0, 2}, {-0.0, 0, -1}}, aabb3d{{-1, -1, -1}, {1, 1, 1}}));
}
CH_TEST(geometry_ray_box_nan_interval) {
    CHECK(!intersect(ray3d{{0, 0, 2}, {0, 0, -1}}, aabb3d{{-1, -1, -1}, {1, 1, 1}}, 0.0,
                     boundary::nan));
}
CH_TEST(geometry_ray_sphere_surface_outward) {
    const auto h = intersect(ray3d{{1, 0, 0}, {1, 0, 0}}, sphered{{}, 1});
    CHECK(h);
    CHECK(h->enter == 0 && h->exit == 0);
}
CH_TEST(geometry_ray_sphere_surface_inward) {
    const auto h = intersect(ray3d{{1, 0, 0}, {-1, 0, 0}}, sphered{{}, 1});
    CHECK(h);
    CHECK(h->enter == 0);
    NEAR(h->exit, 2.0, 1e-14);
}
CH_TEST(geometry_ray_sphere_tiny_scale) {
    const auto h = intersect(ray3d{{0, 0, 3e-200}, {0, 0, -1e-200}}, sphered{{}, 1e-200});
    CHECK(h);
    NEAR(h->enter, 2.0, 1e-13);
    NEAR(h->exit, 4.0, 1e-13);
}
CH_TEST(geometry_ray_sphere_zero_radius_hit) {
    const auto h = intersect(ray3d{{2, 0, 0}, {-1, 0, 0}}, sphered{{}, 0});
    CHECK(h);
    NEAR(h->enter, 2.0, 1e-13);
}
CH_TEST(geometry_sphere_overlap_overflow_separated) {
    CHECK(!overlaps(sphered{{1e308, 0, 0}, 9e307}, sphered{{-1e308, 0, 0}, 9e307}));
}
CH_TEST(geometry_sphere_overlap_overflow_touching) {
    CHECK(overlaps(sphered{{1e308, 0, 0}, 1e308}, sphered{{-1e308, 0, 0}, 1e308}));
}
CH_TEST(geometry_sphere_aabb_corner_miss) {
    CHECK(!overlaps(sphered{{2, 2, 2}, 1}, aabb3d{{-1, -1, -1}, {1, 1, 1}}));
}
CH_TEST(geometry_plane_offset_normalization) {
    const auto p = make_plane(vec3d{0, 0, 10}, -30.0);
    CHECK(p);
    CHECK(p->offset == -3);
    CHECK(p->normal == vec3d{0, 0, 1});
}
CH_TEST(geometry_plane_from_tiny_triangle) {
    const auto p = plane_from_points(vec3d{}, vec3d{1e-200, 0, 0}, vec3d{0, 1e-200, 0});
    CHECK(p);
    NEAR(p->normal, vec3d{0, 0, 1}, 1e-14);
}
CH_TEST(geometry_plane_from_huge_triangle) {
    const auto p = plane_from_points(vec3d{}, vec3d{1e200, 0, 0}, vec3d{0, 1e200, 0});
    CHECK(p);
    NEAR(p->normal, vec3d{0, 0, 1}, 1e-14);
    const auto wide =
        plane_from_points(vec3d{-1e308, 0, 0}, vec3d{1e308, 0, 0}, vec3d{0, 1e308, 0});
    CHECK(wide);
    NEAR(wide->normal, vec3d{0, 0, 1}, 1e-14);
}
CH_TEST(geometry_ray_plane_coplanar_no_unique_hit) {
    CHECK(!intersect(ray3d{{1, 2, 0}, {1, 0, 0}}, plane<double>{{0, 0, 1}, 0}));
}
CH_TEST(geometry_ray_plane_negative_tolerance) {
    CHECK(!intersect(ray3d{{0, 0, 1}, {0, 0, -1}}, plane<double>{{0, 0, 1}, 0}, 0.0, 10.0, -1.0));
}
CH_TEST(geometry_triangle_vertex_hit) {
    const triangle<double> t{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const auto h = intersect(ray3d{{0, 0, 1}, {0, 0, -1}}, t);
    CHECK(h);
    CHECK(h->barycentric == vec3d{1, 0, 0});
}
CH_TEST(geometry_triangle_edge_hit) {
    const triangle<double> t{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const auto h = intersect(ray3d{{0.5, 0.5, 1}, {0, 0, -1}}, t);
    CHECK(h);
    NEAR(h->barycentric, vec3d{0, 0.5, 0.5}, 1e-14);
}
CH_TEST(geometry_triangle_outside_barycentrics) {
    const triangle<double> t{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    CHECK(!intersect(ray3d{{0.6, 0.6, 1}, {0, 0, -1}}, t));
}
CH_TEST(geometry_triangle_winding_flips_normal) {
    const triangle<double> t{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}};
    const auto h = intersect(ray3d{{0.2, 0.2, 1}, {0, 0, -1}}, t);
    CHECK(h);
    CHECK(!h->front_face);
    CHECK(h->normal == vec3d{0, 0, -1});
}
CH_TEST(geometry_triangle_tiny_scale_hit) {
    const triangle<double> t{{0, 0, 0}, {1e-200, 0, 0}, {0, 1e-200, 0}};
    const auto h = intersect(ray3d{{2e-201, 2e-201, 1e-200}, {0, 0, -1e-200}}, t);
    CHECK(h);
    NEAR(h->t, 1.0, 1e-13);
}
CH_TEST(geometry_triangle_closest_collapsed_point) {
    const triangle<double> t{{1, 2, 3}, {1, 2, 3}, {1, 2, 3}};
    CHECK(closest_point(t, vec3d{}) == vec3d{1, 2, 3});
}
CH_TEST(geometry_segment_extreme_endpoints) {
    const auto p = closest_point(segment<double>{{-1e308, 0, 0}, {1e308, 0, 0}}, vec3d{0, 1, 0});
    CHECK(is_finite(p));
    NEAR(p, vec3d{}, 1e-13);
    const segment<double> tiny_segment{{0, 0, 0}, {1e-308, 0, 0}};
    CHECK(closest_point(tiny_segment, vec3d{1e308, 1e308, 0}) == tiny_segment.b);
    CHECK(closest_point(tiny_segment, vec3d{-1e308, 1e308, 0}) == tiny_segment.a);
}
CH_TEST(geometry_obb_nonorthogonal_axes_rejected) {
    obb<double> b{{}, {1, 1, 1}, mat3d::identity()};
    b.axes(0, 1) = 0.1;
    CHECK(!b.valid());
    CHECK(!overlaps(b, b));
}
CH_TEST(geometry_obb_zero_extent_point_contact) {
    const obb<double> a{{}, {}, mat3d::identity()}, b{{}, {1, 1, 1}, mat3d::identity()};
    CHECK(overlaps(a, b));
}
CH_TEST(geometry_obb_negative_margin_rejected) {
    const obb<double> b{{}, {1, 1, 1}, mat3d::identity()};
    CHECK(!overlaps(b, b, -1.0));
}
CH_TEST(geometry_frustum_nonfinite_matrix_rejected) {
    auto m = mat4d::identity();
    m(1, 2) = boundary::inf;
    CHECK(!extract_frustum(m));
}
CH_TEST(geometry_frustum_disabled_infinite_plane) {
    const auto p = perspective(1.0, 1.0, 0.1, boundary::inf);
    CHECK(p);
    const auto f = extract_frustum(*p);
    CHECK(f);
    CHECK(f->active_mask != 63);
    CHECK(classify(*f, vec3d{0, 0, -1e100}) == containment::inside);
}
CH_TEST(geometry_transform_bounds_mirrored_axes) {
    const auto b = transform_bounds(scaling(vec3d{-2, 3, -4}), aabb3d{{1, 2, 3}, {4, 5, 6}});
    CHECK(b.lower == vec3d{-8, 6, -24});
    CHECK(b.upper == vec3d{-2, 15, -12});
}
