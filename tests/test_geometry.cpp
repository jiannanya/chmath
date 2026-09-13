#include "test.hpp"
#include <chmath/geometry/frustum.hpp>
#include <chmath/geometry/intersection.hpp>

using namespace chm;

CH_TEST(bounds_and_closest_points) {
    aabb3d b;
    CHECK(b.empty());
    CHECK(!b.contains(vec3d{}));
    CHECK(volume(b) == 0.0);
    b.expand(vec3d{1, 2, 3});
    b.expand(vec3d{-1, -2, -3});
    CHECK(b.valid());
    CHECK(b.center() == vec3d{});
    CHECK(b.extent() == vec3d{2, 4, 6});
    CHECK(surface_area(b) == 88.0);
    CHECK(volume(b) == 48.0);
    CHECK(b.contains(vec3d{1, 2, 3}));
    CHECK(!b.contains(vec3d{2, 2, 3}));
    CHECK(closest_point(b, vec3d{9, -9, 0}) == vec3d{1, -2, 0});
    CHECK(overlaps(b, aabb3d{{1, 2, 3}, {2, 3, 4}}));
    CHECK(!overlaps(b, aabb3d{}));
    CHECK(!overlaps(b, aabb3d{{5, 5, 5}, {6, 6, 6}}));
    b.expand(aabb3d{});
    CHECK(b.lower == vec3d{-1, -2, -3});
    CHECK(closest_point(segment<double>{{0, 0, 0}, {2, 0, 0}}, vec3d{1, 3, 0}) == vec3d{1, 0, 0});
    CHECK(closest_point(segment<double>{{1, 2, 3}, {1, 2, 3}}, vec3d{}) == vec3d{1, 2, 3});
    CHECK(closest_point(segment<double>{{0, 0, 0}, {2, 0, 0}}, vec3d{4, 1, 0}) == vec3d{2, 0, 0});
    CHECK(overlaps(sphered{{0, 0, 0}, 1}, sphered{{2, 0, 0}, 1}));
    CHECK(!overlaps(sphered{{0, 0, 0}, -1}, sphered{{0, 0, 0}, 1}));
    CHECK(overlaps(sphered{{2, 0, 0}, 1}, b));
    CHECK(!overlaps(sphered{{3, 0, 0}, 1}, b));
    const auto p = make_plane(vec3d{0, 2, 0}, -4.0);
    CHECK(p);
    CHECK(p->signed_distance(vec3d{0, 3, 0}) == 1.0);
    CHECK(closest_point(*p, vec3d{1, 5, 3}) == vec3d{1, 2, 3});
    CHECK(!make_plane(vec3d{}, 1.0));
    const auto plane3 = plane_from_points(vec3d{0, 0, 1}, vec3d{1, 0, 1}, vec3d{0, 1, 1});
    CHECK(plane3);
    NEAR(plane3->signed_distance(vec3d{2, 3, 1}), 0.0, 1e-14);
    CHECK(!plane_from_points(vec3d{}, vec3d{1, 1, 1}, vec3d{2, 2, 2}));
    const triangle<double> tri{{0, 0, 0}, {2, 0, 0}, {0, 2, 0}};
    CHECK(area(tri) == 2.0);
    const auto bc = barycentric(tri, vec3d{0.5, 0.5, 9});
    CHECK(bc);
    NEAR(*bc, vec3d{0.5, 0.25, 0.25}, 1e-14);
    NEAR(closest_point(tri, vec3d{0.5, 0.5, 9}), vec3d{0.5, 0.5, 0}, 1e-14);
    NEAR(closest_point(tri, vec3d{2, 2, 3}), vec3d{1, 1, 0}, 1e-14);
    CHECK(closest_point(tri, vec3d{-1, -1, 0}) == vec3d{});
    const triangle<double> line{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
    CHECK(!barycentric(line, vec3d{}));
    CHECK(closest_point(line, vec3d{0.5, 1, 0}) == vec3d{0.5, 0, 0});
    CHECK(!barycentric(tri, vec3d{std::numeric_limits<double>::quiet_NaN(), 0, 0}));
}
template <class T> void ray_tests() {
    const T tol = T(1e-5), inf = std::numeric_limits<T>::infinity();
    const aabb<T> box{{T(-1), T(-1), T(-1)}, {T(1), T(1), T(1)}};
    const sphere<T> ball{{}, T(1)};
    const ray<T> r{{T(0), T(0), T(3)}, {T(0), T(0), T(-2)}};
    const auto bh = intersect(r, box), sh = intersect(r, ball);
    CHECK(bh && sh);
    NEAR(bh->enter, T(1), tol);
    NEAR(bh->exit, T(2), tol);
    NEAR(sh->enter, T(1), tol);
    NEAR(sh->exit, T(2), tol);
    CHECK(!intersect(r, box, T(0), T(0.5)));
    CHECK(!intersect(r, ball, T(0), T(0.5)));
    CHECK(!intersect(ray<T>{}, box));
    CHECK(!intersect(ray<T>{}, ball));
    CHECK(!intersect(r, aabb<T>{}));
    const ray<T> inside{{}, {T(1), T(0), T(0)}};
    CHECK(intersect(inside, box)->enter == T(0));
    CHECK(intersect(inside, ball)->enter == T(0));
    NEAR(intersect(inside, ball)->exit, T(1), tol);
    const ray<T> edge{{T(1), T(1), T(3)}, {T(0), -T(0), T(-1)}};
    CHECK(intersect(edge, box));
    CHECK(!intersect(ray<T>{{T(2), T(0), T(3)}, {T(0), T(0), T(-1)}}, box));
    const auto tangent = intersect(ray<T>{{T(1), T(0), T(3)}, {T(0), T(0), T(-1)}}, ball);
    CHECK(tangent);
    NEAR(tangent->enter, T(3), tol);
    NEAR(tangent->exit, T(3), tol);
    const auto point = intersect(ray<T>{{}, {T(1), T(0), T(0)}}, sphere<T>{{}, T(0)});
    CHECK(point);
    CHECK(point->enter == T(0));
    CHECK(!intersect(ray<T>{{T(0), T(0), T(3)}, {T(0), T(0), T(1)}}, ball));
    CHECK(!intersect(r, box, T(2), T(1)));
    CHECK(!intersect(r, ball, T(2), T(1)));
    const plane<T> p{{T(0), T(0), T(1)}, T(0)};
    const auto ph = intersect(r, p);
    CHECK(ph);
    NEAR(*ph, T(1.5), tol);
    CHECK(!intersect(inside, p));
    CHECK(!intersect(r, p, T(0), T(1)));
    const triangle<T> tri{{T(-1), T(-1), T(0)}, {T(1), T(-1), T(0)}, {T(0), T(1), T(0)}};
    const auto th = intersect(r, tri);
    CHECK(th);
    NEAR(th->t, T(1.5), tol);
    NEAR(th->barycentric, vec<T, 3>{T(0.25), T(0.25), T(0.5)}, tol);
    CHECK(th->front_face);
    NEAR(r.at(th->t),
         tri.a * th->barycentric[0] + tri.b * th->barycentric[1] + tri.c * th->barycentric[2], tol);
    CHECK(intersect(r, tri, T(0), inf, true));
    const ray<T> back{{T(0), T(0), T(-3)}, {T(0), T(0), T(1)}};
    CHECK(intersect(back, tri));
    CHECK(!intersect(back, tri, T(0), inf, true));
    CHECK(!intersect(r, triangle<T>{}));
    CHECK(!intersect(r, tri, T(0), T(1)));
    CHECK(!intersect(r, tri, T(0), inf, false, T(-1)));
    CHECK(!intersect(ray<T>{{inf, T(0), T(0)}, {T(1), T(0), T(0)}}, box));
}
CH_TEST(ray_intersections_float_and_double) {
    ray_tests<float>();
    ray_tests<double>();
}
CH_TEST(bounds_transform_property) {
    const aabb3d box{{-1, -2, -3}, {4, 5, 6}};
    for (int trial = 0; trial < 500; ++trial) {
        affine3d a;
        for (auto &x : a.matrix.elements)
            x = test::sample<double>();
        const auto transformed = transform_bounds(a, box);
        aabb3d reference;
        for (unsigned bits = 0; bits < 8; ++bits) {
            vec3d corner;
            for (unsigned i = 0; i < 3; ++i)
                corner[i] = (bits & (1U << i)) != 0 ? box.upper[i] : box.lower[i];
            reference.expand(transform_point(a, corner));
        }
        NEAR(transformed.lower, reference.lower, 1e-12);
        NEAR(transformed.upper, reference.upper, 1e-12);
    }
    CHECK(transform_bounds(affine3d{}, aabb3d{}).empty());
}
CH_TEST(obb_sat_and_rays) {
    const obb<double> a{{0, 0, 0}, {1, 2, 3}, mat3d::identity()};
    auto b = a;
    CHECK(a.valid());
    CHECK(overlaps(a, b));
    b.center = {2, 0, 0};
    CHECK(overlaps(a, b));
    b.center = {2.1, 0, 0};
    CHECK(!overlaps(a, b));
    CHECK(overlaps(a, b, 0.11));
    b.axes = to_matrix(from_euler_xyz(vec3d{0, 0, pi<double> / 4}));
    CHECK(overlaps(a, b));
    b.center = {10, 10, 10};
    CHECK(!overlaps(a, b));
    b.half_extent[0] = -1;
    CHECK(!b.valid());
    CHECK(!overlaps(a, b));
    CHECK(!intersect(ray3d{{0, 0, 5}, {0, 0, -1}}, b));
    const auto h = intersect(ray3d{{0, 0, 5}, {0, 0, -1}}, a);
    CHECK(h);
    NEAR(h->enter, 2.0, 1e-14);
    NEAR(h->exit, 8.0, 1e-14);
    for (int i = 0; i < 500; ++i) {
        const auto rot = from_euler_xyz(
            vec3d{test::sample<double>(), test::sample<double>(), test::sample<double>()});
        b = a;
        b.axes = to_matrix(rot);
        b.center = {test::sample<double>(), test::sample<double>(), test::sample<double>()};
        CHECK(overlaps(a, b) == overlaps(b, a));
        const auto local = b.axes * vec3d{0, 0, 5} + b.center;
        const auto hit = intersect(ray3d{local, b.axes * vec3d{0, 0, -1}}, b);
        CHECK(hit);
        NEAR(hit->enter, 2.0, 1e-12);
    }
}
CH_TEST(frustum_conventions_and_infinite_far) {
    for (auto hand : {handedness::right, handedness::left})
        for (auto range : {depth_range::zero_to_one, depth_range::minus_one_to_one})
            for (auto depth : {depth_direction::forward, depth_direction::reverse}) {
                const double s = hand == handedness::right ? -1.0 : 1.0;
                const auto p = perspective(pi<double> / 2, 1.0, 1.0, 10.0, hand, range, depth);
                CHECK(p);
                const auto f = extract_frustum(*p, range);
                CHECK(f);
                CHECK(f->active_mask == 63);
                CHECK(classify(*f, vec3d{0, 0, s * 5}) == containment::inside);
                CHECK(classify(*f, vec3d{0, 0, s * 0.5}) == containment::outside);
                CHECK(classify(*f, vec3d{0, 0, s * 11}) == containment::outside);
                CHECK(classify(*f, vec3d{20, 0, s * 5}) == containment::outside);
                CHECK(classify(*f, sphered{{0, 0, s * 5}, 0.1}) == containment::inside);
                CHECK(classify(*f, sphered{{0, 0, s}, 0.5}) == containment::intersecting);
                CHECK(classify(*f, sphered{{0, 0, s * 20}, 0.5}) == containment::outside);
                const vec3d center{0, 0, s * 5};
                CHECK(classify(*f, aabb3d{center - vec3d{0.1}, center + vec3d{0.1}}) ==
                      containment::inside);
                CHECK(classify(*f, aabb3d{}) == containment::outside);
                const auto inf =
                    perspective(pi<double> / 2, 1.0, 1.0, std::numeric_limits<double>::infinity(),
                                hand, range, depth);
                CHECK(inf);
                const auto fi = extract_frustum(*inf, range);
                CHECK(fi);
                CHECK(classify(*fi, vec3d{0, 0, s * 1e12}) == containment::inside);
                for (int i = 0; i < 300; ++i) {
                    const vec3d point{test::sample<double>(), test::sample<double>(),
                                      test::sample<double>()};
                    const auto clip = *p * vec4d{point[0], point[1], point[2], 1};
                    const bool expected =
                        clip[0] >= -clip[3] && clip[0] <= clip[3] && clip[1] >= -clip[3] &&
                        clip[1] <= clip[3] &&
                        clip[2] >= (range == depth_range::zero_to_one ? 0.0 : -clip[3]) &&
                        clip[2] <= clip[3];
                    CHECK((classify(*f, point) == containment::inside) == expected);
                }
            }
    CHECK(!extract_frustum(mat4d{}));
}
