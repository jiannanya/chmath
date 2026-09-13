#include <atomic>
#include <chmath/chmath.hpp>
#include <cstdio>
#include <cstdlib>
#include <new>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
std::atomic<std::size_t> allocation_count{0};
void *allocate(std::size_t size) {
    allocation_count.fetch_add(1, std::memory_order_relaxed);
    if (void *p = std::malloc(size == 0 ? 1 : size))
        return p;
    throw std::bad_alloc{};
}
void *aligned_allocate(std::size_t size, std::size_t alignment) {
    allocation_count.fetch_add(1, std::memory_order_relaxed);
#if defined(_WIN32)
    void *p = _aligned_malloc(size == 0 ? 1 : size, alignment);
#else
    const std::size_t rounded = ((size == 0 ? 1 : size) + alignment - 1) / alignment * alignment;
    void *p = std::aligned_alloc(alignment, rounded);
#endif
    if (!p)
        throw std::bad_alloc{};
    return p;
}
void aligned_free(void *p) noexcept {
#if defined(_WIN32)
    _aligned_free(p);
#else
    std::free(p);
#endif
}
} // namespace
void *operator new(std::size_t size) {
    return allocate(size);
}
void *operator new[](std::size_t size) {
    return allocate(size);
}
void operator delete(void *p) noexcept {
    std::free(p);
}
void operator delete[](void *p) noexcept {
    std::free(p);
}
void operator delete(void *p, std::size_t) noexcept {
    std::free(p);
}
void operator delete[](void *p, std::size_t) noexcept {
    std::free(p);
}
void *operator new(std::size_t size, std::align_val_t alignment) {
    return aligned_allocate(size, static_cast<std::size_t>(alignment));
}
void *operator new[](std::size_t size, std::align_val_t alignment) {
    return aligned_allocate(size, static_cast<std::size_t>(alignment));
}
void operator delete(void *p, std::align_val_t) noexcept {
    aligned_free(p);
}
void operator delete[](void *p, std::align_val_t) noexcept {
    aligned_free(p);
}
void operator delete(void *p, std::size_t, std::align_val_t) noexcept {
    aligned_free(p);
}
void operator delete[](void *p, std::size_t, std::align_val_t) noexcept {
    aligned_free(p);
}

int main() {
    using namespace chm;
    const auto before = allocation_count.load();
    double checksum = 0;
    for (int i = 0; i < 100; ++i) {
        const double t = static_cast<double>(i) / 100;
        const auto q = from_euler_xyz(vec3d{t, 2 * t, 3 * t});
        const auto model = compose(trs<double>{{1, 2, 3}, q, {2, 3, 4}});
        const auto inv = inverse(model);
        if (!inv)
            return 1;
        const auto eigen = symmetric_eigen(transpose(model.linear()) * model.linear());
        if (!eigen || !eigen->converged)
            return 2;
        const auto qr = factor_qr(model.linear());
        if (!qr)
            return 3;
        const auto lu = factor_lu(model.linear());
        if (!lu)
            return 4;
        const auto solved = solve(model.linear(), vec3d{1, 2, 3});
        if (!solved)
            return 5;
        const auto hit = intersect(ray3d{{0, 0, 4}, {0, 0, -1}}, sphered{{}, 1});
        if (!hit)
            return 6;
        const std::array control{vec3d{0, 0, 0}, vec3d{1, 2, 0}, vec3d{3, 1, 0}, vec3d{4, 0, 0}};
        const std::array<double, 8> knots{0, 0, 0, 0, 1, 1, 1, 1};
        const std::array<double, 4> weights{1, 2, 3, 1};
        const auto spline = bspline_view<double, 3, 3>::create(control, knots);
        const auto rational = nurbs_view<double, 3, 3>::create(control, weights, knots);
        if (!spline || !rational)
            return 7;
        std::array<float, 8> x{}, y{}, z{};
        x.fill(static_cast<float>(t));
        const soa3<float> points{x, y, z};
        if (!transform_points(affine3f{}, points.as_const(), points))
            return 8;
        const auto integral = integrate([](double v) { return v * v; }, 0.0, 1.0);
        if (!integral)
            return 9;
        const auto packed_matrix = mat4f(to_matrix(model));
        const auto matrix_product = packed_matrix * packed_matrix;
        mat4d balanced;
        balanced(0, 0) = balanced(1, 1) = 1e200;
        balanced(2, 2) = balanced(3, 3) = 1e-200;
        const auto extreme_plane =
            plane_from_points(vec3d{}, vec3d{1e200, 0, 0}, vec3d{0, 1e200, 0});
        const auto scaled_integral = integrate([](double) { return 1e308; }, 0.0, 1e-308);
        if (!extreme_plane || !scaled_integral ||
            overlaps(sphered{{1e308, 0, 0}, 9e307}, sphered{{-1e308, 0, 0}, 9e307}))
            return 11;
        const std::array tiny_control{vec2d{}, vec2d{1e-310, 2e-310}};
        const std::array<double, 4> tiny_knots{0, 0, 1e-310, 1e-310};
        const auto tiny_spline = bspline_view<double, 2, 1>::create(tiny_control, tiny_knots);
        if (!tiny_spline)
            return 12;
        const auto tiny_derivative = tiny_spline->derivative(5e-311);
        if (!tiny_derivative)
            return 13;
        const bezier<double, 3, 3> curve{control};
        const auto polynomial = cubic_polynomial<double, 3>::from_bezier(curve);
        const auto joint = curve.evaluate_with_derivative(t);
        const auto multi =
            solve(model.linear(), mat<double, 3, 4>::from_rows(
                                      {vec4d{1, 2, 3, 4}, vec4d{2, 3, 4, 5}, vec4d{3, 4, 5, 6}}));
        const auto prepared = prepared_ray<double>::create(ray3d{{0, 0, 4}, {0, 0, -1}});
        const auto fr = extract_frustum(*perspective(1.0, 1.5, 0.1, 100.0));
        if (!polynomial || !multi || !prepared || !fr)
            return 14;
        const std::array<aabb3d, 2> boxes{aabb3d{{-1, -1, -1}, {1, 1, 1}},
                                          aabb3d{{10, 10, 10}, {11, 11, 11}}};
        std::array<std::optional<ray_interval<double>>, 2> hits{};
        std::array<containment, 2> visible{};
        std::array<double, 9> dx{}, dy{}, dz{};
        dx.fill(t);
        dy.fill(1.0);
        dz.fill(2.0);
        const soa3<double> batch{dx, dy, dz};
        std::array<vec3d, 3> aos{vec3d{t, 1, 2}, vec3d{1, 2, 3}, vec3d{2, 3, 4}};
        if (!intersect_many(*prepared, std::span<const aabb3d>{boxes},
                            std::span<std::optional<ray_interval<double>>>{hits}) ||
            !classify_batch(*fr, std::span<const aabb3d>{boxes}, std::span<containment>{visible}) ||
            !normalize_vectors(batch.as_const(), batch) ||
            !rotate_vectors(q, batch.as_const(), batch) ||
            !transform_vectors(model, batch.as_const(), batch) ||
            !cross_batch(batch.as_const(), batch.as_const(), batch) ||
            !dot_batch(batch.as_const(), batch.as_const(), std::span<double>{dx}) ||
            !normalize_vectors(std::span<const vec3d>{aos}, std::span<vec3d>{aos}) ||
            !transform_vectors(model, std::span<const vec3d>{aos}, std::span<vec3d>{aos}))
            return 15;
        checksum += polynomial->evaluate_with_derivative(t).position[0] + joint.derivative[0] +
                    (*multi)(0, 0) + dx[0] + aos[0][0];
        checksum += double(matrix_product(0, 0)) + determinant(balanced) +
                    extreme_plane->normal[2] + scaled_integral->value + (*tiny_derivative)[0];
        checksum += length(transform_point(*inv, vec3d{t, 1, 2})) + eigen->values[0] +
                    (*spline->evaluate(t))[0] + (*rational->evaluate(t))[1] + integral->value +
                    (*solved)[0] + hit->enter + x[0];
    }
    const auto allocations = allocation_count.load() - before;
    std::printf("core operations: %zu heap allocations, checksum=%.9f\n", allocations, checksum);
    return allocations == 0 && is_finite(checksum) && checksum > 0 ? 0 : 10;
}
