#pragma once

#include <chmath/geometry/primitives.hpp>
#include <chmath/numeric/numeric.hpp>

namespace chm {

template <floating T> struct ray_interval {
    T enter{}, exit{};
};
namespace detail {
template <floating T> [[nodiscard]] constexpr bool valid_interval(T lo, T hi) noexcept {
    return is_finite(lo) && hi == hi && lo <= hi;
}
template <floating T> [[nodiscard]] inline T difference_quotient(T a, T b, T divisor) noexcept {
    const T difference = a - b;
    return is_finite(difference) ? difference / divisor : a / divisor - b / divisor;
}
} // namespace detail
// Closed intervals: touching counts as a hit. Direction need not be normalized.
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<ray_interval<T>>
intersect(const ray<T, N> &r, const aabb<T, N> &b, T t_min = T(0),
          T t_max = std::numeric_limits<T>::infinity()) noexcept {
    if (!r.valid() || !b.valid() || !detail::valid_interval(t_min, t_max))
        return std::nullopt;
    for (std::size_t i = 0; i < N; ++i) {
        if (r.direction[i] == T(0)) {
            if (r.origin[i] < b.lower[i] || r.origin[i] > b.upper[i])
                return std::nullopt;
            continue;
        }
        T a = detail::difference_quotient(b.lower[i], r.origin[i], r.direction[i]),
          c = detail::difference_quotient(b.upper[i], r.origin[i], r.direction[i]);
        if (a > c)
            std::swap(a, c);
        t_min = std::max(t_min, a);
        t_max = std::min(t_max, c);
        if (t_min > t_max)
            return std::nullopt;
    }
    if (!is_finite(t_min))
        return std::nullopt;
    return ray_interval<T>{t_min, t_max};
}
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<ray_interval<T>>
intersect(const ray<T, N> &r, const sphere<T, N> &s, T t_min = T(0),
          T t_max = std::numeric_limits<T>::infinity()) noexcept {
    if (!r.valid() || !s.valid() || !detail::valid_interval(t_min, t_max))
        return std::nullopt;
    auto oc = r.origin - s.center;
    T spatial = s.radius, ds{};
    const bool scaled_difference = !is_finite(oc);
    for (std::size_t i = 0; i < N; ++i) {
        spatial = scaled_difference
                      ? std::max({spatial, std::abs(r.origin[i]), std::abs(s.center[i])})
                      : std::max(spatial, std::abs(oc[i]));
        ds = std::max(ds, std::abs(r.direction[i]));
    }
    if (spatial == T(0)) {
        if (t_min <= T(0) && t_max >= T(0))
            return ray_interval<T>{T(0), T(0)};
        return std::nullopt;
    }
    const auto d = r.direction / ds;
    oc = scaled_difference ? r.origin / spatial - s.center / spatial : oc / spatial;
    const T radius = s.radius / spatial;
    const auto roots = solve_quadratic(dot(d, d), T(2) * dot(oc, d), dot(oc, oc) - radius * radius);
    if (!roots || roots->count == 0)
        return std::nullopt;
    const T a = multiply_divide(roots->values[0], spatial, ds),
            b = multiply_divide(roots->values[roots->count - 1], spatial, ds);
    t_min = std::max(t_min, a);
    t_max = std::min(t_max, b);
    if (t_min > t_max || !is_finite(t_min))
        return std::nullopt;
    return ray_interval<T>{t_min, t_max};
}
template <floating T>
[[nodiscard]] inline std::optional<T> intersect(const ray<T> &r, const plane<T> &p, T t_min = T(0),
                                                T t_max = std::numeric_limits<T>::infinity(),
                                                T tolerance = epsilon<T>) noexcept {
    if (!r.valid() || !is_finite(p.normal) || !is_finite(p.offset) ||
        !detail::valid_interval(t_min, t_max) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    const T ds =
        std::max({std::abs(r.direction[0]), std::abs(r.direction[1]), std::abs(r.direction[2])});
    const auto d = r.direction / ds;
    const T denom = dot(p.normal, d);
    if (std::abs(denom) <= tolerance * length(p.normal) * length(d))
        return std::nullopt;
    const T t = (-p.signed_distance(r.origin) / denom) / ds;
    if (!is_finite(t) || t < t_min || t > t_max)
        return std::nullopt;
    return t;
}
template <floating T> struct triangle_hit {
    T t{};
    vec<T, 3> barycentric{}; // Weights of a,b,c.
    vec<T, 3> normal{};
    bool front_face = false;
};
template <floating T>
[[nodiscard]] inline std::optional<triangle_hit<T>>
intersect(const ray<T> &r, const triangle<T> &tri, T t_min = T(0),
          T t_max = std::numeric_limits<T>::infinity(), bool cull_back_face = false,
          T tolerance = epsilon<T>) noexcept {
    if (!r.valid() || !is_finite(tri.a) || !is_finite(tri.b) || !is_finite(tri.c) ||
        !detail::valid_interval(t_min, t_max) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    auto e1 = tri.b - tri.a, e2 = tri.c - tri.a;
    T es{}, ds{};
    for (std::size_t i = 0; i < 3; ++i) {
        es = std::max({es, std::abs(e1[i]), std::abs(e2[i])});
        ds = std::max(ds, std::abs(r.direction[i]));
    }
    if (es == T(0) || !is_finite(es))
        return std::nullopt;
    e1 /= es;
    e2 /= es;
    const auto d = r.direction / ds, normal = cross(e1, e2), p = cross(d, e2);
    const T det = dot(e1, p), threshold = tolerance * length(normal) * length(d);
    if (cull_back_face ? det <= threshold : std::abs(det) <= threshold)
        return std::nullopt;
    const auto tv = (r.origin - tri.a) / es;
    const T u = dot(tv, p) / det;
    if (!(u >= T(0) && u <= T(1)))
        return std::nullopt;
    const auto q = cross(tv, e1);
    const T v = dot(d, q) / det;
    if (!(v >= T(0) && u + v <= T(1)))
        return std::nullopt;
    const T t = multiply_divide(dot(e2, q) / det, es, ds);
    if (!is_finite(t) || t < t_min || t > t_max)
        return std::nullopt;
    return triangle_hit<T>{t, {T(1) - u - v, u, v}, normalize(normal), det > T(0)};
}
template <floating T>
[[nodiscard]] inline std::optional<ray_interval<T>>
intersect(const ray<T> &r, const obb<T> &box, T t_min = T(0),
          T t_max = std::numeric_limits<T>::infinity()) noexcept {
    if (!box.valid())
        return std::nullopt;
    const auto inv = transpose(box.axes);
    return intersect(ray<T>{inv * (r.origin - box.center), inv * r.direction},
                     aabb<T>{-box.half_extent, box.half_extent}, t_min, t_max);
}
// Separating-axis theorem over 3+3+9 axes. Margin is a world-space contact tolerance.
template <floating T>
[[nodiscard]] inline bool overlaps(const obb<T> &a, const obb<T> &b, T margin = T(0)) noexcept {
    if (!a.valid() || !b.valid() || !is_finite(margin) || margin < T(0))
        return false;
    const auto delta = b.center - a.center;
    const auto separated = [&](const vec<T, 3> &axis) {
        const auto n = try_normalize(axis);
        if (!n)
            return false;
        T radius{};
        for (std::size_t i = 0; i < 3; ++i)
            radius += a.half_extent[i] * std::abs(dot(a.axes.column(i), *n)) +
                      b.half_extent[i] * std::abs(dot(b.axes.column(i), *n));
        return std::abs(dot(delta, *n)) > radius + margin;
    };
    for (std::size_t i = 0; i < 3; ++i)
        if (separated(a.axes.column(i)) || separated(b.axes.column(i)))
            return false;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            if (separated(cross(a.axes.column(i), b.axes.column(j))))
                return false;
    return true;
}

} // namespace chm
