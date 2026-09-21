#pragma once

#include <chmath/transform/transform.hpp>
#include <numeric>

namespace chm {

template <floating T, std::size_t N = 3> struct ray {
    vec<T, N> origin{}, direction{};
    [[nodiscard]] constexpr vec<T, N> at(T t) const noexcept { return origin + direction * t; }
    [[nodiscard]] constexpr bool valid() const noexcept {
        if (!is_finite(origin) || !is_finite(direction))
            return false;
        for (T v : direction.elements)
            if (v != T(0))
                return true;
        return false;
    }
};
template <floating T, std::size_t N = 3> struct segment {
    vec<T, N> a{}, b{};
};
template <floating T, std::size_t N = 3> struct aabb {
    vec<T, N> lower{std::numeric_limits<T>::max()}, upper{std::numeric_limits<T>::lowest()};
    [[nodiscard]] constexpr bool empty() const noexcept {
        for (std::size_t i = 0; i < N; ++i)
            if (!(lower[i] <= upper[i]))
                return true;
        return false;
    }
    [[nodiscard]] constexpr bool valid() const noexcept {
        return !empty() && is_finite(lower) && is_finite(upper);
    }
    constexpr void expand(const vec<T, N> &p) noexcept {
        lower = min(lower, p);
        upper = max(upper, p);
    }
    constexpr void expand(const aabb &b) noexcept {
        if (!b.empty()) {
            lower = min(lower, b.lower);
            upper = max(upper, b.upper);
        }
    }
    [[nodiscard]] constexpr vec<T, N> center() const noexcept {
        vec<T, N> r;
        for (std::size_t i = 0; i < N; ++i)
            r[i] = std::midpoint(lower[i], upper[i]);
        return r;
    }
    [[nodiscard]] constexpr vec<T, N> extent() const noexcept {
        return empty() ? vec<T, N>{} : upper - lower;
    }
    [[nodiscard]] constexpr bool contains(const vec<T, N> &p) const noexcept {
        for (std::size_t i = 0; i < N; ++i)
            if (!(p[i] >= lower[i] && p[i] <= upper[i]))
                return false;
        return true;
    }
};
template <floating T, std::size_t N = 3> struct sphere {
    vec<T, N> center{};
    T radius{};
    [[nodiscard]] constexpr bool valid() const noexcept {
        return is_finite(center) && is_finite(radius) && radius >= T(0);
    }
};
template <floating T> struct plane {
    vec<T, 3> normal{T(0), T(1), T(0)};
    T offset{}; // dot(normal,p)+offset=0; unit normal.
    [[nodiscard]] constexpr T signed_distance(const vec<T, 3> &p) const noexcept {
        return dot(normal, p) + offset;
    }
};
template <floating T> struct triangle {
    vec<T, 3> a{}, b{}, c{};
};
template <floating T> struct obb {
    vec<T, 3> center{};
    vec<T, 3> half_extent{};
    mat<T, 3, 3> axes = mat<T, 3, 3>::identity(); // Orthonormal columns.
    [[nodiscard]] constexpr bool valid(T tolerance = T(4) * epsilon<T>) const noexcept {
        return is_finite(center) && is_finite(half_extent) && is_finite(axes) &&
               half_extent[0] >= T(0) && half_extent[1] >= T(0) && half_extent[2] >= T(0) &&
               almost_equal(transpose(axes) * axes, mat<T, 3, 3>::identity(), tolerance, tolerance);
    }
};
template <floating T>
[[nodiscard]] inline std::optional<plane<T>> make_plane(const vec<T, 3> &normal,
                                                        T offset) noexcept {
    if (!is_finite(normal) || !is_finite(offset))
        return std::nullopt;
    const T s = std::max({std::abs(normal[0]), std::abs(normal[1]), std::abs(normal[2])});
    if (s == T(0))
        return std::nullopt;
    const T l = length(normal / s), d = (offset / s) / l;
    if (!is_finite(d))
        return std::nullopt;
    return plane<T>{(normal / s) / l, d};
}
template <floating T>
[[nodiscard]] inline std::optional<plane<T>>
plane_from_points(const vec<T, 3> &a, const vec<T, 3> &b, const vec<T, 3> &c) noexcept {
    if (!is_finite(a) || !is_finite(b) || !is_finite(c))
        return std::nullopt;
    auto ab = b - a, ac = c - a;
    if (!is_finite(ab) || !is_finite(ac)) {
        T scale{};
        for (std::size_t i = 0; i < 3; ++i)
            scale = std::max({scale, std::abs(a[i]), std::abs(b[i]), std::abs(c[i])});
        ab = b / scale - a / scale;
        ac = c / scale - a / scale;
    }
    const auto u = try_normalize(ab), v = try_normalize(ac);
    if (!u || !v)
        return std::nullopt;
    // Collinear points have no unique plane. The magnitude test is the same
    // relative degeneracy criterion used by barycentric(): with unit u and v it
    // compares sin^2(angle) against the tolerance. A plain exact-zero test on the
    // cross product is not sufficient because a fused multiply-subtract leaves a
    // non-zero rounding residual for two identical unit vectors.
    const auto n0 = cross(*u, *v);
    if (dot(n0, n0) <= epsilon<T>)
        return std::nullopt;
    const auto n = try_normalize(n0);
    if (!n)
        return std::nullopt;
    return make_plane(*n, -dot(*n, a));
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr bool overlaps(const aabb<T, N> &a, const aabb<T, N> &b) noexcept {
    if (a.empty() || b.empty())
        return false;
    for (std::size_t i = 0; i < N; ++i)
        if (a.upper[i] < b.lower[i] || b.upper[i] < a.lower[i])
            return false;
    return true;
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> closest_point(const aabb<T, N> &box,
                                                const vec<T, N> &p) noexcept {
    return clamp(p, box.lower, box.upper);
}
template <floating T, std::size_t N>
[[nodiscard]] inline vec<T, N> closest_point(const segment<T, N> &s, const vec<T, N> &p) noexcept {
    auto d = s.b - s.a, v = p - s.a;
    if ((!is_finite(d) || !is_finite(v)) && is_finite(s.a) && is_finite(s.b) && is_finite(p)) {
        T coordinate_scale{};
        for (std::size_t i = 0; i < N; ++i)
            coordinate_scale =
                std::max({coordinate_scale, std::abs(s.a[i]), std::abs(s.b[i]), std::abs(p[i])});
        d = s.b / coordinate_scale - s.a / coordinate_scale;
        v = p / coordinate_scale - s.a / coordinate_scale;
    }
    T scale{};
    for (T x : d.elements)
        scale = std::max(scale, std::abs(x));
    if (scale == T(0))
        return s.a;
    const auto dn = d / scale;
    const auto vn = v / scale;
    T parameter{};
    if (is_finite(vn)) {
        parameter = dot(vn, dn) / dot(dn, dn);
    } else {
        // A distant query must not create infinity*zero in perpendicular axes.
        T query_scale{};
        for (T x : v.elements)
            query_scale = std::max(query_scale, std::abs(x));
        parameter = multiply_divide(dot(v / query_scale, dn) / dot(dn, dn), query_scale, scale);
    }
    const T t = saturate(parameter);
    return lerp(s.a, s.b, t);
}
template <floating T>
[[nodiscard]] constexpr vec<T, 3> closest_point(const plane<T> &p,
                                                const vec<T, 3> &point) noexcept {
    return point - p.normal * p.signed_distance(point);
}
template <floating T, std::size_t N>
[[nodiscard]] inline bool overlaps(const sphere<T, N> &a, const sphere<T, N> &b) noexcept {
    if (!a.valid() || !b.valid())
        return false;
    const auto delta = a.center - b.center;
    const T radii = a.radius + b.radius;
    if (is_finite(delta) && is_finite(radii))
        return length(delta) <= radii;
    T scale = std::max(a.radius, b.radius);
    for (std::size_t i = 0; i < N; ++i)
        scale = std::max({scale, std::abs(a.center[i]), std::abs(b.center[i])});
    return distance(a.center / scale, b.center / scale) <= a.radius / scale + b.radius / scale;
}
template <floating T, std::size_t N>
[[nodiscard]] inline bool overlaps(const sphere<T, N> &s, const aabb<T, N> &b) noexcept {
    return s.valid() && b.valid() && distance(s.center, closest_point(b, s.center)) <= s.radius;
}
template <floating T> [[nodiscard]] constexpr T surface_area(const aabb<T, 3> &b) noexcept {
    const auto e = b.extent();
    return T(2) * (e[0] * e[1] + e[1] * e[2] + e[2] * e[0]);
}
template <floating T> [[nodiscard]] constexpr T volume(const aabb<T, 3> &b) noexcept {
    const auto e = b.extent();
    return e[0] * e[1] * e[2];
}
template <floating T> [[nodiscard]] inline T area(const triangle<T> &t) noexcept {
    return length(cross(t.b - t.a, t.c - t.a)) / T(2);
}
template <floating T>
[[nodiscard]] inline std::optional<vec<T, 3>>
barycentric(const triangle<T> &tri, const vec<T, 3> &p, T tolerance = epsilon<T>) noexcept {
    if (!is_finite(tri.a) || !is_finite(tri.b) || !is_finite(tri.c) || !is_finite(p) ||
        !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    auto a = tri.b - tri.a, b = tri.c - tri.a, v = p - tri.a;
    const T scale = std::max(length(a), length(b));
    if (scale == T(0) || !is_finite(scale))
        return std::nullopt;
    a /= scale;
    b /= scale;
    v /= scale;
    const auto normal = cross(a, b);
    const T denom = dot(normal, normal);
    if (denom <= tolerance * dot(a, a) * dot(b, b))
        return std::nullopt;
    const T y = dot(cross(v, b), normal) / denom, z = dot(cross(a, v), normal) / denom;
    const vec<T, 3> r{T(1) - y - z, y, z};
    if (!is_finite(r))
        return std::nullopt;
    return r;
}
template <floating T>
[[nodiscard]] inline vec<T, 3> closest_point(const triangle<T> &tri, const vec<T, 3> &p) noexcept {
    const auto bc = barycentric(tri, p, T(0));
    if (bc && (*bc)[0] >= T(0) && (*bc)[1] >= T(0) && (*bc)[2] >= T(0))
        return tri.a * (*bc)[0] + tri.b * (*bc)[1] + tri.c * (*bc)[2];
    const auto a = closest_point(segment<T>{tri.a, tri.b}, p),
               b = closest_point(segment<T>{tri.b, tri.c}, p),
               c = closest_point(segment<T>{tri.c, tri.a}, p);
    const T da = distance(a, p), db = distance(b, p), dc = distance(c, p);
    return da <= db && da <= dc ? a : (db <= dc ? b : c);
}
// Transform the bounds by interval arithmetic, including reflection and shear.
template <floating T>
[[nodiscard]] constexpr aabb<T, 3> transform_bounds(const affine3<T> &a,
                                                    const aabb<T, 3> &b) noexcept {
    if (b.empty())
        return {};
    aabb<T, 3> r{a.translation(), a.translation()};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const T x = a.matrix(i, j) * b.lower[j], y = a.matrix(i, j) * b.upper[j];
            r.lower[i] += std::min(x, y);
            r.upper[i] += std::max(x, y);
        }
    return r;
}
using ray3f = ray<float>;
using ray3d = ray<double>;
using aabb3f = aabb<float>;
using aabb3d = aabb<double>;
using spheref = sphere<float>;
using sphered = sphere<double>;

} // namespace chm
