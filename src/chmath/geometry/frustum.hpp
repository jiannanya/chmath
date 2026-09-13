#pragma once

#include <chmath/geometry/primitives.hpp>
#include <cstdint>

namespace chm {

enum class containment { outside, intersecting, inside };
template <floating T> struct frustum {
    std::array<plane<T>, 6> planes{}; // left,right,bottom,top,lower-depth,upper-depth
    std::uint8_t active_mask{};       // Infinite projections may have an inactive depth plane.
};
template <floating T>
[[nodiscard]] inline std::optional<frustum<T>>
extract_frustum(const mat<T, 4, 4> &view_projection,
                depth_range range = depth_range::zero_to_one) noexcept {
    if (!is_finite(view_projection))
        return std::nullopt;
    const auto x = view_projection.row(0), y = view_projection.row(1), z = view_projection.row(2),
               w = view_projection.row(3);
    const std::array<vec<T, 4>, 6> coefficients{
        w + x, w - x, w + y, w - y, range == depth_range::zero_to_one ? z : w + z, w - z};
    frustum<T> result;
    for (std::size_t i = 0; i < 6; ++i) {
        const auto &c = coefficients[i];
        const vec<T, 3> n{c[0], c[1], c[2]};
        if (n == vec<T, 3>{}) {
            if (c[3] < T(0))
                return std::nullopt;
            continue;
        }
        const auto p = make_plane(n, c[3]);
        if (!p)
            return std::nullopt;
        result.planes[i] = *p;
        result.active_mask |= static_cast<std::uint8_t>(1U << i);
    }
    if (result.active_mask == 0)
        return std::nullopt;
    return result;
}
template <floating T>
[[nodiscard]] constexpr containment classify(const frustum<T> &f, const vec<T, 3> &p) noexcept {
    if (!is_finite(p))
        return containment::outside;
    for (std::size_t i = 0; i < 6; ++i)
        if ((f.active_mask & (1U << i)) != 0 && f.planes[i].signed_distance(p) < T(0))
            return containment::outside;
    return containment::inside;
}
template <floating T>
[[nodiscard]] constexpr containment classify(const frustum<T> &f, const sphere<T> &s) noexcept {
    if (!s.valid())
        return containment::outside;
    bool inside = true;
    for (std::size_t i = 0; i < 6; ++i)
        if ((f.active_mask & (1U << i)) != 0) {
            const T d = f.planes[i].signed_distance(s.center);
            if (d < -s.radius)
                return containment::outside;
            if (d < s.radius)
                inside = false;
        }
    return inside ? containment::inside : containment::intersecting;
}
template <floating T>
[[nodiscard]] constexpr containment classify(const frustum<T> &f, const aabb<T> &b) noexcept {
    if (!b.valid())
        return containment::outside;
    bool inside = true;
    for (std::size_t i = 0; i < 6; ++i)
        if ((f.active_mask & (1U << i)) != 0) {
            vec<T, 3> p, n;
            for (std::size_t j = 0; j < 3; ++j) {
                const bool positive = f.planes[i].normal[j] >= T(0);
                p[j] = positive ? b.upper[j] : b.lower[j];
                n[j] = positive ? b.lower[j] : b.upper[j];
            }
            if (f.planes[i].signed_distance(p) < T(0))
                return containment::outside;
            if (f.planes[i].signed_distance(n) < T(0))
                inside = false;
        }
    return inside ? containment::inside : containment::intersecting;
}

} // namespace chm
