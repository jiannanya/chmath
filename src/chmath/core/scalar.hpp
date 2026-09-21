#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numbers>
#include <type_traits>

namespace chm {

template <class T>
concept scalar = std::is_arithmetic_v<T> && !std::same_as<T, bool>;

template <class T>
concept floating = std::floating_point<T>;

template <floating T> inline constexpr T pi = std::numbers::pi_v<T>;
template <floating T> inline constexpr T tau = T(2) * pi<T>;
template <floating T> inline constexpr T epsilon = T(32) * std::numeric_limits<T>::epsilon();

template <scalar T> [[nodiscard]] constexpr T square(T x) noexcept {
    return x * x;
}
template <scalar T> [[nodiscard]] constexpr T clamp(T x, T lo, T hi) noexcept {
    return std::min(std::max(x, lo), hi);
}
template <scalar T> [[nodiscard]] constexpr T saturate(T x) noexcept {
    return clamp(x, T(0), T(1));
}
template <floating T> [[nodiscard]] constexpr T radians(T degrees) noexcept {
    return degrees * (pi<T> / T(180));
}
template <floating T> [[nodiscard]] constexpr T degrees(T radians_) noexcept {
    return radians_ * (T(180) / pi<T>);
}
template <floating T> [[nodiscard]] constexpr T lerp(T a, T b, T t) noexcept {
    return std::lerp(a, b, t);
}
template <floating T> [[nodiscard]] constexpr T smoothstep(T lo, T hi, T x) noexcept {
    if (lo == hi)
        return x < lo ? T(0) : T(1);
    const T width = hi - lo;
    const T t =
        saturate(width <= std::numeric_limits<T>::max() && width >= std::numeric_limits<T>::lowest()
                     ? (x - lo) / width
                     : (x / T(2) - lo / T(2)) / (hi / T(2) - lo / T(2)));
    return t * t * (T(3) - T(2) * t);
}
// Three comparisons classify NaN and both infinities as non-finite. This form is
// kept deliberately: on the measured targets it vectorizes and schedules better
// than a library classification call in the validation loops that dominate the
// geometry and factorization paths.
template <floating T> [[nodiscard]] constexpr bool is_finite(T x) noexcept {
    return x == x && x <= std::numeric_limits<T>::max() && x >= std::numeric_limits<T>::lowest();
}
namespace detail {
// Interpolation used by internal hot loops (de Casteljau, de Boor). For finite
// operands whose difference is representable this is exactly a + t * (b - a),
// which is also the form std::lerp reduces to in that range; the guard keeps the
// public function's overflow protection when the difference itself overflows,
// and the t == 1 shortcut keeps the exact endpoint that std::lerp guarantees, so
// curve and patch corners still reproduce their control points bit for bit.
// Avoiding the remaining branches of std::lerp matters because this runs
// O(Degree^2) times per evaluation on very short vectors.
template <floating T> [[nodiscard]] constexpr T lerp_finite(T a, T b, T t) noexcept {
    if (t == T(1))
        return b;
    const T difference = b - a;
    return is_finite(difference) ? a + t * difference : std::lerp(a, b, t);
}
} // namespace detail
// Relative AND absolute tolerances are caller-selectable; infinities only equal themselves.
namespace detail {
// Body of almost_equal without the tolerance validation, so aggregate overloads
// can validate the tolerances once instead of once per component.
template <floating T>
[[nodiscard]] constexpr bool almost_equal_value(T a, T b, T rel, T abs) noexcept {
    if (a == b)
        return true;
    if (!is_finite(a) || !is_finite(b))
        return false;
    const T difference = std::abs(a - b), scale = std::max(std::abs(a), std::abs(b));
    const T relative_limit = rel * scale;
    if (is_finite(difference) && is_finite(relative_limit))
        return difference <= std::max(abs, relative_limit);
    return difference <= abs || std::abs(a / scale - b / scale) <= rel;
}
} // namespace detail
template <floating T> [[nodiscard]] constexpr bool tolerance_valid(T rel, T abs) noexcept {
    return is_finite(rel) && is_finite(abs) && rel >= T(0) && abs >= T(0);
}
template <floating T>
[[nodiscard]] constexpr bool almost_equal(T a, T b, T rel = epsilon<T>,
                                          T abs = epsilon<T>) noexcept {
    if (!tolerance_valid(rel, abs))
        return false;
    return detail::almost_equal_value(a, b, rel, abs);
}
template <floating T> [[nodiscard]] inline T wrap_angle(T x) noexcept {
    T r = std::remainder(x, tau<T>);
    return r == pi<T> ? -pi<T> : r;
}

// Compute a*b/c without avoidable overflow or underflow in the intermediate product.
template <floating T> [[nodiscard]] inline T multiply_divide(T a, T b, T c) noexcept {
    const T product = a * b;
    if (!is_finite(a) || !is_finite(b) || !is_finite(c) || c == T(0) || a == T(0) || b == T(0) ||
        (is_finite(product) && std::abs(product) >= std::numeric_limits<T>::min()))
        return product / c;
    int ea{}, eb{}, ec{};
    const T ma = std::frexp(a, &ea), mb = std::frexp(b, &eb), mc = std::frexp(c, &ec);
    return std::ldexp((ma * mb) / mc, ea + eb - ec);
}

inline constexpr int version_major = 2;
inline constexpr int version_minor = 0;
inline constexpr int version_patch = 0;

} // namespace chm
