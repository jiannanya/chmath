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
template <floating T> [[nodiscard]] constexpr bool is_finite(T x) noexcept {
    return x == x && x <= std::numeric_limits<T>::max() && x >= std::numeric_limits<T>::lowest();
}
// Relative AND absolute tolerances are caller-selectable; infinities only equal themselves.
template <floating T>
[[nodiscard]] constexpr bool almost_equal(T a, T b, T rel = epsilon<T>,
                                          T abs = epsilon<T>) noexcept {
    if (!is_finite(rel) || !is_finite(abs) || rel < T(0) || abs < T(0))
        return false;
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
