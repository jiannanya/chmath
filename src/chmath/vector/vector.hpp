#pragma once

#include <array>
#include <cassert>
#include <chmath/core/scalar.hpp>
#include <optional>

namespace chm {

// Packed, contiguous, zero-initialized; no SIMD padding or pointer aliasing tricks.
template <scalar T, std::size_t N>
    requires(N > 0)
struct vec {
    std::array<T, N> elements{};
    using value_type = T;
    static constexpr std::size_t size = N;

    constexpr vec() noexcept = default;
    explicit constexpr vec(T value) noexcept { elements.fill(value); }
    template <class... U>
        requires(sizeof...(U) == N && N > 1 && (std::convertible_to<U, T> && ...))
    constexpr vec(U... values) noexcept : elements{static_cast<T>(values)...} {}
    explicit constexpr vec(const std::array<T, N> &values) noexcept : elements(values) {}
    template <scalar U> explicit constexpr vec(const vec<U, N> &other) noexcept {
        for (std::size_t i = 0; i < N; ++i)
            elements[i] = static_cast<T>(other[i]);
    }

    [[nodiscard]] constexpr T *data() noexcept { return elements.data(); }
    [[nodiscard]] constexpr const T *data() const noexcept { return elements.data(); }
    [[nodiscard]] constexpr T &operator[](std::size_t i) noexcept {
        assert(i < N);
        return elements[i];
    }
    [[nodiscard]] constexpr const T &operator[](std::size_t i) const noexcept {
        assert(i < N);
        return elements[i];
    }
    [[nodiscard]] constexpr T &x() noexcept { return elements[0]; }
    [[nodiscard]] constexpr const T &x() const noexcept { return elements[0]; }
    [[nodiscard]] constexpr T &y() noexcept
        requires(N >= 2)
    {
        return elements[1];
    }
    [[nodiscard]] constexpr const T &y() const noexcept
        requires(N >= 2)
    {
        return elements[1];
    }
    [[nodiscard]] constexpr T &z() noexcept
        requires(N >= 3)
    {
        return elements[2];
    }
    [[nodiscard]] constexpr const T &z() const noexcept
        requires(N >= 3)
    {
        return elements[2];
    }
    [[nodiscard]] constexpr T &w() noexcept
        requires(N >= 4)
    {
        return elements[3];
    }
    [[nodiscard]] constexpr const T &w() const noexcept
        requires(N >= 4)
    {
        return elements[3];
    }

    constexpr vec &operator+=(const vec &b) noexcept {
        for (std::size_t i = 0; i < N; ++i)
            elements[i] += b[i];
        return *this;
    }
    constexpr vec &operator-=(const vec &b) noexcept {
        for (std::size_t i = 0; i < N; ++i)
            elements[i] -= b[i];
        return *this;
    }
    constexpr vec &operator*=(T s) noexcept {
        for (auto &x_ : elements)
            x_ *= s;
        return *this;
    }
    constexpr vec &operator/=(T s) noexcept {
        for (auto &x_ : elements)
            x_ /= s;
        return *this;
    }
    [[nodiscard]] constexpr bool operator==(const vec &) const noexcept = default;
};

template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> operator+(vec<T, N> a, const vec<T, N> &b) noexcept {
    return a += b;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> operator-(vec<T, N> a, const vec<T, N> &b) noexcept {
    return a -= b;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> operator-(vec<T, N> a) noexcept {
    for (auto &x : a.elements)
        x = -x;
    return a;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> operator*(vec<T, N> a, T s) noexcept {
    return a *= s;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> operator*(T s, vec<T, N> a) noexcept {
    return a *= s;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> operator/(vec<T, N> a, T s) noexcept {
    return a /= s;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> hadamard(vec<T, N> a, const vec<T, N> &b) noexcept {
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b[i];
    return a;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr T dot(const vec<T, N> &a, const vec<T, N> &b) noexcept {
    T result{};
    for (std::size_t i = 0; i < N; ++i)
        result += a[i] * b[i];
    return result;
}
template <scalar T>
[[nodiscard]] constexpr vec<T, 3> cross(const vec<T, 3> &a, const vec<T, 3> &b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
template <scalar T>
[[nodiscard]] constexpr T cross(const vec<T, 2> &a, const vec<T, 2> &b) noexcept {
    return a[0] * b[1] - a[1] * b[0];
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr T length_squared(const vec<T, N> &a) noexcept {
    return dot(a, a);
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr bool is_finite(const vec<T, N> &a) noexcept {
    for (T x : a.elements)
        if (!is_finite(x))
            return false;
    return true;
}
// Scaling avoids overflow/underflow in the norm and normalization of finite vectors.
template <floating T, std::size_t N> [[nodiscard]] inline T length(const vec<T, N> &a) noexcept {
    const T squared = dot(a, a);
    if (squared >= std::numeric_limits<T>::min() && is_finite(squared))
        return std::sqrt(squared);
    T scale{};
    for (T x : a.elements) {
        if (x != x)
            return x;
        scale = std::max(scale, std::abs(x));
    }
    if (scale == T(0) || !is_finite(scale))
        return scale;
    return scale * std::sqrt(length_squared(a / scale));
}
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<vec<T, N>> try_normalize(const vec<T, N> &a) noexcept {
    // A finite normal sum of squares proves all components finite and avoids
    // component scaling for the overwhelmingly common representable range.
    const T squared = dot(a, a);
    if (squared >= std::numeric_limits<T>::min() && is_finite(squared))
        return a * (T(1) / std::sqrt(squared));
    T scale{};
    for (T x : a.elements) {
        if (!is_finite(x))
            return std::nullopt;
        scale = std::max(scale, std::abs(x));
    }
    if (scale == T(0))
        return std::nullopt;
    const auto scaled = a / scale;
    return scaled * (T(1) / std::sqrt(dot(scaled, scaled)));
}
template <floating T, std::size_t N>
[[nodiscard]] inline vec<T, N> normalize(const vec<T, N> &a) noexcept {
    return try_normalize(a).value_or(vec<T, N>{});
}
template <floating T, std::size_t N>
[[nodiscard]] inline T distance(const vec<T, N> &a, const vec<T, N> &b) noexcept {
    return length(a - b);
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr T distance_squared(const vec<T, N> &a, const vec<T, N> &b) noexcept {
    return length_squared(a - b);
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> lerp(const vec<T, N> &a, const vec<T, N> &b, T t) noexcept {
    vec<T, N> r;
    for (std::size_t i = 0; i < N; ++i)
        r[i] = chm::lerp(a[i], b[i], t);
    return r;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> min(const vec<T, N> &a, const vec<T, N> &b) noexcept {
    vec<T, N> r;
    for (std::size_t i = 0; i < N; ++i)
        r[i] = std::min(a[i], b[i]);
    return r;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> max(const vec<T, N> &a, const vec<T, N> &b) noexcept {
    vec<T, N> r;
    for (std::size_t i = 0; i < N; ++i)
        r[i] = std::max(a[i], b[i]);
    return r;
}
template <scalar T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> clamp(const vec<T, N> &a, const vec<T, N> &lo,
                                        const vec<T, N> &hi) noexcept {
    return min(max(a, lo), hi);
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr bool almost_equal(const vec<T, N> &a, const vec<T, N> &b,
                                          T rel = epsilon<T>, T abs = epsilon<T>) noexcept {
    for (std::size_t i = 0; i < N; ++i)
        if (!almost_equal(a[i], b[i], rel, abs))
            return false;
    return true;
}
// n must be a unit normal.
template <floating T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> reflect(const vec<T, N> &incident, const vec<T, N> &n) noexcept {
    return incident - T(2) * dot(incident, n) * n;
}
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<vec<T, N>> refract(const vec<T, N> &incident, const vec<T, N> &n,
                                                      T eta) noexcept {
    if (!is_finite(incident) || !is_finite(n) || !is_finite(eta) || eta <= T(0))
        return std::nullopt;
    // Work with the tangential component to avoid cancellation at the critical angle
    // and infinity*zero for finite, very large eta at normal incidence.
    const auto tangent = (incident - n * dot(n, incident)) * eta;
    const T k = T(1) - dot(tangent, tangent);
    if (!is_finite(k) || k < T(0))
        return std::nullopt;
    const auto result = tangent - std::sqrt(k) * n;
    if (!is_finite(result))
        return std::nullopt;
    return result;
}
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<vec<T, N>> project(const vec<T, N> &a,
                                                      const vec<T, N> &onto) noexcept {
    if (!is_finite(a))
        return std::nullopt;
    const auto n = try_normalize(onto);
    if (!n)
        return std::nullopt;
    const auto result = *n * dot(a, *n);
    if (!is_finite(result))
        return std::nullopt;
    return result;
}
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<T> angle_between(const vec<T, N> &a,
                                                    const vec<T, N> &b) noexcept {
    const auto na = try_normalize(a), nb = try_normalize(b);
    if (!na || !nb)
        return std::nullopt;
    // More accurate near parallel/antiparallel than acos(dot).
    return T(2) * std::atan2(length(*na - *nb), length(*na + *nb));
}

template <class T> using vec2 = vec<T, 2>;
template <class T> using vec3 = vec<T, 3>;
template <class T> using vec4 = vec<T, 4>;
using vec2f = vec2<float>;
using vec3f = vec3<float>;
using vec4f = vec4<float>;
using vec2d = vec2<double>;
using vec3d = vec3<double>;
using vec4d = vec4<double>;
using vec2i = vec2<int>;
using vec3i = vec3<int>;
using vec4i = vec4<int>;

} // namespace chm
