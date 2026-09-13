#pragma once

#include <chmath/matrix/matrix.hpp>

namespace chm {

// Hamilton product, xyzw storage, identity by default. Rotation operations require unit
// quaternions.
template <floating T> struct quat {
    vec<T, 4> coefficients{T(0), T(0), T(0), T(1)};
    constexpr quat() noexcept = default;
    constexpr quat(T x, T y, T z, T w) noexcept : coefficients(x, y, z, w) {}
    constexpr quat(const vec<T, 3> &v, T w) noexcept : coefficients(v[0], v[1], v[2], w) {}
    explicit constexpr quat(const vec<T, 4> &v) noexcept : coefficients(v) {}
    [[nodiscard]] constexpr T &x() noexcept { return coefficients[0]; }
    [[nodiscard]] constexpr T x() const noexcept { return coefficients[0]; }
    [[nodiscard]] constexpr T &y() noexcept { return coefficients[1]; }
    [[nodiscard]] constexpr T y() const noexcept { return coefficients[1]; }
    [[nodiscard]] constexpr T &z() noexcept { return coefficients[2]; }
    [[nodiscard]] constexpr T z() const noexcept { return coefficients[2]; }
    [[nodiscard]] constexpr T &w() noexcept { return coefficients[3]; }
    [[nodiscard]] constexpr T w() const noexcept { return coefficients[3]; }
    [[nodiscard]] constexpr vec<T, 3> imaginary() const noexcept { return {x(), y(), z()}; }
    [[nodiscard]] constexpr T *data() noexcept { return coefficients.data(); }
    [[nodiscard]] constexpr const T *data() const noexcept { return coefficients.data(); }
    [[nodiscard]] constexpr bool operator==(const quat &) const noexcept = default;
};
template <floating T>
[[nodiscard]] constexpr quat<T> operator+(const quat<T> &a, const quat<T> &b) noexcept {
    return quat<T>(a.coefficients + b.coefficients);
}
template <floating T>
[[nodiscard]] constexpr quat<T> operator-(const quat<T> &a, const quat<T> &b) noexcept {
    return quat<T>(a.coefficients - b.coefficients);
}
template <floating T> [[nodiscard]] constexpr quat<T> operator-(const quat<T> &q) noexcept {
    return quat<T>(-q.coefficients);
}
template <floating T> [[nodiscard]] constexpr quat<T> operator*(const quat<T> &q, T s) noexcept {
    return quat<T>(q.coefficients * s);
}
template <floating T> [[nodiscard]] constexpr quat<T> operator*(T s, const quat<T> &q) noexcept {
    return q * s;
}
template <floating T> [[nodiscard]] constexpr quat<T> operator/(const quat<T> &q, T s) noexcept {
    return quat<T>(q.coefficients / s);
}
template <floating T>
[[nodiscard]] constexpr quat<T> operator*(const quat<T> &a, const quat<T> &b) noexcept {
    return {a.w() * b.imaginary() + b.w() * a.imaginary() + cross(a.imaginary(), b.imaginary()),
            a.w() * b.w() - dot(a.imaginary(), b.imaginary())};
}
template <floating T> [[nodiscard]] constexpr T dot(const quat<T> &a, const quat<T> &b) noexcept {
    return dot(a.coefficients, b.coefficients);
}
template <floating T> [[nodiscard]] constexpr quat<T> conjugate(const quat<T> &q) noexcept {
    return {-q.imaginary(), q.w()};
}
template <floating T> [[nodiscard]] constexpr bool is_finite(const quat<T> &q) noexcept {
    return is_finite(q.coefficients);
}
template <floating T>
[[nodiscard]] inline std::optional<quat<T>> try_normalize(const quat<T> &q) noexcept {
    const auto v = try_normalize(q.coefficients);
    if (!v)
        return std::nullopt;
    return quat<T>(*v);
}
template <floating T> [[nodiscard]] inline quat<T> normalize(const quat<T> &q) noexcept {
    return try_normalize(q).value_or(quat<T>{});
}
template <floating T>
[[nodiscard]] inline std::optional<quat<T>> inverse(const quat<T> &q) noexcept {
    if (!is_finite(q))
        return std::nullopt;
    T s{};
    for (T v : q.coefficients.elements)
        s = std::max(s, std::abs(v));
    if (s == T(0))
        return std::nullopt;
    const quat<T> p = q / s;
    const auto result = (conjugate(p) / dot(p, p)) / s;
    if (!is_finite(result))
        return std::nullopt;
    return result;
}
template <floating T>
[[nodiscard]] constexpr vec<T, 3> rotate(const quat<T> &q, const vec<T, 3> &v) noexcept {
    const auto t = T(2) * cross(q.imaginary(), v);
    return v + q.w() * t + cross(q.imaginary(), t);
}
template <floating T>
[[nodiscard]] inline std::optional<quat<T>> from_axis_angle(const vec<T, 3> &axis,
                                                            T angle) noexcept {
    const auto n = try_normalize(axis);
    if (!n || !is_finite(angle))
        return std::nullopt;
    return quat<T>(*n * std::sin(angle / T(2)), std::cos(angle / T(2)));
}
template <floating T> struct axis_angle {
    vec<T, 3> axis{T(1), T(0), T(0)};
    T angle{};
};
template <floating T> [[nodiscard]] inline axis_angle<T> to_axis_angle(quat<T> q) noexcept {
    q = normalize(q);
    if (q.w() < T(0))
        q = -q;
    const T s = length(q.imaginary());
    if (s == T(0))
        return {};
    return {q.imaginary() / s, T(2) * std::atan2(s, q.w())};
}
// Fixed axes X, then Y, then Z: q = qz*qy*qx (radians).
template <floating T>
[[nodiscard]] inline quat<T> from_euler_xyz(const vec<T, 3> &angles) noexcept {
    const T cx = std::cos(angles[0] / T(2)), sx = std::sin(angles[0] / T(2));
    const T cy = std::cos(angles[1] / T(2)), sy = std::sin(angles[1] / T(2));
    const T cz = std::cos(angles[2] / T(2)), sz = std::sin(angles[2] / T(2));
    return {sx * cy * cz - cx * sy * sz, cx * sy * cz + sx * cy * sz, cx * cy * sz - sx * sy * cz,
            cx * cy * cz + sx * sy * sz};
}
template <floating T> [[nodiscard]] constexpr mat<T, 3, 3> to_matrix(const quat<T> &q) noexcept {
    const T x = q.x(), y = q.y(), z = q.z(), w = q.w();
    return mat<T, 3, 3>::from_rows(
        {vec<T, 3>{T(1) - T(2) * (y * y + z * z), T(2) * (x * y - z * w), T(2) * (x * z + y * w)},
         vec<T, 3>{T(2) * (x * y + z * w), T(1) - T(2) * (x * x + z * z), T(2) * (y * z - x * w)},
         vec<T, 3>{T(2) * (x * z - y * w), T(2) * (y * z + x * w), T(1) - T(2) * (x * x + y * y)}});
}
template <floating T>
[[nodiscard]] inline std::optional<quat<T>> from_matrix(const mat<T, 3, 3> &m,
                                                        T tolerance = T(4) * epsilon<T>) noexcept {
    if (!is_finite(m) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    if (!almost_equal(transpose(m) * m, mat<T, 3, 3>::identity(), tolerance, tolerance) ||
        !(dot(m.column(0), cross(m.column(1), m.column(2))) > T(0)))
        return std::nullopt;
    quat<T> q;
    const T tr = trace(m);
    if (tr > T(0)) {
        const T s = T(2) * std::sqrt(tr + T(1));
        q = {(m(2, 1) - m(1, 2)) / s, (m(0, 2) - m(2, 0)) / s, (m(1, 0) - m(0, 1)) / s, s / T(4)};
    } else {
        std::size_t i = 0;
        if (m(1, 1) > m(0, 0))
            i = 1;
        if (m(2, 2) > m(i, i))
            i = 2;
        const std::size_t j = (i + 1) % 3, k = (i + 2) % 3;
        const T s = T(2) * std::sqrt(std::max(T(0), T(1) + m(i, i) - m(j, j) - m(k, k)));
        q.coefficients[i] = s / T(4);
        q.coefficients[j] = (m(j, i) + m(i, j)) / s;
        q.coefficients[k] = (m(k, i) + m(i, k)) / s;
        q.w() = (m(k, j) - m(j, k)) / s;
    }
    return try_normalize(q);
}
template <floating T> [[nodiscard]] inline quat<T> nlerp(quat<T> a, quat<T> b, T t) noexcept {
    a = normalize(a);
    b = normalize(b);
    if (dot(a, b) < T(0))
        b = -b;
    return normalize(a * (T(1) - t) + b * t);
}
template <floating T> [[nodiscard]] inline quat<T> slerp(quat<T> a, quat<T> b, T t) noexcept {
    a = normalize(a);
    b = normalize(b);
    T d = dot(a, b);
    if (d < T(0)) {
        b = -b;
        d = -d;
    }
    d = clamp(d, T(0), T(1));
    if (d > T(0.9995))
        return normalize(a * (T(1) - t) + b * t);
    const T angle = std::acos(d);
    // The common positive 1/sin(angle) factor cancels on normalization.
    return normalize(a * std::sin((T(1) - t) * angle) + b * std::sin(t * angle));
}
template <floating T>
[[nodiscard]] inline std::optional<quat<T>> rotation_between(const vec<T, 3> &from,
                                                             const vec<T, 3> &to) noexcept {
    const auto a = try_normalize(from), b = try_normalize(to);
    if (!a || !b)
        return std::nullopt;
    const T d = clamp(dot(*a, *b), T(-1), T(1));
    const auto c = cross(*a, *b);
    if (d < T(0) && length_squared(c) <= square(epsilon<T>)) {
        std::size_t i = 0;
        if (std::abs((*a)[1]) < std::abs((*a)[i]))
            i = 1;
        if (std::abs((*a)[2]) < std::abs((*a)[i]))
            i = 2;
        vec<T, 3> basis;
        basis[i] = T(1);
        return quat<T>(normalize(cross(*a, basis)), T(0));
    }
    // atan2 remains accurate near 180 degrees where 1+dot loses precision.
    const T s = length(c);
    if (s == T(0))
        return quat<T>{};
    const T half_angle = std::atan2(s, d) / T(2);
    return quat<T>((c / s) * std::sin(half_angle), std::cos(half_angle));
}
using quatf = quat<float>;
using quatd = quat<double>;

} // namespace chm
