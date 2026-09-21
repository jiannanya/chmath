#pragma once

#include <chmath/rotation/quaternion.hpp>

namespace chm {

enum class handedness { right, left };
enum class depth_range { minus_one_to_one, zero_to_one };
enum class depth_direction { forward, reverse };

// Affine 3x4: twelve scalars, implicit bottom row [0,0,0,1].
template <floating T> struct affine3 {
    mat<T, 3, 4> matrix{T(1)};
    constexpr affine3() noexcept = default;
    explicit constexpr affine3(const mat<T, 3, 4> &m) noexcept : matrix(m) {}
    constexpr affine3(const mat<T, 3, 3> &linear, const vec<T, 3> &translation) noexcept {
        for (std::size_t j = 0; j < 3; ++j)
            matrix.set_column(j, linear.column(j));
        matrix.set_column(3, translation);
    }
    [[nodiscard]] constexpr mat<T, 3, 3> linear() const noexcept {
        mat<T, 3, 3> r;
        for (std::size_t j = 0; j < 3; ++j)
            r.set_column(j, matrix.column(j));
        return r;
    }
    [[nodiscard]] constexpr vec<T, 3> translation() const noexcept { return matrix.column(3); }
    [[nodiscard]] constexpr bool operator==(const affine3 &) const noexcept = default;
};
template <floating T>
[[nodiscard]] constexpr vec<T, 3> transform_point(const affine3<T> &a,
                                                  const vec<T, 3> &p) noexcept {
    const auto &m = a.matrix;
    return {((m(0, 0) * p[0] + m(0, 1) * p[1]) + m(0, 2) * p[2]) + m(0, 3),
            ((m(1, 0) * p[0] + m(1, 1) * p[1]) + m(1, 2) * p[2]) + m(1, 3),
            ((m(2, 0) * p[0] + m(2, 1) * p[1]) + m(2, 2) * p[2]) + m(2, 3)};
}
template <floating T>
[[nodiscard]] constexpr vec<T, 3> transform_vector(const affine3<T> &a,
                                                   const vec<T, 3> &p) noexcept {
    const auto &m = a.matrix;
    return {(m(0, 0) * p[0] + m(0, 1) * p[1]) + m(0, 2) * p[2],
            (m(1, 0) * p[0] + m(1, 1) * p[1]) + m(1, 2) * p[2],
            (m(2, 0) * p[0] + m(2, 1) * p[1]) + m(2, 2) * p[2]};
}
template <floating T>
[[nodiscard]] constexpr affine3<T> operator*(const affine3<T> &a, const affine3<T> &b) noexcept {
    // Compose directly into the 3x4 result instead of materializing a.linear()
    // and then multiplying two temporary 3x3 matrices. The expression grouping
    // matches the generic small-matrix product and transform_point, so results
    // are unchanged.
    const auto &am = a.matrix;
    affine3<T> r;
    for (std::size_t j = 0; j < 3; ++j) {
        const auto column = b.matrix.column(j);
        for (std::size_t i = 0; i < 3; ++i)
            r.matrix(i, j) = (am(i, 0) * column[0] + am(i, 1) * column[1]) + am(i, 2) * column[2];
    }
    const auto translation = b.matrix.column(3);
    for (std::size_t i = 0; i < 3; ++i)
        r.matrix(i, 3) = ((am(i, 0) * translation[0] + am(i, 1) * translation[1]) +
                          am(i, 2) * translation[2]) +
                         am(i, 3);
    return r;
}
template <floating T> [[nodiscard]] constexpr mat<T, 4, 4> to_matrix(const affine3<T> &a) noexcept {
    auto r = mat<T, 4, 4>::identity();
    for (std::size_t j = 0; j < 4; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            r(i, j) = a.matrix(i, j);
    return r;
}
template <floating T>
[[nodiscard]] constexpr std::optional<affine3<T>> to_affine(const mat<T, 4, 4> &a,
                                                            T tolerance = epsilon<T>) noexcept {
    if (!is_finite(a) || !is_finite(tolerance) || tolerance < T(0) ||
        !almost_equal(a.row(3), vec<T, 4>{T(0), T(0), T(0), T(1)}, T(0), tolerance))
        return std::nullopt;
    affine3<T> r;
    for (std::size_t j = 0; j < 4; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            r.matrix(i, j) = a(i, j);
    return r;
}
template <floating T>
[[nodiscard]] constexpr std::optional<affine3<T>> inverse(const affine3<T> &a,
                                                          T tolerance = epsilon<T>) noexcept {
    if (!is_finite(a.matrix))
        return std::nullopt;
    const auto inv = inverse(a.linear(), tolerance);
    if (!inv)
        return std::nullopt;
    const affine3<T> r{*inv, -(*inv * a.translation())};
    if (!is_finite(r.matrix))
        return std::nullopt;
    return r;
}
template <floating T>
[[nodiscard]] constexpr std::optional<mat<T, 3, 3>>
normal_matrix(const affine3<T> &a, T tolerance = epsilon<T>) noexcept {
    const auto inv = inverse(a.linear(), tolerance);
    if (!inv)
        return std::nullopt;
    return transpose(*inv);
}
template <floating T>
[[nodiscard]] inline std::optional<vec<T, 3>> transform_normal(const affine3<T> &a,
                                                               const vec<T, 3> &n) noexcept {
    const auto nm = normal_matrix(a);
    if (!nm)
        return std::nullopt;
    return try_normalize(*nm * n);
}
template <floating T> [[nodiscard]] constexpr affine3<T> translation(const vec<T, 3> &v) noexcept {
    return {mat<T, 3, 3>::identity(), v};
}
template <floating T> [[nodiscard]] constexpr affine3<T> scaling(const vec<T, 3> &v) noexcept {
    mat<T, 3, 3> m;
    for (std::size_t i = 0; i < 3; ++i)
        m(i, i) = v[i];
    return {m, {}};
}
template <floating T> [[nodiscard]] constexpr affine3<T> rotation(const quat<T> &q) noexcept {
    return {to_matrix(q), {}};
}
template <floating T> struct trs {
    vec<T, 3> position{};
    quat<T> orientation{};
    vec<T, 3> scale{T(1)};
};
template <floating T> [[nodiscard]] constexpr affine3<T> compose(const trs<T> &t) noexcept {
    auto m = to_matrix(t.orientation);
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            m(i, j) *= t.scale[j];
    return {m, t.position};
}
// Rejects shear and zero scale. Reflection is represented by a negative X scale.
template <floating T>
[[nodiscard]] inline std::optional<trs<T>> decompose(const affine3<T> &a,
                                                     T tolerance = T(4) * epsilon<T>) noexcept {
    if (!is_finite(a.matrix))
        return std::nullopt;
    trs<T> r;
    r.position = a.translation();
    auto m = a.linear();
    for (std::size_t j = 0; j < 3; ++j) {
        r.scale[j] = length(m.column(j));
        if (r.scale[j] == T(0) || !is_finite(r.scale[j]))
            return std::nullopt;
        m.set_column(j, m.column(j) / r.scale[j]);
    }
    if (dot(m.column(0), cross(m.column(1), m.column(2))) < T(0)) {
        r.scale[0] = -r.scale[0];
        m.set_column(0, -m.column(0));
    }
    const auto q = from_matrix(m, tolerance);
    if (!q)
        return std::nullopt;
    r.orientation = *q;
    return r;
}
template <floating T>
[[nodiscard]] constexpr std::optional<vec<T, 3>> transform_point(const mat<T, 4, 4> &m,
                                                                 const vec<T, 3> &p) noexcept {
    const auto h = m * vec<T, 4>{p[0], p[1], p[2], T(1)};
    if (!is_finite(h) || h[3] == T(0))
        return std::nullopt;
    const vec<T, 3> r{h[0] / h[3], h[1] / h[3], h[2] / h[3]};
    if (!is_finite(r))
        return std::nullopt;
    return r;
}
template <floating T>
[[nodiscard]] constexpr vec<T, 3> transform_vector(const mat<T, 4, 4> &m,
                                                   const vec<T, 3> &p) noexcept {
    const auto h = m * vec<T, 4>{p[0], p[1], p[2], T(0)};
    return {h[0], h[1], h[2]};
}
template <floating T>
[[nodiscard]] inline std::optional<affine3<T>>
look_at(const vec<T, 3> &eye, const vec<T, 3> &target, const vec<T, 3> &up,
        handedness hand = handedness::right) noexcept {
    if (!is_finite(eye) || !is_finite(target))
        return std::nullopt;
    auto direction = hand == handedness::right ? eye - target : target - eye;
    if (!is_finite(direction)) {
        T scale{};
        for (std::size_t i = 0; i < 3; ++i)
            scale = std::max({scale, std::abs(eye[i]), std::abs(target[i])});
        direction =
            hand == handedness::right ? eye / scale - target / scale : target / scale - eye / scale;
    }
    const auto z = try_normalize(direction);
    if (!z)
        return std::nullopt;
    const auto x = try_normalize(cross(up, *z));
    if (!x)
        return std::nullopt;
    const auto y = cross(*z, *x);
    const auto m = mat<T, 3, 3>::from_rows({*x, y, *z});
    const affine3<T> r{m, -(m * eye)};
    if (!is_finite(r.matrix))
        return std::nullopt;
    return r;
}
template <floating T>
[[nodiscard]] inline std::optional<mat<T, 4, 4>>
perspective(T fovy, T aspect, T near_z, T far_z, handedness hand = handedness::right,
            depth_range range = depth_range::zero_to_one,
            depth_direction direction = depth_direction::forward) noexcept {
    if (!is_finite(fovy) || !is_finite(aspect) || !is_finite(near_z) ||
        !(fovy > T(0) && fovy < pi<T>) || aspect <= T(0) || near_z <= T(0) || !(far_z > near_z))
        return std::nullopt;
    T dn = range == depth_range::zero_to_one ? T(0) : T(-1), df = T(1);
    if (direction == depth_direction::reverse)
        std::swap(dn, df);
    const T s = hand == handedness::right ? T(-1) : T(1), f = T(1) / std::tan(fovy / T(2));
    mat<T, 4, 4> m;
    m(0, 0) = f / aspect;
    m(1, 1) = f;
    m(3, 2) = s;
    if (far_z == std::numeric_limits<T>::infinity()) {
        m(2, 2) = s * df;
        m(2, 3) = near_z * (dn - df);
    } else {
        const T ratio = near_z / far_z;
        m(2, 2) = s * (df - ratio * dn) / (T(1) - ratio);
        m(2, 3) = near_z * (dn - df) / (T(1) - ratio);
    }
    if (!is_finite(m))
        return std::nullopt;
    return m;
}
template <floating T>
[[nodiscard]] inline std::optional<mat<T, 4, 4>>
orthographic(T left, T right, T bottom, T top, T near_z, T far_z,
             handedness hand = handedness::right, depth_range range = depth_range::zero_to_one,
             depth_direction direction = depth_direction::forward) noexcept {
    if (!is_finite(left) || !is_finite(right) || !is_finite(bottom) || !is_finite(top) ||
        !is_finite(near_z) || !is_finite(far_z) ||
        !(right > left && top > bottom && far_z > near_z))
        return std::nullopt;
    T dn = range == depth_range::zero_to_one ? T(0) : T(-1), df = T(1);
    if (direction == depth_direction::reverse)
        std::swap(dn, df);
    const T s = hand == handedness::right ? T(-1) : T(1);
    auto m = mat<T, 4, 4>::identity();
    m(0, 0) = T(2) / (right - left);
    m(1, 1) = T(2) / (top - bottom);
    m(2, 2) = s * (df - dn) / (far_z - near_z);
    m(0, 3) = -(right + left) / (right - left);
    m(1, 3) = -(top + bottom) / (top - bottom);
    m(2, 3) = (far_z * dn - near_z * df) / (far_z - near_z);
    if (!is_finite(m))
        return std::nullopt;
    return m;
}
template <floating T> struct viewport {
    T x{}, y{}, width{}, height{}, min_depth{}, max_depth = T(1);
    [[nodiscard]] constexpr bool valid() const noexcept {
        return is_finite(x) && is_finite(y) && is_finite(width) && is_finite(height) &&
               is_finite(min_depth) && is_finite(max_depth) && width > T(0) && height > T(0) &&
               max_depth > min_depth;
    }
};
// Viewport Y increases upward. For top-left window systems, flip Y in the caller.
template <floating T>
[[nodiscard]] constexpr std::optional<vec<T, 3>>
project(const vec<T, 3> &p, const mat<T, 4, 4> &mvp, const viewport<T> &vp,
        depth_range range = depth_range::zero_to_one) noexcept {
    if (!vp.valid())
        return std::nullopt;
    const auto ndc = transform_point(mvp, p);
    if (!ndc)
        return std::nullopt;
    const T z = range == depth_range::zero_to_one ? (*ndc)[2] : ((*ndc)[2] + T(1)) / T(2);
    const vec<T, 3> r{vp.x + ((*ndc)[0] + T(1)) * vp.width / T(2),
                      vp.y + ((*ndc)[1] + T(1)) * vp.height / T(2),
                      vp.min_depth + z * (vp.max_depth - vp.min_depth)};
    if (!is_finite(r))
        return std::nullopt;
    return r;
}
// Supply inverse(mvp) once for repeated unprojections.
template <floating T>
[[nodiscard]] constexpr std::optional<vec<T, 3>>
unproject(const vec<T, 3> &p, const mat<T, 4, 4> &inverse_mvp, const viewport<T> &vp,
          depth_range range = depth_range::zero_to_one) noexcept {
    // Keep divisor checks visible here for constant invalid viewports as well.
    if (vp.width == T(0) || vp.height == T(0) || vp.max_depth == vp.min_depth || !vp.valid())
        return std::nullopt;
    T z = (p[2] - vp.min_depth) / (vp.max_depth - vp.min_depth);
    if (range == depth_range::minus_one_to_one)
        z = T(2) * z - T(1);
    return transform_point(inverse_mvp, vec<T, 3>{T(2) * (p[0] - vp.x) / vp.width - T(1),
                                                  T(2) * (p[1] - vp.y) / vp.height - T(1), z});
}
using affine3f = affine3<float>;
using affine3d = affine3<double>;

} // namespace chm
