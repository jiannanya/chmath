#pragma once

#include <chmath/matrix/matrix.hpp>

namespace chm {

// Symmetric positive-definite factorization A=L*transpose(L).
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<mat<T, N, N>> cholesky(const mat<T, N, N> &a,
                                                          T tolerance = epsilon<T>) noexcept {
    if (!is_finite(a) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    T scale{};
    for (T v : a.elements)
        scale = std::max(scale, std::abs(v));
    if (scale == T(0))
        return std::nullopt;
    // Normalize once. transpose(a/scale) is bit-identical to transpose(a)/scale
    // because both divide the same element by the same scale, so sharing the
    // scaled copy removes two full-matrix divisions and two stack temporaries.
    // Symmetry is then checked element-wise, which also avoids materializing a
    // transpose.
    const auto scaled = a / scale;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < i; ++j)
            if (!detail::almost_equal_value(scaled(i, j), scaled(j, i), tolerance, tolerance))
                return std::nullopt;
    mat<T, N, N> l;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            T sum = scaled(i, j);
            for (std::size_t k = 0; k < j; ++k)
                sum -= l(i, k) * l(j, k);
            if (i == j) {
                if (sum <= tolerance)
                    return std::nullopt;
                l(i, j) = std::sqrt(sum);
            } else
                l(i, j) = sum / l(j, j);
        }
    l *= std::sqrt(scale);
    if (!is_finite(l))
        return std::nullopt;
    return l;
}
template <floating T, std::size_t R, std::size_t C>
    requires(R >= C)
struct qr_factorization {
    mat<T, R, R> q = mat<T, R, R>::identity();
    mat<T, R, C> r{};
};
namespace detail {
// Apply each reflector to R and let the caller accumulate Q or transform a RHS.
// The latter avoids the R-by-R workspace and O(R*R*C) work of forming full Q.
template <floating T, std::size_t R, std::size_t C, class Apply>
inline void householder_reduce(mat<T, R, C> &r, Apply apply) noexcept {
    for (std::size_t k = 0; k < C; ++k) {
        vec<T, R> v;
        T scale{};
        for (std::size_t i = k; i < R; ++i)
            scale = std::max(scale, std::abs(r(i, k)));
        if (scale == T(0))
            continue;
        for (std::size_t i = k; i < R; ++i)
            v[i] = r(i, k) / scale;
        const T norm = length(v);
        v[k] += std::copysign(norm, v[k]);
        v = normalize(v);
        for (std::size_t j = k; j < C; ++j) {
            T d{};
            for (std::size_t i = k; i < R; ++i)
                d += v[i] * r(i, j);
            for (std::size_t i = k; i < R; ++i)
                r(i, j) -= T(2) * v[i] * d;
        }
        apply(k, v);
        for (std::size_t i = k + 1; i < R; ++i)
            r(i, k) = T(0);
    }
}
} // namespace detail

// Householder QR. Full Q is appropriate for small fixed matrices; no heap workspace.
template <floating T, std::size_t R, std::size_t C>
    requires(R >= C)
[[nodiscard]] inline std::optional<qr_factorization<T, R, C>>
factor_qr(const mat<T, R, C> &a) noexcept {
    if (!is_finite(a))
        return std::nullopt;
    qr_factorization<T, R, C> result;
    result.r = a;
    detail::householder_reduce(result.r, [&](std::size_t k, const vec<T, R> &v) {
        if constexpr (R <= 4) {
            for (std::size_t i = 0; i < R; ++i) {
                T d{};
                for (std::size_t j = k; j < R; ++j)
                    d += result.q(i, j) * v[j];
                for (std::size_t j = k; j < R; ++j)
                    result.q(i, j) -= T(2) * d * v[j];
            }
        } else {
            // Traverse column-major Q contiguously. Each row retains the same
            // accumulation order, while independent rows can be vectorized.
            vec<T, R> d;
            for (std::size_t j = k; j < R; ++j)
                for (std::size_t i = 0; i < R; ++i)
                    d[i] += result.q(i, j) * v[j];
            for (std::size_t j = k; j < R; ++j)
                for (std::size_t i = 0; i < R; ++i)
                    result.q(i, j) -= T(2) * d[i] * v[j];
        }
    });
    if (!is_finite(result.q) || !is_finite(result.r))
        return std::nullopt;
    return result;
}
template <floating T, std::size_t R, std::size_t C>
    requires(R >= C)
[[nodiscard]] inline std::optional<vec<T, C>>
least_squares(const mat<T, R, C> &a, const vec<T, R> &b, T tolerance = epsilon<T>) noexcept {
    if (!is_finite(a) || !is_finite(b) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    auto r = a;
    auto rhs = b;
    detail::householder_reduce(r, [&](std::size_t k, const vec<T, R> &v) {
        T d{};
        for (std::size_t i = k; i < R; ++i)
            d += v[i] * rhs[i];
        if (!is_finite(d)) {
            // A projection can overflow even when the reflected RHS is finite.
            T rhs_scale{};
            for (std::size_t i = k; i < R; ++i)
                rhs_scale = std::max(rhs_scale, std::abs(rhs[i]));
            d = T(0);
            for (std::size_t i = k; i < R; ++i)
                d += v[i] * (rhs[i] / rhs_scale);
            for (std::size_t i = k; i < R; ++i)
                rhs[i] = (rhs[i] / rhs_scale - T(2) * v[i] * d) * rhs_scale;
        } else {
            for (std::size_t i = k; i < R; ++i) {
                const T value = rhs[i] - T(2) * v[i] * d;
                rhs[i] = is_finite(value) ? value : (rhs[i] / T(2) - v[i] * d) * T(2);
            }
        }
    });
    if (!is_finite(r))
        return std::nullopt;
    vec<T, C> x;
    T scale{};
    for (T v : a.elements)
        scale = std::max(scale, std::abs(v));
    if (scale == T(0))
        return std::nullopt;
    for (std::size_t i = C; i-- > 0;) {
        if (std::abs(r(i, i)) / scale <= tolerance)
            return std::nullopt;
        x[i] = rhs[i];
        for (std::size_t j = i + 1; j < C; ++j)
            x[i] -= r(i, j) * x[j];
        x[i] /= r(i, i);
    }
    if (!is_finite(x))
        return std::nullopt;
    return x;
}
template <floating T, std::size_t N> struct symmetric_eigen_result {
    vec<T, N> values{}; // Ascending eigenvalues; corresponding eigenvectors in columns.
    mat<T, N, N> vectors = mat<T, N, N>::identity();
    std::size_t sweeps{};
    bool converged = false;
};
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<symmetric_eigen_result<T, N>>
symmetric_eigen(const mat<T, N, N> &input, T tolerance = epsilon<T>,
                std::size_t max_sweeps = 32) noexcept {
    if (!is_finite(input) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    T scale{};
    for (T v : input.elements)
        scale = std::max(scale, std::abs(v));
    symmetric_eigen_result<T, N> result;
    if (scale == T(0)) {
        result.converged = true;
        return result;
    }
    auto a = input / scale;
    // Symmetry validation and symmetrization in one pass: element (i,j) of
    // (a + transpose(a))/2 is (a(i,j) + a(j,i))/2, and checking the upper
    // triangle against its mirror is exactly the almost_equal(a, transpose(a))
    // test the general form would perform. This removes two full-matrix
    // temporaries and passes.
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j) {
            if (!detail::almost_equal_value(a(i, j), a(j, i), tolerance, tolerance))
                return std::nullopt;
            const T average = (a(i, j) + a(j, i)) / T(2);
            a(i, j) = average;
            a(j, i) = average;
        }
    for (std::size_t sweep = 0; sweep <= max_sweeps; ++sweep) {
        T off{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = i + 1; j < N; ++j)
                off = std::max(off, std::abs(a(i, j)));
        if (off <= tolerance) {
            result.converged = true;
            break;
        }
        if (sweep == max_sweeps)
            break;
        for (std::size_t p = 0; p < N; ++p)
            for (std::size_t q = p + 1; q < N; ++q) {
                const T apq = a(p, q);
                if (std::abs(apq) <= tolerance)
                    continue;
                // Rotation angle from the Jacobi tangent formula: one hypot and
                // one square root replace atan2 + cos + sin. hypot keeps the
                // tangent finite for very large off-diagonal ratios.
                const T theta = (a(q, q) - a(p, p)) / (T(2) * apq);
                const T t = std::copysign(T(1), theta) /
                            (std::abs(theta) + std::hypot(theta, T(1)));
                const T c = T(1) / std::sqrt(T(1) + t * t), s = t * c;
                const T app = a(p, p), aqq = a(q, q);
                for (std::size_t k = 0; k < N; ++k)
                    if (k != p && k != q) {
                        const T akp = a(k, p), akq = a(k, q);
                        a(k, p) = a(p, k) = c * akp - s * akq;
                        a(k, q) = a(q, k) = s * akp + c * akq;
                    }
                a(p, p) = c * c * app - T(2) * s * c * apq + s * s * aqq;
                a(q, q) = s * s * app + T(2) * s * c * apq + c * c * aqq;
                a(p, q) = a(q, p) = T(0);
                for (std::size_t k = 0; k < N; ++k) {
                    const T vp = result.vectors(k, p), vq = result.vectors(k, q);
                    result.vectors(k, p) = c * vp - s * vq;
                    result.vectors(k, q) = s * vp + c * vq;
                }
            }
        ++result.sweeps;
    }
    for (std::size_t i = 0; i < N; ++i)
        result.values[i] = a(i, i) * scale;
    // Ascending sort by exchanging the eigenvalue together with its eigenvector
    // column element by element; the previous form rebuilt three column
    // temporaries per exchange.
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
            if (result.values[j] < result.values[i]) {
                std::swap(result.values[i], result.values[j]);
                for (std::size_t k = 0; k < N; ++k)
                    std::swap(result.vectors(k, i), result.vectors(k, j));
            }
    if (!is_finite(result.values))
        return std::nullopt;
    return result;
}

} // namespace chm
