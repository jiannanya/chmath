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
    if (!almost_equal(a / scale, transpose(a) / scale, tolerance, tolerance))
        return std::nullopt;
    const auto scaled = a / scale;
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
// Householder QR. Full Q is appropriate for small fixed matrices; no heap workspace.
template <floating T, std::size_t R, std::size_t C>
    requires(R >= C)
[[nodiscard]] inline std::optional<qr_factorization<T, R, C>>
factor_qr(const mat<T, R, C> &a) noexcept {
    if (!is_finite(a))
        return std::nullopt;
    qr_factorization<T, R, C> result;
    result.r = a;
    for (std::size_t k = 0; k < C; ++k) {
        vec<T, R> v;
        T scale{};
        for (std::size_t i = k; i < R; ++i)
            scale = std::max(scale, std::abs(result.r(i, k)));
        if (scale == T(0))
            continue;
        for (std::size_t i = k; i < R; ++i)
            v[i] = result.r(i, k) / scale;
        const T norm = length(v);
        v[k] += std::copysign(norm, v[k]);
        v = normalize(v);
        for (std::size_t j = k; j < C; ++j) {
            T d{};
            for (std::size_t i = k; i < R; ++i)
                d += v[i] * result.r(i, j);
            for (std::size_t i = k; i < R; ++i)
                result.r(i, j) -= T(2) * v[i] * d;
        }
        for (std::size_t i = 0; i < R; ++i) {
            T d{};
            for (std::size_t j = k; j < R; ++j)
                d += result.q(i, j) * v[j];
            for (std::size_t j = k; j < R; ++j)
                result.q(i, j) -= T(2) * d * v[j];
        }
        for (std::size_t i = k + 1; i < R; ++i)
            result.r(i, k) = T(0);
    }
    if (!is_finite(result.q) || !is_finite(result.r))
        return std::nullopt;
    return result;
}
template <floating T, std::size_t R, std::size_t C>
    requires(R >= C)
[[nodiscard]] inline std::optional<vec<T, C>>
least_squares(const mat<T, R, C> &a, const vec<T, R> &b, T tolerance = epsilon<T>) noexcept {
    if (!is_finite(b) || !is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    const auto qr = factor_qr(a);
    if (!qr)
        return std::nullopt;
    const auto rhs = transpose(qr->q) * b;
    vec<T, C> x;
    T scale{};
    for (T v : a.elements)
        scale = std::max(scale, std::abs(v));
    if (scale == T(0))
        return std::nullopt;
    for (std::size_t i = C; i-- > 0;) {
        if (std::abs(qr->r(i, i)) / scale <= tolerance)
            return std::nullopt;
        x[i] = rhs[i];
        for (std::size_t j = i + 1; j < C; ++j)
            x[i] -= qr->r(i, j) * x[j];
        x[i] /= qr->r(i, i);
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
    if (!almost_equal(a, transpose(a), tolerance, tolerance))
        return std::nullopt;
    a = (a + transpose(a)) / T(2);
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
                const T theta = T(0.5) * std::atan2(T(2) * apq, a(q, q) - a(p, p)),
                        c = std::cos(theta), s = std::sin(theta);
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
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
            if (result.values[j] < result.values[i]) {
                std::swap(result.values[i], result.values[j]);
                const auto v = result.vectors.column(i);
                result.vectors.set_column(i, result.vectors.column(j));
                result.vectors.set_column(j, v);
            }
    if (!is_finite(result.values))
        return std::nullopt;
    return result;
}

} // namespace chm
