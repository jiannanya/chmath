#pragma once

#include <chmath/simd/detail/matrix_kernel.hpp>
#include <chmath/vector/vector.hpp>

namespace chm {

// Column-major memory, column vectors: (A*B)*v applies B first.
template <scalar T, std::size_t R, std::size_t C>
    requires(R > 0 && C > 0)
struct mat {
    std::array<T, R * C> elements{};
    using value_type = T;
    static constexpr std::size_t rows = R, columns = C;
    constexpr mat() noexcept = default;
    explicit constexpr mat(T diagonal) noexcept {
        for (std::size_t i = 0; i < std::min(R, C); ++i)
            (*this)(i, i) = diagonal;
    }
    explicit constexpr mat(const std::array<T, R * C> &column_major) noexcept
        : elements(column_major) {}
    template <scalar U> explicit constexpr mat(const mat<U, R, C> &other) noexcept {
        for (std::size_t i = 0; i < R * C; ++i)
            elements[i] = static_cast<T>(other.elements[i]);
    }
    [[nodiscard]] static constexpr mat identity() noexcept
        requires(R == C)
    {
        return mat(T(1));
    }
    [[nodiscard]] static constexpr mat from_rows(const std::array<vec<T, C>, R> &values) noexcept {
        mat r;
        for (std::size_t i = 0; i < R; ++i)
            for (std::size_t j = 0; j < C; ++j)
                r(i, j) = values[i][j];
        return r;
    }
    [[nodiscard]] constexpr T &operator()(std::size_t row, std::size_t col) noexcept {
        assert(row < R && col < C);
        return elements[col * R + row];
    }
    [[nodiscard]] constexpr const T &operator()(std::size_t row, std::size_t col) const noexcept {
        assert(row < R && col < C);
        return elements[col * R + row];
    }
    [[nodiscard]] constexpr T *data() noexcept { return elements.data(); }
    [[nodiscard]] constexpr const T *data() const noexcept { return elements.data(); }
    [[nodiscard]] constexpr vec<T, C> row(std::size_t i) const noexcept {
        vec<T, C> r;
        for (std::size_t j = 0; j < C; ++j)
            r[j] = (*this)(i, j);
        return r;
    }
    [[nodiscard]] constexpr vec<T, R> column(std::size_t j) const noexcept {
        vec<T, R> r;
        for (std::size_t i = 0; i < R; ++i)
            r[i] = (*this)(i, j);
        return r;
    }
    constexpr void set_column(std::size_t j, const vec<T, R> &v) noexcept {
        for (std::size_t i = 0; i < R; ++i)
            (*this)(i, j) = v[i];
    }
    constexpr mat &operator+=(const mat &b) noexcept {
        for (std::size_t i = 0; i < R * C; ++i)
            elements[i] += b.elements[i];
        return *this;
    }
    constexpr mat &operator-=(const mat &b) noexcept {
        for (std::size_t i = 0; i < R * C; ++i)
            elements[i] -= b.elements[i];
        return *this;
    }
    constexpr mat &operator*=(T s) noexcept {
        for (auto &x : elements)
            x *= s;
        return *this;
    }
    constexpr mat &operator/=(T s) noexcept {
        for (auto &x : elements)
            x /= s;
        return *this;
    }
    [[nodiscard]] constexpr bool operator==(const mat &) const noexcept = default;
};

template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator+(mat<T, R, C> a, const mat<T, R, C> &b) noexcept {
    return a += b;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator-(mat<T, R, C> a, const mat<T, R, C> &b) noexcept {
    return a -= b;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator-(mat<T, R, C> a) noexcept {
    for (auto &x : a.elements)
        x = -x;
    return a;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator*(mat<T, R, C> a, T s) noexcept {
    return a *= s;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator*(T s, mat<T, R, C> a) noexcept {
    return a *= s;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator/(mat<T, R, C> a, T s) noexcept {
    return a /= s;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr vec<T, R> operator*(const mat<T, R, C> &a, const vec<T, C> &b) noexcept {
    if constexpr (std::same_as<T, float> && R == 4 && C == 4 && detail::native_mat4_kernel) {
        if (!std::is_constant_evaluated()) {
            vec<T, R> r;
            detail::multiply_mat4_vec4_float(a.data(), b.data(), r.data());
            return r;
        }
    }
    // One local accumulator per row, seeded with the first product so the whole
    // dot product is a single dependency chain assigned once. This avoids the
    // read-modify-write of a zero-initialized result element for every column.
    vec<T, R> r;
    for (std::size_t i = 0; i < R; ++i) {
        T sum = a(i, 0) * b[0];
        for (std::size_t j = 1; j < C; ++j)
            sum += a(i, j) * b[j];
        r[i] = sum;
    }
    return r;
}
template <scalar T, std::size_t R, std::size_t K, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> operator*(const mat<T, R, K> &a,
                                               const mat<T, K, C> &b) noexcept {
    mat<T, R, C> r;
    if constexpr (std::same_as<T, float> && R == 4 && K == 4 && C == 4 &&
                  detail::native_mat4_kernel) {
        if (!std::is_constant_evaluated()) {
            detail::multiply_mat4_float(a.data(), b.data(), r.data());
            return r;
        }
    }
    // Small products keep the whole dot product in one local accumulator and
    // assign once. Seeding it with the first product (instead of zero) lets the
    // remaining terms contract into a single fused chain, which is both faster
    // and gives the product a fixed evaluation shape independent of how the
    // caller consumes it.
    if constexpr (R <= 4 && K <= 4 && C <= 4) {
        for (std::size_t j = 0; j < C; ++j)
            for (std::size_t i = 0; i < R; ++i) {
                T sum = a(i, 0) * b(0, j);
                for (std::size_t k = 1; k < K; ++k)
                    sum += a(i, k) * b(k, j);
                r(i, j) = sum;
            }
    } else {
        for (std::size_t j = 0; j < C; ++j)
            for (std::size_t k = 0; k < K; ++k)
                for (std::size_t i = 0; i < R; ++i)
                    r(i, j) += a(i, k) * b(k, j);
    }
    return r;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, C, R> transpose(const mat<T, R, C> &a) noexcept {
    mat<T, C, R> r;
    for (std::size_t j = 0; j < C; ++j)
        for (std::size_t i = 0; i < R; ++i)
            r(j, i) = a(i, j);
    return r;
}
template <scalar T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr mat<T, R, C> outer_product(const vec<T, R> &a,
                                                   const vec<T, C> &b) noexcept {
    mat<T, R, C> r;
    for (std::size_t j = 0; j < C; ++j)
        for (std::size_t i = 0; i < R; ++i)
            r(i, j) = a[i] * b[j];
    return r;
}
template <scalar T, std::size_t N> [[nodiscard]] constexpr T trace(const mat<T, N, N> &a) noexcept {
    T r{};
    for (std::size_t i = 0; i < N; ++i)
        r += a(i, i);
    return r;
}
template <floating T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr bool is_finite(const mat<T, R, C> &a) noexcept {
    for (T x : a.elements)
        if (!is_finite(x))
            return false;
    return true;
}
template <floating T, std::size_t R, std::size_t C>
[[nodiscard]] constexpr bool almost_equal(const mat<T, R, C> &a, const mat<T, R, C> &b,
                                          T rel = epsilon<T>, T abs = epsilon<T>) noexcept {
    if (!tolerance_valid(rel, abs))
        return false;
    for (std::size_t i = 0; i < R * C; ++i)
        if (!detail::almost_equal_value(a.elements[i], b.elements[i], rel, abs))
            return false;
    return true;
}

// Scaled partial pivoting; scale-independent singularity policy, fixed stack storage.
namespace detail {
template <floating T> [[nodiscard]] constexpr T matrix_magnitude(T value) noexcept {
    if (std::is_constant_evaluated())
        return value < T(0) ? -value : value;
    return std::abs(value);
}
template <floating T, std::size_t N>
[[nodiscard]] inline T scaled_diagonal_product(const mat<T, N, N> &factors, int parity,
                                               T scale = T(1)) noexcept {
    T mantissa = static_cast<T>(parity);
    long long exponent = 0;
    int scale_exponent{};
    const T scale_fraction = std::frexp(scale, &scale_exponent);
    for (std::size_t i = 0; i < N; ++i) {
        int e{}, adjustment{};
        const T fraction = std::frexp(factors(i, i), &e);
        mantissa = std::frexp(mantissa * fraction * scale_fraction, &adjustment);
        exponent += static_cast<long long>(e) + scale_exponent + adjustment;
    }
    if (exponent > std::numeric_limits<T>::max_exponent + 1LL)
        return std::copysign(std::numeric_limits<T>::infinity(), mantissa);
    if (exponent < std::numeric_limits<T>::min_exponent - std::numeric_limits<T>::digits - 1LL)
        return std::copysign(T(0), mantissa);
    return std::ldexp(mantissa, static_cast<int>(exponent));
}
} // namespace detail
template <floating T, std::size_t N> struct lu_factorization {
    mat<T, N, N> factors{};
    std::array<std::size_t, N> permutation{};
    int parity = 1;

    [[nodiscard]] constexpr vec<T, N> solve(const vec<T, N> &b) const noexcept {
        vec<T, N> x;
        for (std::size_t i = 0; i < N; ++i) {
            x[i] = b[permutation[i]];
            for (std::size_t j = 0; j < i; ++j)
                x[i] -= factors(i, j) * x[j];
        }
        for (std::size_t i = N; i-- > 0;) {
            for (std::size_t j = i + 1; j < N; ++j)
                x[i] -= factors(i, j) * x[j];
            x[i] /= factors(i, i);
        }
        return x;
    }
    // Multiple right-hand sides reuse the diagonal reciprocals. A reciprocal
    // that cannot be represented falls back to division (including 0/tiny).
    template <std::size_t M>
    [[nodiscard]] constexpr mat<T, N, M> solve(const mat<T, N, M> &b) const noexcept {
        mat<T, N, M> x;
        std::array<T, N> reciprocal{};
        for (std::size_t i = 0; i < N; ++i)
            reciprocal[i] = std::is_constant_evaluated() ? T(0) : T(1) / factors(i, i);
        for (std::size_t col = 0; col < M; ++col) {
            for (std::size_t i = 0; i < N; ++i) {
                T value = b(permutation[i], col);
                for (std::size_t j = 0; j < i; ++j)
                    value -= factors(i, j) * x(j, col);
                x(i, col) = value;
            }
            for (std::size_t i = N; i-- > 0;) {
                T value = x(i, col);
                for (std::size_t j = i + 1; j < N; ++j)
                    value -= factors(i, j) * x(j, col);
                x(i, col) = std::is_constant_evaluated() || !is_finite(reciprocal[i])
                                ? value / factors(i, i)
                                : value * reciprocal[i];
            }
        }
        return x;
    }
    [[nodiscard]] constexpr T determinant() const noexcept {
        T r = static_cast<T>(parity);
        for (std::size_t i = 0; i < N; ++i) {
            r *= factors(i, i);
            if (!is_finite(r) || std::abs(r) < std::numeric_limits<T>::min())
                return detail::scaled_diagonal_product(factors, parity);
        }
        return r;
    }
};

template <floating T, std::size_t N>
[[nodiscard]] constexpr std::optional<lu_factorization<T, N>>
factor_lu(const mat<T, N, N> &a, T tolerance = epsilon<T>) noexcept {
    if (!is_finite(tolerance) || tolerance < T(0))
        return std::nullopt;
    lu_factorization<T, N> lu;
    lu.factors = a;
    std::array<T, N> scales{}, reciprocal_scales{};
    for (std::size_t i = 0; i < N; ++i) {
        lu.permutation[i] = i;
        for (std::size_t j = 0; j < N; ++j) {
            if (!is_finite(a(i, j)))
                return std::nullopt;
            scales[i] = std::max(scales[i], detail::matrix_magnitude(a(i, j)));
        }
        if (scales[i] == T(0))
            return std::nullopt;
        reciprocal_scales[i] = std::is_constant_evaluated() ? T(0) : T(1) / scales[i];
    }
    for (std::size_t k = 0; k < N; ++k) {
        std::size_t pivot = k;
        T best = T(-1);
        for (std::size_t i = k; i < N; ++i) {
            T ratio{};
            if (std::is_constant_evaluated() || !is_finite(reciprocal_scales[i]))
                ratio = detail::matrix_magnitude(lu.factors(i, k)) / scales[i];
            else
                ratio = detail::matrix_magnitude(lu.factors(i, k)) * reciprocal_scales[i];
            if (ratio > best) {
                best = ratio;
                pivot = i;
            }
        }
        if (best <= tolerance)
            return std::nullopt;
        if (pivot != k) {
            for (std::size_t j = 0; j < N; ++j)
                std::swap(lu.factors(k, j), lu.factors(pivot, j));
            std::swap(scales[k], scales[pivot]);
            std::swap(reciprocal_scales[k], reciprocal_scales[pivot]);
            std::swap(lu.permutation[k], lu.permutation[pivot]);
            lu.parity = -lu.parity;
        }
        T reciprocal_pivot{};
        if (!std::is_constant_evaluated())
            reciprocal_pivot = T(1) / lu.factors(k, k);
        for (std::size_t i = k + 1; i < N; ++i) {
            lu.factors(i, k) = std::is_constant_evaluated() || !is_finite(reciprocal_pivot)
                                   ? lu.factors(i, k) / lu.factors(k, k)
                                   : lu.factors(i, k) * reciprocal_pivot;
            for (std::size_t j = k + 1; j < N; ++j)
                lu.factors(i, j) -= lu.factors(i, k) * lu.factors(k, j);
        }
    }
    if (!is_finite(lu.factors))
        return std::nullopt;
    return lu;
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr std::optional<vec<T, N>> solve(const mat<T, N, N> &a, const vec<T, N> &b,
                                                       T tolerance = epsilon<T>) noexcept {
    if (!is_finite(b))
        return std::nullopt;
    const auto lu = factor_lu(a, tolerance);
    if (!lu)
        return std::nullopt;
    const auto x = lu->solve(b);
    if (!is_finite(x))
        return std::nullopt;
    return x;
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr std::optional<mat<T, N, N>> inverse(const mat<T, N, N> &a,
                                                            T tolerance = epsilon<T>) noexcept {
    const auto lu = factor_lu(a, tolerance);
    if (!lu)
        return std::nullopt;
    const auto r = lu->solve(mat<T, N, N>::identity());
    if (!is_finite(r))
        return std::nullopt;
    return r;
}

template <floating T, std::size_t N, std::size_t M>
[[nodiscard]] constexpr std::optional<mat<T, N, M>>
solve(const mat<T, N, N> &a, const mat<T, N, M> &b, T tolerance = epsilon<T>) noexcept {
    if (!is_finite(b))
        return std::nullopt;
    const auto lu = factor_lu(a, tolerance);
    if (!lu)
        return std::nullopt;
    const auto x = lu->solve(b);
    return is_finite(x) ? std::optional{x} : std::nullopt;
}
template <floating T, std::size_t N>
[[nodiscard]] constexpr T determinant(const mat<T, N, N> &a) noexcept {
    // A zero pivot means singular; no tolerance truncation of small determinants.
    if (!is_finite(a))
        return std::numeric_limits<T>::quiet_NaN();
    const auto lu = factor_lu(a, T(0));
    if (lu)
        return lu->determinant();
    // Retry in a bounded range when elimination itself overflowed. This also
    // preserves true singular matrices as zero instead of fabricating an inverse.
    T scale{};
    for (T value : a.elements)
        scale = std::max(scale, std::abs(value));
    if (scale == T(0))
        return T(0);
    const auto scaled = factor_lu(a / scale, T(0));
    return scaled ? detail::scaled_diagonal_product(scaled->factors, scaled->parity, scale) : T(0);
}

template <class T> using mat2 = mat<T, 2, 2>;
template <class T> using mat3 = mat<T, 3, 3>;
template <class T> using mat4 = mat<T, 4, 4>;
using mat2f = mat2<float>;
using mat3f = mat3<float>;
using mat4f = mat4<float>;
using mat2d = mat2<double>;
using mat3d = mat3<double>;
using mat4d = mat4<double>;

} // namespace chm
