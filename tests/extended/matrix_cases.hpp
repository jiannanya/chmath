#pragma once
#include "common.hpp"
namespace extended_matrix {
using namespace chm;
template <class T, std::size_t N> void product_oracle() {
    auto a = extended::positive_matrix<T, N>(), b = a;
    for (T &v : b.elements)
        v = test::sample<T>(T(-2), T(2));
    const auto c = a * b;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) {
            long double sum{};
            for (std::size_t k = 0; k < N; ++k)
                sum += static_cast<long double>(a(i, k)) * b(k, j);
            NEAR(c(i, j), T(sum), extended::tolerance<T>);
        }
}
template <class T, std::size_t N> void transpose_contract() {
    auto a = extended::positive_matrix<T, N>();
    a.set_column(0, extended::random_vector<T, N>());
    const auto t = transpose(a);
    CHECK(transpose(t) == a);
    for (std::size_t i = 0; i < N; ++i)
        CHECK(t.row(i) == a.column(i));
}
template <class T, std::size_t N> void arithmetic() {
    const auto a = extended::positive_matrix<T, N>();
    auto b = a;
    b += a;
    b -= a;
    CHECK(a == b);
    NEAR((a * T(3) - a) / T(2), a, extended::tolerance<T>);
    CHECK(a - a == mat<T, N, N>{});
}
template <class T, std::size_t N> void triangular_determinant() {
    mat<T, N, N> a;
    T expected = T(1);
    for (std::size_t i = 0; i < N; ++i) {
        a(i, i) = T(i + 1);
        expected *= T(i + 1);
        for (std::size_t j = i + 1; j < N; ++j)
            a(i, j) = T(i + j + 2);
    }
    NEAR(determinant(a), expected, extended::tolerance<T>);
    if constexpr (N > 1) {
        for (std::size_t j = 0; j < N; ++j)
            std::swap(a(0, j), a(N - 1, j));
        NEAR(determinant(a), -expected, extended::tolerance<T>);
    }
}
template <class T, std::size_t N> void solve_known_solution() {
    const auto a = extended::positive_matrix<T, N>();
    const auto x = extended::random_vector<T, N>();
    const auto b = a * x;
    const auto actual = solve(a, b);
    CHECK(actual);
    NEAR(*actual, x, extended::tolerance<T>);
    NEAR(a * *actual, b, extended::tolerance<T>);
}
template <class T, std::size_t N> void multiple_rhs() {
    const auto a = extended::positive_matrix<T, N>();
    mat<T, N, 7> x;
    for (T &v : x.elements)
        v = test::sample<T>(T(-2), T(2));
    const auto rhs = a * x;
    const auto lu = factor_lu(a);
    CHECK(lu);
    const auto actual = solve(a, rhs);
    CHECK(actual);
    NEAR(*actual, x, extended::tolerance<T>);
    NEAR(lu->solve(rhs), x, extended::tolerance<T>);
    for (std::size_t j = 0; j < 7; ++j)
        NEAR(actual->column(j), lu->solve(rhs.column(j)), extended::tolerance<T>);
}
template <class T, std::size_t N> void two_sided_inverse() {
    const auto a = extended::positive_matrix<T, N>();
    const auto inv = inverse(a);
    CHECK(inv);
    NEAR(a * *inv, (mat<T, N, N>::identity()), extended::tolerance<T>);
    NEAR(*inv * a, (mat<T, N, N>::identity()), extended::tolerance<T>);
}
template <class T, std::size_t N> void small_scale_lu() {
    const T scale = std::same_as<T, float> ? T(1e-30L) : T(1e-250L);
    const auto a = extended::positive_matrix<T, N>() * scale;
    const vec<T, N> x{T(0.5)};
    const auto got = solve(a, a * x);
    CHECK(got);
    NEAR(*got, x, extended::tolerance<T>);
    const auto inv = inverse(a);
    CHECK(inv);
    NEAR(a * *inv, (mat<T, N, N>::identity()), extended::tolerance<T>);
}
template <class T, std::size_t N> void qr_reconstruction() {
    auto a = extended::positive_matrix<T, N>();
    a.set_column(0, -a.column(0));
    const auto qr = factor_qr(a);
    CHECK(qr);
    NEAR(qr->q * qr->r, a, T(4) * extended::tolerance<T>);
    NEAR(transpose(qr->q) * qr->q, (mat<T, N, N>::identity()), T(4) * extended::tolerance<T>);
}
template <class T, std::size_t N> void cholesky_reconstruction() {
    const auto p = extended::positive_matrix<T, N>();
    const auto a = transpose(p) * p + mat<T, N, N>::identity();
    const auto l = cholesky(a);
    CHECK(l);
    NEAR(*l * transpose(*l), a, extended::tolerance<T>);
    for (std::size_t i = 0; i < N; ++i) {
        CHECK((*l)(i, i) > T(0));
        for (std::size_t j = i + 1; j < N; ++j)
            CHECK((*l)(i, j) == T(0));
    }
}
template <class T, std::size_t N> void eigen_residuals() {
    const auto a = extended::positive_matrix<T, N>();
    const auto e = symmetric_eigen(a, epsilon<T>, 64);
    CHECK(e && e->converged);
    for (std::size_t i = 0; i < N; ++i) {
        NEAR(a * e->vectors.column(i), e->vectors.column(i) * e->values[i],
             T(8) * extended::tolerance<T>);
        if (i)
            CHECK(e->values[i] >= e->values[i - 1]);
    }
    NEAR(transpose(e->vectors) * e->vectors, (mat<T, N, N>::identity()),
         T(8) * extended::tolerance<T>);
}
template <class T, std::size_t N> void self_assignment() {
    auto a = extended::positive_matrix<T, N>();
    const auto expected = a * a;
    a = a * a;
    CHECK(a == expected);
    const auto before = a;
    a = a * mat<T, N, N>::identity();
    CHECK(a == before);
}
} // namespace extended_matrix
