#pragma once
#include "common.hpp"
namespace extended_vector {
using namespace chm;
template <class T, std::size_t N> void addition() {
    const auto a = extended::random_vector<T, N>(), b = extended::random_vector<T, N>();
    auto c = a;
    c += b;
    CHECK(c == a + b);
    for (std::size_t i = 0; i < N; ++i)
        NEAR(c[i], T(static_cast<long double>(a[i]) + b[i]), extended::tolerance<T>);
}
template <class T, std::size_t N> void subtraction() {
    const auto a = extended::random_vector<T, N>(), b = extended::random_vector<T, N>();
    auto c = a;
    c -= b;
    CHECK(c == a - b);
    NEAR(a - c, b, extended::tolerance<T>);
    CHECK(a - a == vec<T, N>{});
}
template <class T, std::size_t N> void scaling() {
    const auto a = extended::random_vector<T, N>();
    auto c = a;
    c *= T(-2.5);
    c /= T(-2.5);
    NEAR(c, a, extended::tolerance<T>);
    CHECK(T(2) * a == a * T(2));
    CHECK(a * T(0) == vec<T, N>{});
}
template <class T, std::size_t N> void component_product() {
    const auto a = extended::random_vector<T, N>(), b = extended::random_vector<T, N>();
    const auto c = hadamard(a, b);
    for (std::size_t i = 0; i < N; ++i)
        CHECK(c[i] == a[i] * b[i]);
    CHECK(hadamard(a, vec<T, N>{T(1)}) == a);
}
template <class T, std::size_t N> void dot_reference() {
    for (int trial = 0; trial < 32; ++trial) {
        const auto a = extended::random_vector<T, N>(), b = extended::random_vector<T, N>();
        long double oracle{};
        for (std::size_t i = 0; i < N; ++i)
            oracle += static_cast<long double>(a[i]) * b[i];
        NEAR(dot(a, b), T(oracle), extended::tolerance<T>);
        CHECK(dot(a, b) == dot(b, a));
    }
}
template <class T, std::size_t N> void normalize_reference() {
    for (int trial = 0; trial < 32; ++trial) {
        const auto a = extended::random_vector<T, N>();
        const auto n = try_normalize(a);
        CHECK(n);
        NEAR(*n, (extended::normalization_reference(a)), extended::tolerance<T>);
        NEAR(length(*n), T(1), extended::tolerance<T>);
    }
}
template <class T, std::size_t N> void projection_residual() {
    const auto a = extended::random_vector<T, N>(), b = extended::random_vector<T, N>();
    const auto p = project(a, b);
    CHECK(p);
    NEAR(dot(a - *p, b), T(0), T(4) * extended::tolerance<T>);
    const auto again = project(*p, b);
    CHECK(again);
    NEAR(*again, *p, extended::tolerance<T>);
}
template <class T, std::size_t N> void angle_invariants() {
    const auto a = extended::random_vector<T, N>(), b = extended::random_vector<T, N>();
    const auto ab = angle_between(a, b), ba = angle_between(b, a),
               scaled = angle_between(a * T(3), b * T(0.25));
    CHECK(ab && ba && scaled);
    CHECK(*ab >= T(0) && *ab <= pi<T>);
    NEAR(*ab, *ba, extended::tolerance<T>);
    NEAR(*ab, *scaled, extended::tolerance<T>);
}
template <class T, std::size_t N> void length_extremes() {
    const auto a = extended::random_vector<T, N>();
    const auto n = normalize(a);
    const int exponent = std::same_as<T, float> ? 100 : 800;
    for (int e : {-exponent, exponent}) {
        const T scale = std::ldexp(T(1), e);
        const auto v = a * scale;
        NEAR(length(v) / scale, length(a), extended::tolerance<T>);
        NEAR(normalize(v), n, extended::tolerance<T>);
    }
}
template <class T, std::size_t N> void packed_storage() {
    using V = vec<T, N>;
    static_assert(sizeof(V) == sizeof(T) * N);
    static_assert(std::is_trivially_copyable_v<V>);
    constexpr V uniform{T(3)};
    static_assert(uniform[0] == T(3));
    auto a = extended::random_vector<T, N>();
    const auto bytes = std::bit_cast<std::array<std::byte, sizeof(V)>>(a);
    CHECK(std::bit_cast<V>(bytes) == a);
    a.data()[N - 1] = T(7);
    CHECK(a[N - 1] == T(7));
    CHECK(alignof(V) == alignof(T));
}
} // namespace extended_vector
