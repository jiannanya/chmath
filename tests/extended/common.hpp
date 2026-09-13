#pragma once
#include "../test.hpp"
#include <bit>
#include <chmath/chmath.hpp>
#include <thread>
#include <vector>

namespace extended {
using namespace chm;
template <class T> constexpr T tolerance = std::same_as<T, float> ? T(2e-4) : T(2e-11);
template <class T, std::size_t N> vec<T, N> random_vector() {
    vec<T, N> v;
    for (T &x : v.elements)
        x = test::sample<T>(T(-2), T(2));
    return v;
}
template <class T, std::size_t N> mat<T, N, N> positive_matrix() {
    mat<T, N, N> a;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            a(i, j) = i == j ? T(N + 1) : T(1) / T(i + j + 2);
    return a;
}
template <class T, std::size_t N> vec<T, N> normalization_reference(const vec<T, N> &v) {
    long double norm{};
    for (T x : v.elements)
        norm = std::hypot(norm, static_cast<long double>(x));
    vec<T, N> out;
    if (!std::isfinite(norm) || norm == 0)
        return out;
    for (std::size_t i = 0; i < N; ++i)
        out[i] = T(static_cast<long double>(v[i]) / norm);
    return out;
}
} // namespace extended
