#pragma once
#include <chmath/chmath.hpp>
#include <limits>

namespace boundary {
inline constexpr double inf = std::numeric_limits<double>::infinity();
inline constexpr double nan = std::numeric_limits<double>::quiet_NaN();
inline constexpr double maximum = std::numeric_limits<double>::max();
inline constexpr double tiny = std::numeric_limits<double>::denorm_min();
template <class T, std::size_t N> chm::vec<T, N> sequence() {
    chm::vec<T, N> v;
    for (std::size_t i = 0; i < N; ++i)
        v[i] = static_cast<T>(i + 1);
    return v;
}
template <class T, std::size_t N> chm::mat<T, N, N> well_conditioned() {
    chm::mat<T, N, N> a;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            a(i, j) = i == j ? T(N + 2) : T(1) / T(i + j + 2);
    return a;
}
} // namespace boundary
