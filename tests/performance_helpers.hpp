#pragma once
#include "boundary_helpers.hpp"
#include "test.hpp"
#include <atomic>
#include <chrono>
#include <vector>

#if defined(_MSC_VER)
#define CH_TEST_NOINLINE __declspec(noinline)
#else
#define CH_TEST_NOINLINE __attribute__((noinline))
#endif
#ifndef CHMATH_PERFORMANCE_CHECKS
#define CHMATH_PERFORMANCE_CHECKS 0
#endif

namespace perf {
inline volatile double observed{};
template <class F> double time(const char *name, std::size_t items, F work) {
    constexpr std::size_t samples = 5, rounds = 8;
    double checksum = work();
    std::array<double, samples> elapsed{};
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t round = 0; round < rounds; ++round) {
            std::atomic_signal_fence(std::memory_order_seq_cst);
            checksum += work();
        }
        elapsed[sample] =
            std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start)
                .count() /
            double(items * rounds);
    }
    observed = checksum;
    std::sort(elapsed.begin(), elapsed.end());
    CHECK(std::isfinite(checksum));
    CHECK(elapsed[samples / 2] > 0);
    std::printf("PERF %-28s median=%10.3f ns/item min=%10.3f max=%10.3f checksum=%.9g\n", name,
                elapsed[2], elapsed[0], elapsed[4], checksum);
    return elapsed[2];
}
inline void relative(double actual, double reference, double limit) {
    const double ratio = actual / reference;
    std::printf("PERF ratio=%.3f ceiling=%.1f enforced=%d\n", ratio, limit,
                CHMATH_PERFORMANCE_CHECKS);
    CHECK(std::isfinite(ratio) && ratio > 0);
#if CHMATH_PERFORMANCE_CHECKS
    CHECK(ratio <= limit);
#endif
}
template <class T, std::size_t N>
CH_TEST_NOINLINE double matrix_kernel(const chm::mat<T, N, N> &a,
                                      std::span<const chm::mat<T, N, N>> in,
                                      std::span<chm::mat<T, N, N>> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = a * in[i];
    return double(out[in.size() / 2](0, 0));
}
template <class T, std::size_t N>
CH_TEST_NOINLINE double matrix_reference(const chm::mat<T, N, N> &a,
                                         std::span<const chm::mat<T, N, N>> in,
                                         std::span<chm::mat<T, N, N>> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        for (std::size_t row = 0; row < N; ++row)
            for (std::size_t col = 0; col < N; ++col) {
                T sum{};
                for (std::size_t k = 0; k < N; ++k)
                    sum += a(row, k) * in[i](k, col);
                out[i](row, col) = sum;
            }
    return double(out[in.size() / 2](0, 0));
}
template <class T, std::size_t N> void matrix_case() {
    using M = chm::mat<T, N, N>;
    constexpr std::size_t n = 256;
    const auto a = boundary::well_conditioned<T, N>();
    std::vector<M> in(n), out(n), reference(n);
    for (auto &m : in)
        for (auto &v : m.elements)
            v = test::sample<T>(T(-1), T(1));
    const double rt =
        time("matrix reference", n, [&] { return matrix_reference<T, N>(a, in, reference); });
    const double at = time("matrix product", n, [&] { return matrix_kernel<T, N>(a, in, out); });
    for (std::size_t i = 0; i < n; ++i)
        NEAR(out[i], reference[i], T(1e-5));
    relative(at, rt, 3.0);
}
template <class T>
CH_TEST_NOINLINE double normalize_kernel(std::span<const chm::vec<T, 3>> in,
                                         std::span<chm::vec<T, 3>> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = chm::normalize(in[i]);
    return double(out[in.size() / 2][0]);
}
template <class T>
CH_TEST_NOINLINE double normalize_reference(std::span<const chm::vec<T, 3>> in,
                                            std::span<chm::vec<T, 3>> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = in[i] / std::hypot(in[i][0], in[i][1], in[i][2]);
    return double(out[in.size() / 2][0]);
}
template <class T> void normalize_case(bool extreme = false) {
    constexpr std::size_t n = 1024;
    std::vector<chm::vec<T, 3>> in(n), out(n), ref(n);
    for (std::size_t i = 0; i < n; ++i) {
        const T scale = extreme ? (i % 2 ? std::numeric_limits<T>::max() / T(16)
                                         : std::numeric_limits<T>::min())
                                : T(1);
        in[i] = {test::sample<T>(T(1), T(2)) * scale, test::sample<T>(T(1), T(2)) * scale, scale};
    }
    const double rt =
        time("hypot normalization", n, [&] { return normalize_reference<T>(in, ref); });
    const double at = time("stable normalization", n, [&] { return normalize_kernel<T>(in, out); });
    for (std::size_t i = 0; i < n; ++i) {
        NEAR(out[i], ref[i], T(2e-6));
        NEAR(chm::length(out[i]), T(1), T(2e-6));
    }
    relative(at, rt, 4.0);
}
CH_TEST_NOINLINE inline double soa_kernel(const chm::affine3f &a, chm::const_soa3<float> in,
                                          chm::soa3<float> out) {
    if (!chm::transform_points(a, in, out))
        return std::numeric_limits<double>::quiet_NaN();
    return out.x[in.size() / 2];
}
CH_TEST_NOINLINE inline double soa_reference(const chm::affine3f &a, chm::const_soa3<float> in,
                                             chm::soa3<float> out) {
    for (std::size_t i = 0; i < in.size(); ++i) {
        const auto p = chm::transform_point(a, chm::vec3f{in.x[i], in.y[i], in.z[i]});
        out.x[i] = p[0];
        out.y[i] = p[1];
        out.z[i] = p[2];
    }
    return out.x[in.size() / 2];
}
inline void soa_case(std::size_t n) {
    using namespace chm;
    std::vector<float> x(n), y(n), z(n), ox(n), oy(n), oz(n), rx(n), ry(n), rz(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = test::sample<float>();
        y[i] = test::sample<float>();
        z[i] = test::sample<float>();
    }
    const auto a =
        compose(trs<float>{{1, 2, 3}, from_euler_xyz(vec3f{0.3f, 0.5f, 0.7f}), {2, 3, 4}});
    const const_soa3<float> in{x, y, z};
    const soa3<float> out{ox, oy, oz}, ref{rx, ry, rz};
    const double rt = time("SoA reference", n, [&] { return soa_reference(a, in, ref); });
    const double at = time("SoA transform", n, [&] { return soa_kernel(a, in, out); });
    for (std::size_t i = 0; i < n; ++i)
        NEAR((vec3f{ox[i], oy[i], oz[i]}), (vec3f{rx[i], ry[i], rz[i]}), 1e-5f);
    relative(at, rt, 5.0);
}
CH_TEST_NOINLINE inline double aos_kernel(const chm::affine3f &a, std::span<const chm::vec3f> in,
                                          std::span<chm::vec3f> out) {
    if (!chm::transform_points(a, in, out))
        return std::numeric_limits<double>::quiet_NaN();
    return out[in.size() / 2][0];
}
inline void aos_case(std::size_t n) {
    using namespace chm;
    std::vector<vec3f> in(n), out(n);
    for (auto &p : in)
        p = {test::sample<float>(), test::sample<float>(), test::sample<float>()};
    const auto a =
        compose(trs<float>{{1, 2, 3}, from_euler_xyz(vec3f{0.3f, 0.5f, 0.7f}), {2, 3, 4}});
    time("AoS transform", n, [&] { return aos_kernel(a, in, out); });
    for (std::size_t i = 0; i < n; ++i)
        NEAR(out[i], transform_point(a, in[i]), 1e-5f);
}
CH_TEST_NOINLINE inline double dot_kernel(chm::const_soa3<float> in, std::span<float> out) {
    if (!chm::dot_batch(in, in, out))
        return std::numeric_limits<double>::quiet_NaN();
    return out[in.size() / 2];
}
inline void dot_case(std::size_t n) {
    std::vector<float> x(n), y(n), z(n), out(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = test::sample<float>();
        y[i] = test::sample<float>();
        z[i] = test::sample<float>();
    }
    time("SoA dot", n, [&] { return dot_kernel({x, y, z}, out); });
    for (std::size_t i = 0; i < n; ++i)
        NEAR(out[i], x[i] * x[i] + y[i] * y[i] + z[i] * z[i], 1e-5f);
}
} // namespace perf
