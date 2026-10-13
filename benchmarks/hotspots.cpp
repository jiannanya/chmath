#include <atomic>
#include <charconv>
#include <chmath/chmath.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string_view>
#include <vector>

#if defined(_MSC_VER)
#define CH_HOT_NOINLINE __declspec(noinline)
#else
#define CH_HOT_NOINLINE __attribute__((noinline))
#endif

using namespace chm;
template <class F> CH_HOT_NOINLINE double run_hotspot(F &work) {
    return work();
}
template <class F>
void measure_hotspot(const char *name, std::size_t count, F work, double &checksum) {
    for (int i = 0; i < 3; ++i)
        checksum += run_hotspot(work);
    std::array<double, 9> times{};
    for (double &time : times) {
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < 12; ++i) {
            std::atomic_signal_fence(std::memory_order_seq_cst);
            checksum += run_hotspot(work);
        }
        time = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start)
                   .count() /
               double(count * 12);
    }
    std::sort(times.begin(), times.end());
    std::printf("%-32s %10.3f ns/item  min=%10.3f max=%10.3f\n", name, times[4], times[0],
                times[8]);
}

template <floating T>
void batch_hotspots(std::size_t n, const char *transform_name, const char *normalize_name,
                    const char *mixed_name, double &checksum) {
    std::vector<vec<T, 3>> input(n), output(n);
    std::mt19937 random(20261013);
    std::uniform_real_distribution<T> sample(T(-10), T(10));
    for (auto &v : input)
        v = {sample(random), sample(random), sample(random)};
    const auto a = compose(trs<T>{{T(1), T(2), T(3)},
                                  from_euler_xyz(vec<T, 3>{T(0.3), T(0.7), T(-0.9)}),
                                  {T(2), T(3), T(4)}});
    const std::span<const vec<T, 3>> in{input};
    const std::span<vec<T, 3>> out{output};
    measure_hotspot(
        transform_name, n,
        [&] {
            if (!transform_points(a, in, out))
                std::abort();
            return double(output[n / 2][0]);
        },
        checksum);
    measure_hotspot(
        normalize_name, n,
        [&] {
            if (!normalize_vectors(in, out))
                std::abort();
            return double(output[n / 2][0]);
        },
        checksum);
    // Occasional zero/extreme vectors exercise the stable fallback in a packet.
    for (std::size_t i = 0; i < n; i += 17)
        input[i] = i % 2 == 0 ? vec<T, 3>{} : vec<T, 3>{std::numeric_limits<T>::max()};
    measure_hotspot(
        mixed_name, n,
        [&] {
            if (!normalize_vectors(in, out))
                std::abort();
            return double(output[n / 2][0]);
        },
        checksum);
    for (std::size_t i = 0; i < n; ++i)
        if (!almost_equal(output[i], normalize(input[i]), T(1e-5), T(1e-5)))
            std::abort();
}

template <floating T, std::size_t R, std::size_t C>
void least_squares_hotspot(const char *name, std::size_t count, double &checksum) {
    std::vector<mat<T, R, C>> matrices(count);
    std::vector<vec<T, R>> observations(count);
    std::vector<vec<T, C>> output(count);
    std::mt19937 random(20261013);
    std::uniform_real_distribution<T> sample(T(-1), T(1));
    for (std::size_t i = 0; i < count; ++i) {
        for (T &v : matrices[i].elements)
            v = sample(random);
        for (T &v : observations[i].elements)
            v = sample(random);
    }
    measure_hotspot(
        name, count,
        [&] {
            for (std::size_t i = 0; i < count; ++i) {
                const auto x = least_squares(matrices[i], observations[i]);
                if (!x)
                    std::abort();
                output[i] = *x;
            }
            return double(output[count / 2][0]);
        },
        checksum);
    for (std::size_t i = 0; i < count; ++i)
        if (!almost_equal(transpose(matrices[i]) * (matrices[i] * output[i] - observations[i]),
                          vec<T, C>{}, T(1e-3), T(1e-3)))
            std::abort();
}

template <std::size_t R, std::size_t C>
void qr_hotspot(const char *name, std::size_t count, double &checksum) {
    std::vector<mat<double, R, C>> input(count);
    std::vector<qr_factorization<double, R, C>> output(count);
    std::mt19937 random(20261013);
    std::uniform_real_distribution<double> sample(-1, 1);
    for (auto &a : input)
        for (double &v : a.elements)
            v = sample(random);
    measure_hotspot(
        name, count,
        [&] {
            for (std::size_t i = 0; i < count; ++i) {
                const auto qr = factor_qr(input[i]);
                if (!qr)
                    std::abort();
                output[i] = *qr;
            }
            return output[count / 2].q(0, 0) + output[count / 2].r(0, 0);
        },
        checksum);
    if (!almost_equal(output[0].q * output[0].r, input[0], 1e-10, 1e-10))
        std::abort();
}

template <std::size_t N> void lu_hotspot(const char *name, std::size_t count, double &checksum) {
    std::vector<mat<double, N, N>> input(count);
    std::vector<lu_factorization<double, N>> output(count);
    std::mt19937 random(20261013);
    std::uniform_real_distribution<double> sample(-1, 1);
    for (auto &a : input) {
        for (double &v : a.elements)
            v = sample(random);
        for (std::size_t i = 0; i < N; ++i)
            a(i, i) += double(N);
    }
    measure_hotspot(
        name, count,
        [&] {
            for (std::size_t i = 0; i < count; ++i) {
                const auto lu = factor_lu(input[i]);
                if (!lu)
                    std::abort();
                output[i] = *lu;
            }
            return output[count / 2].factors(0, 0);
        },
        checksum);
    const vec<double, N> rhs{1.0};
    if (!almost_equal(input[0] * output[0].solve(rhs), rhs, 1e-10, 1e-10))
        std::abort();
}

int main(int argc, char **argv) {
    std::size_t n = 65536;
    if (argc > 1) {
        const std::string_view arg{argv[1]};
        const auto parsed = std::from_chars(arg.data(), arg.data() + arg.size(), n);
        if (parsed.ec != std::errc{} || parsed.ptr != arg.data() + arg.size() || n < 16 ||
            n > (1U << 24))
            return 1;
    }
    std::printf("hotspots backend=%s points=%zu samples=9 rounds=12 warmup=3\n", batch_backend(),
                n);
    double checksum{};
    batch_hotspots<float>(n, "AoS float transform", "AoS float normalize",
                          "AoS float normalize mixed", checksum);
    batch_hotspots<double>(n, "AoS double transform", "AoS double normalize",
                           "AoS double normalize mixed", checksum);
    const auto small = std::min(n, std::size_t(256));
    least_squares_hotspot<double, 4, 2>("least squares double 4x2", small, checksum);
    least_squares_hotspot<float, 8, 3>("least squares float 8x3", small, checksum);
    least_squares_hotspot<double, 32, 4>("least squares double 32x4", small, checksum);
    least_squares_hotspot<double, 64, 8>("least squares double 64x8", small, checksum);
    qr_hotspot<8, 3>("QR double 8x3", small, checksum);
    qr_hotspot<32, 4>("QR double 32x4", small, checksum);
    lu_hotspot<8>("LU double 8x8", small, checksum);
    lu_hotspot<16>("LU double 16x16", small, checksum);
    lu_hotspot<32>("LU double 32x32", small, checksum);
    std::printf("checksum=%.12g\n", checksum);
    return is_finite(checksum) ? 0 : 2;
}
