// Parallel batch benchmark. This executable is compiled with
// CHMATH_ENABLE_PARALLEL=1 so the worker-count overloads actually spawn
// threads; every row is measured with workers=1 and with the hardware worker
// count on the same data, and the parallel results are verified against the
// serial ones.
#include <algorithm>
#include <array>
#include <charconv>
#include <chmath/chmath.hpp>
#include <chrono>
#include <cstdio>
#include <limits>
#include <random>
#include <span>
#include <string_view>
#include <vector>

#if defined(_MSC_VER)
#define CH_PARALLEL_BENCH_NOINLINE __declspec(noinline)
#else
#define CH_PARALLEL_BENCH_NOINLINE __attribute__((noinline))
#endif
using namespace chm;

#if !CHMATH_ENABLE_PARALLEL
#error "chmath_bench_parallel must be compiled with CHMATH_ENABLE_PARALLEL=1"
#endif

namespace {
template <class F> CH_PARALLEL_BENCH_NOINLINE double invoke_work(F &f, unsigned workers) {
    return f(workers);
}

template <class F>
double measure(const char *name, std::size_t items, unsigned workers, F &work, double &checksum) {
    for (int i = 0; i < 2; ++i)
        checksum += invoke_work(work, workers);
    std::array<double, 9> times{};
    for (auto &sample : times) {
        const auto start = std::chrono::steady_clock::now();
        constexpr int rounds = 6;
        for (int j = 0; j < rounds; ++j)
            checksum += invoke_work(work, workers);
        sample = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start)
                     .count() /
                 double(items * rounds);
    }
    std::sort(times.begin(), times.end());
    std::printf("%-30s workers=%-3u %9.3f ns/item  %8.3f Mitems/s  min=%9.3f max=%9.3f\n", name,
                workers, times[4], 1000.0 / times[4], times[0], times[8]);
    return times[4];
}

template <class F>
void compare(const char *name, std::size_t items, unsigned workers, F work, double &checksum) {
    const double serial = measure(name, items, 1U, work, checksum);
    const double parallel = measure(name, items, workers, work, checksum);
    std::printf("%-30s speedup %.2fx over workers=1\n", name, serial / parallel);
}
} // namespace

int main(int argc, char **argv) {
    std::size_t n = std::size_t(1) << 20;
    if (argc > 1) {
        const std::string_view text = argv[1];
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), n);
        if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
            n < (std::size_t(1) << 16) || n > (std::size_t(1) << 26)) {
            std::fprintf(stderr, "usage: chmath_bench_parallel [item_count: 65536..67108864]\n");
            return 1;
        }
    }
    const unsigned workers = parallel_default_workers();
    std::mt19937 random(20260921);
    std::uniform_real_distribution<double> dist(-4, 4);
    std::vector<double> x(n), y(n), z(n), ox(n), oy(n), oz(n), rx(n), ry(n), rz(n), dots(n);
    std::vector<vec3d> packed(n), packed_ref(n), packed_out(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = dist(random);
        y[i] = dist(random);
        z[i] = dist(random);
        packed[i] = {x[i], y[i], z[i]};
    }
    std::printf("parallel benchmark backend=%s points=%zu hardware_threads=%u min_chunk=%zu\n",
                batch_backend(), n, workers, parallel_min_chunk);
    if (workers < 2U)
        std::printf("NOTE: only one hardware thread reported; parallel rows stay serial.\n");
    const auto model =
        compose(trs<double>{{1, 2, 3}, from_euler_xyz(vec3d{0.3, 0.7, -0.9}), {2, 3, 4}});
    const const_soa3<double> input{x, y, z};
    const soa3<double> output{ox, oy, oz}, reference{rx, ry, rz};
    double checksum{};

    if (!transform_points(model, input, reference, 1U))
        return 2;
    compare(
        "SoA double transform", n, workers,
        [&](unsigned w) {
            if (!transform_points(model, input, output, w))
                std::abort();
            return ox[n / 2];
        },
        checksum);
    for (std::size_t i = 0; i < n; ++i)
        if (!almost_equal(vec3d{ox[i], oy[i], oz[i]}, vec3d{rx[i], ry[i], rz[i]}, 1e-12, 1e-12)) {
            std::fprintf(stderr, "parallel SoA transform mismatch at %zu\n", i);
            return 3;
        }
    compare(
        "SoA double dot", n, workers,
        [&](unsigned w) {
            if (!dot_batch(input, input, std::span<double>{dots}, w))
                std::abort();
            return dots[n / 2];
        },
        checksum);
    if (!normalize_vectors(std::span<const vec3d>{packed}, std::span<vec3d>{packed_ref}, 1U))
        return 4;
    compare(
        "AoS double normalize", n, workers,
        [&](unsigned w) {
            if (!normalize_vectors(std::span<const vec3d>{packed},
                                   std::span<vec3d>{packed_out}, w))
                std::abort();
            return packed_out[n / 2][0];
        },
        checksum);
    for (std::size_t i = 0; i < n; ++i)
        if (!almost_equal(packed_out[i], packed_ref[i], 1e-12, 1e-12)) {
            std::fprintf(stderr, "parallel AoS normalize mismatch at %zu\n", i);
            return 5;
        }
    const ray3d ray{{0, 0, 10}, {0.2, 0.3, -1}};
    const auto prepared = *prepared_ray<double>::create(ray);
    std::vector<aabb3d> boxes(n);
    std::vector<std::optional<ray_interval<double>>> hits(n);
    for (std::size_t i = 0; i < n; ++i)
        boxes[i] = {packed[i] - vec3d{1}, packed[i] + vec3d{1}};
    compare(
        "prepared ray AABB queries", n, workers,
        [&](unsigned w) {
            if (!intersect_many(prepared, std::span<const aabb3d>{boxes},
                                std::span<std::optional<ray_interval<double>>>{hits}, 0.0,
                                std::numeric_limits<double>::infinity(), w))
                std::abort();
            return hits[n / 2] ? hits[n / 2]->enter : -1.0;
        },
        checksum);
    const auto frustum_view = *extract_frustum(*perspective(1.0, 1.5, 0.1, 100.0));
    std::vector<aabb3d> visible_boxes(n);
    std::vector<containment> visible(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = double(i % 17) - 8;
        const double depth = double(i % 103);
        visible_boxes[i] = {vec3d{t, t, -depth}, vec3d{t + 1, t + 1, 1 - depth}};
    }
    compare(
        "frustum classification", n, workers,
        [&](unsigned w) {
            if (!classify_batch(frustum_view, std::span<const aabb3d>{visible_boxes},
                                std::span<containment>{visible}, w))
                std::abort();
            return double(int(visible[n / 2]));
        },
        checksum);
    std::printf("checksum=%.12g\n", checksum);
    return is_finite(checksum) ? 0 : 6;
}
