#include <algorithm>
#include <array>
#include <charconv>
#include <chmath/chmath.hpp>
#include <chrono>
#include <cstdio>
#include <random>
#include <string_view>
#include <vector>

#if defined(_MSC_VER)
#define BENCH_NOINLINE __declspec(noinline)
#else
#define BENCH_NOINLINE __attribute__((noinline))
#endif
using namespace chm;

// No-inline call boundaries and consumed outputs prevent loop deletion/hoisting.
BENCH_NOINLINE void baseline_soa(const affine3f &a, const_soa3<float> input, soa3<float> output) {
    for (std::size_t i = 0; i < input.size(); ++i) {
        const float x = input.x[i], y = input.y[i], z = input.z[i];
        output.x[i] =
            ((a.matrix(0, 0) * x + a.matrix(0, 1) * y) + a.matrix(0, 2) * z) + a.matrix(0, 3);
        output.y[i] =
            ((a.matrix(1, 0) * x + a.matrix(1, 1) * y) + a.matrix(1, 2) * z) + a.matrix(1, 3);
        output.z[i] =
            ((a.matrix(2, 0) * x + a.matrix(2, 1) * y) + a.matrix(2, 2) * z) + a.matrix(2, 3);
    }
}
BENCH_NOINLINE bool library_soa(const affine3f &a, const_soa3<float> in, soa3<float> out) {
    return transform_points(a, in, out);
}
BENCH_NOINLINE bool library_aos(const affine3f &a, std::span<const vec3f> in,
                                std::span<vec3f> out) {
    return transform_points(a, in, out);
}
BENCH_NOINLINE bool library_dot(const_soa3<float> in, std::span<float> out) {
    return dot_batch(in, in, out);
}
BENCH_NOINLINE void library_normalize(std::span<const vec3f> in, std::span<vec3f> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = normalize(in[i]);
}
BENCH_NOINLINE void library_rotate(const quatf &q, std::span<const vec3f> in,
                                   std::span<vec3f> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = rotate(q, in[i]);
}
BENCH_NOINLINE void library_multiply(const mat4f &a, std::span<const mat4f> in,
                                     std::span<mat4f> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = a * in[i];
}
BENCH_NOINLINE bool library_inverse(std::span<const mat4f> in, std::span<mat4f> out) {
    for (std::size_t i = 0; i < in.size(); ++i) {
        const auto r = inverse(in[i]);
        if (!r)
            return false;
        out[i] = *r;
    }
    return true;
}

template <class F>
double measure(const char *label, std::size_t items, F &&work, double &checksum) {
    constexpr std::size_t samples = 9, rounds = 12;
    for (int warmup = 0; warmup < 3; ++warmup)
        checksum += work();
    std::array<double, samples> times{};
    for (std::size_t i = 0; i < samples; ++i) {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t j = 0; j < rounds; ++j)
            checksum += work();
        const auto elapsed = std::chrono::steady_clock::now() - start;
        times[i] =
            std::chrono::duration<double, std::nano>(elapsed).count() / double(items * rounds);
    }
    std::sort(times.begin(), times.end());
    const double median = times[samples / 2];
    std::printf("%-25s %10.3f ns/item  %10.3f Mitems/s  min=%8.3f max=%8.3f\n", label, median,
                1000.0 / median, times.front(), times.back());
    return median;
}
int main(int argc, char **argv) {
    std::size_t count = 65536;
    if (argc > 1) {
        const std::string_view text = argv[1];
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), count);
        if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || count < 16 ||
            count > (1U << 24)) {
            std::fprintf(stderr, "usage: chmath_bench [item_count: 16..16777216]\n");
            return 1;
        }
    }
#if defined(__clang__)
    std::printf("compiler: Clang %s\n", __clang_version__);
#elif defined(_MSC_VER)
    std::printf("compiler: MSVC %d\n", _MSC_VER);
#elif defined(__GNUC__)
    std::printf("compiler: GCC %s\n", __VERSION__);
#endif
#if !defined(NDEBUG)
    std::printf("NOTE: Debug build; configure Release for meaningful timings.\n");
#endif
    std::printf("backend=%s, points=%zu, samples=9, rounds=12, warmup=3\n", batch_backend(), count);
    std::printf(
        "sizeof(vec3f)=%zu align=%zu sizeof(mat4f)=%zu sizeof(quatf)=%zu sizeof(affine3f)=%zu\n",
        sizeof(vec3f), alignof(vec3f), sizeof(mat4f), sizeof(quatf), sizeof(affine3f));
    std::vector<float> x(count), y(count), z(count), ox(count), oy(count), oz(count), dots(count);
    std::vector<vec3f> packed(count), output(count);
    std::mt19937 random(20260908);
    std::uniform_real_distribution<float> dist(-10, 10);
    for (std::size_t i = 0; i < count; ++i) {
        packed[i] = {dist(random), dist(random), dist(random)};
        x[i] = packed[i][0];
        y[i] = packed[i][1];
        z[i] = packed[i][2];
    }
    const auto q = from_euler_xyz(vec3f{0.3f, 0.7f, -0.9f});
    const auto model = compose(trs<float>{{1, 2, 3}, q, {2, 3, 4}});
    const const_soa3<float> input{x, y, z};
    const soa3<float> out{ox, oy, oz};
    if (!library_soa(model, input, out) || !library_aos(model, packed, output))
        return 2;
    for (std::size_t i = 0; i < count; ++i)
        if (!almost_equal(vec3f{ox[i], oy[i], oz[i]}, output[i], 1e-5f, 1e-5f))
            return 3;
    double checksum{};
    const double baseline = measure(
        "baseline SoA transform", count,
        [&] {
            baseline_soa(model, input, out);
            return double(ox[count / 2]);
        },
        checksum);
    const double simd = measure(
        "chmath SoA transform", count,
        [&] {
            if (!library_soa(model, input, out))
                std::abort();
            return double(ox[count / 2]);
        },
        checksum);
    measure(
        "chmath AoS transform", count,
        [&] {
            if (!library_aos(model, packed, output))
                std::abort();
            return double(output[count / 2][0]);
        },
        checksum);
    measure(
        "chmath SoA dot", count,
        [&] {
            if (!library_dot(input, dots))
                std::abort();
            return double(dots[count / 2]);
        },
        checksum);
    measure(
        "chmath stable normalize", count,
        [&] {
            library_normalize(packed, output);
            return double(output[count / 2][0]);
        },
        checksum);
    measure(
        "chmath quaternion rotate", count,
        [&] {
            library_rotate(q, packed, output);
            return double(output[count / 2][0]);
        },
        checksum);
    const std::size_t matrices = std::min(count, std::size_t(4096));
    std::vector<mat4f> ms(matrices), mo(matrices);
    for (std::size_t i = 0; i < matrices; ++i)
        ms[i] = to_matrix(compose(trs<float>{packed[i], q, {1, 2, 3}}));
    measure(
        "chmath mat4 multiply", matrices,
        [&] {
            library_multiply(to_matrix(model), ms, mo);
            return double(mo[matrices / 2](0, 0));
        },
        checksum);
    measure(
        "chmath mat4 LU inverse", matrices,
        [&] {
            if (!library_inverse(ms, mo))
                std::abort();
            return double(mo[matrices / 2](0, 0));
        },
        checksum);
    std::printf("SoA speedup over compiler-optimized baseline: %.3fx; checksum=%.9f\n",
                baseline / simd, checksum);
    return is_finite(checksum) ? 0 : 4;
}
