#include <atomic>
#include <charconv>
#include <chmath/chmath.hpp>
#include <chrono>
#include <cstdio>
#include <random>
#include <string_view>
#include <vector>
#if defined(_MSC_VER)
#define CH_BENCH_NOINLINE __declspec(noinline)
#else
#define CH_BENCH_NOINLINE __attribute__((noinline))
#endif
using namespace chm;
template <class F> CH_BENCH_NOINLINE double invoke_work(F &f) {
    return f();
}
template <class F>
void measure_extended(const char *name, std::size_t n, F work, double &checksum) {
    for (int i = 0; i < 3; ++i)
        checksum += invoke_work(work);
    std::array<double, 9> times{};
    for (auto &sample : times) {
        const auto start = std::chrono::steady_clock::now();
        for (int j = 0; j < 12; ++j) {
            std::atomic_signal_fence(std::memory_order_seq_cst);
            checksum += invoke_work(work);
        }
        sample = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start)
                     .count() /
                 double(n * 12);
    }
    std::sort(times.begin(), times.end());
    std::printf("%-30s %10.3f ns/item  min=%10.3f max=%10.3f\n", name, times[4], times[0],
                times[8]);
}
template <class T, std::size_t N>
void matrix_benchmark(const char *name, std::size_t count, bool invert, double &checksum) {
    std::vector<mat<T, N, N>> in(count), out(count);
    for (std::size_t p = 0; p < count; ++p)
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                in[p](i, j) = i == j ? T(N) + T(p % 7) * T(0.01) : T(1) / T(i + j + 2);
    const auto a = in[0];
    if (invert)
        measure_extended(
            name, count,
            [&] {
                for (std::size_t i = 0; i < count; ++i) {
                    const auto r = inverse(in[i]);
                    if (!r)
                        std::abort();
                    out[i] = *r;
                }
                return double(out[count / 2](0, 0));
            },
            checksum);
    else
        measure_extended(
            name, count,
            [&] {
                for (std::size_t i = 0; i < count; ++i)
                    out[i] = a * in[i];
                return double(out[count / 2](0, 0));
            },
            checksum);
    if (!is_finite(out[count / 2]))
        std::abort();
}
int main(int argc, char **argv) {
    std::size_t n = 65536;
    if (argc > 1) {
        const std::string_view arg = argv[1];
        const auto parsed = std::from_chars(arg.data(), arg.data() + arg.size(), n);
        if (parsed.ec != std::errc{} || parsed.ptr != arg.data() + arg.size() || n < 16 ||
            n > (1U << 24))
            return 1;
    }
    std::printf("extended benchmark backend=%s points=%zu samples=9 rounds=12 warmup=3\n",
                batch_backend(), n);
    double checksum{};
    const auto small = std::min(n, std::size_t(2048));
    matrix_benchmark<float, 2>("mat2f multiply", small, false, checksum);
    matrix_benchmark<float, 3>("mat3f multiply", small, false, checksum);
    matrix_benchmark<double, 4>("mat4d multiply", small, false, checksum);
    matrix_benchmark<double, 3>("mat3d inverse", small, true, checksum);
    matrix_benchmark<double, 4>("mat4d inverse", small, true, checksum);
    std::vector<vec3d> input(n), output(n);
    std::vector<double> x(n), y(n), z(n), ox(n), oy(n), oz(n), dots(n);
    std::mt19937 rng(20260912);
    std::uniform_real_distribution<double> random(-2, 2);
    for (std::size_t i = 0; i < n; ++i) {
        input[i] = {random(rng), random(rng), random(rng)};
        x[i] = input[i][0];
        y[i] = input[i][1];
        z[i] = input[i][2];
    }
    measure_extended(
        "normalize vec3d", n,
        [&] {
            for (std::size_t i = 0; i < n; ++i)
                output[i] = normalize(input[i]);
            return output[n / 2][0];
        },
        checksum);
    std::vector<quatd> qa(small), qb(small), qo(small);
    std::vector<mat3d> rotations(small);
    for (std::size_t i = 0; i < small; ++i) {
        qa[i] = from_euler_xyz(input[i]);
        qb[i] = from_euler_xyz(input[i] + vec3d{0.001});
        rotations[i] = to_matrix(qa[i]);
    }
    measure_extended(
        "quaternion nlerp", small,
        [&] {
            for (std::size_t i = 0; i < small; ++i)
                qo[i] = nlerp(qa[i], qb[i], 0.3);
            return qo[small / 2].w();
        },
        checksum);
    measure_extended(
        "quaternion slerp near", small,
        [&] {
            for (std::size_t i = 0; i < small; ++i)
                qo[i] = slerp(qa[i], qb[i], 0.3);
            return qo[small / 2].w();
        },
        checksum);
    for (std::size_t i = 0; i < small; ++i)
        qb[i] = from_euler_xyz(input[i] + vec3d{0.7});
    measure_extended(
        "quaternion slerp wide", small,
        [&] {
            for (std::size_t i = 0; i < small; ++i)
                qo[i] = slerp(qa[i], qb[i], 0.3);
            return qo[small / 2].w();
        },
        checksum);
    measure_extended(
        "quaternion from matrix", small,
        [&] {
            for (std::size_t i = 0; i < small; ++i) {
                const auto r = from_matrix(rotations[i]);
                if (!r)
                    std::abort();
                qo[i] = *r;
            }
            return qo[small / 2].w();
        },
        checksum);
    measure_extended(
        "rotation between", small,
        [&] {
            for (std::size_t i = 0; i < small; ++i) {
                const auto r = rotation_between(input[i], vec3d{1, 2, 3});
                if (!r)
                    std::abort();
                qo[i] = *r;
            }
            return qo[small / 2].w();
        },
        checksum);
    const auto a = compose(trs<double>{{1, 2, 3}, qa[0], {2, 3, 4}});
    const const_soa3<double> in{x, y, z};
    const soa3<double> out{ox, oy, oz};
    measure_extended(
        "SoA double transform", n,
        [&] {
            if (!transform_points(a, in, out))
                std::abort();
            return ox[n / 2];
        },
        checksum);
    measure_extended(
        "SoA double dot", n,
        [&] {
            if (!dot_batch(in, in, std::span<double>{dots}))
                std::abort();
            return dots[n / 2];
        },
        checksum);
    const bezier<double, 3, 3> curve{
        {vec3d{0, 0, 0}, vec3d{1, 2, 0}, vec3d{2, 2, 1}, vec3d{3, 0, 1}}};
    measure_extended(
        "Bezier value and derivative", small,
        [&] {
            for (std::size_t i = 0; i < small; ++i) {
                const double t = double(i) / double(small);
                output[i] = curve.evaluate(t) + curve.derivative().evaluate(t);
            }
            return output[small / 2][0];
        },
        checksum);
    const ray3d ray{{0, 0, 10}, {0.2, 0.3, -1}};
    std::vector<aabb3d> boxes(n);
    std::vector<ray_interval<double>> hits(n);
    for (std::size_t i = 0; i < n; ++i)
        boxes[i] = {input[i] - vec3d{1}, input[i] + vec3d{1}};
    measure_extended(
        "ray AABB queries", n,
        [&] {
            for (std::size_t i = 0; i < n; ++i)
                hits[i] = intersect(ray, boxes[i]).value_or(ray_interval<double>{-1, -1});
            return hits[n / 2].enter;
        },
        checksum);
    std::printf("checksum=%.12g\n", checksum);
    return is_finite(checksum) ? 0 : 2;
}
