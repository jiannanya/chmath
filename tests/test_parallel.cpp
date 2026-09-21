// Parallel execution tests. This executable is compiled with
// CHMATH_ENABLE_PARALLEL=1 (see CMakeLists.txt) so both the serial and the
// worker-spawning paths are exercised in every configuration, including
// sanitizer builds. Each parallel result is compared against the workers=1 run
// of the same operation.
#include "guarded_buffer.hpp"
#include "test.hpp"

#include <algorithm>
#include <atomic>
#include <chmath/chmath.hpp>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

using namespace chm;

CH_TEST(parallel_configuration_is_reported) {
    const unsigned hardware = parallel_hardware_threads();
    CHECK(hardware >= 1U);
    CHECK(parallel_default_workers() == hardware);
    CHECK(parallel_min_chunk >= 1U);
    // Short batches and explicit serial requests never split.
    CHECK(parallel_workers(0U) == 0U);
    CHECK(parallel_workers(parallel_min_chunk - 1U, 8U) == 1U);
    CHECK(parallel_workers(parallel_min_chunk * 8U, 1U) == 1U);
    // Requests are clamped to the hardware and to the available item budget.
    CHECK(parallel_workers(parallel_min_chunk * 8U, 1000U) <= hardware);
    if (hardware >= 2U) {
        CHECK(parallel_workers(parallel_min_chunk * 8U, 2U) == 2U);
        CHECK(parallel_workers(parallel_min_chunk * 2U + 1U, 64U) <= hardware);
    }
}

CH_TEST(parallel_split_respects_cache_lines) {
    using detail::aligned_split;
    const std::size_t line = detail::parallel_line_bytes;
    CHECK(aligned_split(0U, 4U, sizeof(float), 0U) == 0U);
    CHECK(aligned_split(1000U, 4U, sizeof(float), 0U) == 0U);
    CHECK(aligned_split(1000U, 4U, sizeof(float), 4U) == 1000U);
    CHECK(aligned_split(1000U, 4U, sizeof(float), 7U) == 1000U);
    if (line % sizeof(float) == 0U) {
        const std::size_t per_line = line / sizeof(float);
        for (unsigned index = 1U; index < 4U; ++index) {
            const std::size_t split = aligned_split(100000U, 4U, sizeof(float), index);
            CHECK(split % per_line == 0U);
        }
    }
    // Monotonic in the index, so every chunk end is also the next chunk begin.
    for (unsigned index = 0U; index + 1U <= 8U; ++index)
        CHECK(aligned_split(65537U, 8U, sizeof(double), index) <=
              aligned_split(65537U, 8U, sizeof(double), index + 1U));
}

CH_TEST(parallel_partition_covers_each_index_once) {
    const std::size_t counts[] = {0U,
                                  1U,
                                  parallel_min_chunk - 1U,
                                  parallel_min_chunk,
                                  parallel_min_chunk + 1U,
                                  parallel_min_chunk * 3U + 7U};
    const unsigned requests[] = {1U, 2U, 3U, 5U, 8U, 64U};
    for (std::size_t count : counts)
        for (unsigned requested : requests) {
            std::vector<unsigned> visits(count, 0U);
            std::atomic<bool> malformed{false};
            const unsigned active =
                parallel_for(count, requested, sizeof(unsigned),
                             [&](std::size_t begin, std::size_t end) {
                                 if (begin > end || end > count)
                                     malformed.store(true, std::memory_order_relaxed);
                                 for (std::size_t i = begin; i < end; ++i)
                                     ++visits[i];
                             });
            CHECK(!malformed.load());
            CHECK(active == parallel_workers(count, requested));
            std::size_t total = 0U;
            for (unsigned visit : visits) {
                CHECK(visit == 1U);
                total += visit;
            }
            CHECK(total == count);
        }
}

CH_TEST(parallel_matches_serial_values) {
    const std::size_t count = parallel_min_chunk * 2U + 123U;
    const auto fill = [](double *out, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i)
            out[i] = static_cast<double>(i) * 0.5 + 3.0;
    };
    std::vector<double> serial(count), parallel(count);
    parallel_for(count, 1U, sizeof(double),
                 [&](std::size_t begin, std::size_t end) { fill(serial.data(), begin, end); });
    parallel_for(count, 8U, sizeof(double),
                 [&](std::size_t begin, std::size_t end) { fill(parallel.data(), begin, end); });
    CHECK(parallel == serial);
}

CH_TEST(parallel_repeated_calls_are_stable) {
    const std::size_t count = parallel_min_chunk * 2U + 11U;
    std::vector<double> first(count), second(count);
    for (int repeat = 0; repeat < 3; ++repeat) {
        std::vector<double> &out = repeat == 0 ? first : second;
        std::fill(out.begin(), out.end(), -1.0);
        parallel_for(count, 4U, sizeof(double), [&](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i)
                out[i] = std::sin(static_cast<double>(i)) * 2.0 + static_cast<double>(i % 7U);
        });
        if (repeat > 0)
            CHECK(out == first);
    }
}

CH_TEST(parallel_zero_count_runs_no_work) {
    std::atomic<unsigned> calls{0U};
    const unsigned active =
        parallel_for(0U, 8U, sizeof(int), [&](std::size_t, std::size_t) { ++calls; });
    CHECK(active == 0U);
    CHECK(calls.load() == 0U);
}

CH_TEST(parallel_clamps_requested_workers) {
    const std::size_t count = parallel_min_chunk * 16U;
    std::atomic<unsigned> ranges{0U};
    const unsigned active = parallel_for(count, 1000U, sizeof(int), [&](std::size_t, std::size_t) {
        ++ranges;
    });
    CHECK(active == parallel_workers(count, 1000U));
    CHECK(active <= parallel_hardware_threads());
    CHECK(ranges.load() == active);
    CHECK(parallel_for(count, 0U, sizeof(int), [](std::size_t, std::size_t) {}) == 1U);
}

CH_TEST(parallel_propagates_body_exceptions) {
    if (parallel_workers(parallel_min_chunk * 4U, 4U) <= 1U)
        return; // Single-thread hardware; the worker path cannot be reached.
    const std::size_t count = parallel_min_chunk * 4U;
    bool caught_caller = false;
    try {
        parallel_for(count, 4U, sizeof(int), [](std::size_t begin, std::size_t) {
            if (begin == 0U)
                throw std::runtime_error{"caller range"};
            throw std::runtime_error{"worker range"};
        });
    } catch (const std::runtime_error &error) {
        caught_caller = std::string_view{error.what()} == "caller range";
    }
    CHECK(caught_caller);
    bool caught_worker = false;
    try {
        parallel_for(count, 4U, sizeof(int), [](std::size_t begin, std::size_t) {
            if (begin != 0U)
                throw std::runtime_error{"worker range"};
        });
    } catch (const std::runtime_error &error) {
        caught_worker = std::string_view{error.what()} == "worker range";
    }
    CHECK(caught_worker);
}

CH_TEST(parallel_batch_transform_matches_serial) {
    const std::size_t n = parallel_min_chunk * 2U + 7U;
    std::vector<float> x(n), y(n), z(n), px(n), py(n), pz(n), rx(n), ry(n), rz(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = test::sample<float>(-2.0F, 2.0F);
        y[i] = test::sample<float>(-2.0F, 2.0F);
        z[i] = test::sample<float>(-2.0F, 2.0F);
    }
    const auto a =
        compose(trs<float>{{1, 2, 3}, from_euler_xyz(vec3f{0.3F, 0.5F, 0.7F}), {2, 3, 4}});
    const const_soa3<float> input{x, y, z};
    const soa3<float> reference{rx, ry, rz};
    CHECK(transform_points(a, input, reference, 1U));
    CHECK(transform_points(a, input, soa3<float>{px, py, pz}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR((vec3f{px[i], py[i], pz[i]}), (vec3f{rx[i], ry[i], rz[i]}), 1e-5F);
    // Exact in-place is allowed and must match the out-of-place serial result.
    std::vector<float> ix(x), iy(y), iz(z);
    CHECK(transform_points(a, const_soa3<float>{ix, iy, iz}, soa3<float>{ix, iy, iz}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR((vec3f{ix[i], iy[i], iz[i]}), (vec3f{rx[i], ry[i], rz[i]}), 1e-5F);
    // Packed AoS overload shares the same contract.
    std::vector<vec3f> packed(n), aos_out(n), aos_ref(n);
    for (std::size_t i = 0; i < n; ++i)
        packed[i] = {x[i], y[i], z[i]};
    CHECK(transform_points(a, std::span<const vec3f>{packed}, std::span<vec3f>{aos_ref}, 1U));
    CHECK(transform_points(a, std::span<const vec3f>{packed}, std::span<vec3f>{aos_out}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR(aos_out[i], aos_ref[i], 1e-5F);
    CHECK(transform_vectors(a, std::span<const vec3f>{packed}, std::span<vec3f>{aos_out}, 8U));
    CHECK(transform_vectors(a, std::span<const vec3f>{packed}, std::span<vec3f>{aos_ref}, 1U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR(aos_out[i], aos_ref[i], 1e-5F);
}

CH_TEST(parallel_batch_rotation_and_normalization_match_serial) {
    const std::size_t n = parallel_min_chunk * 2U + 3U;
    std::vector<double> x(n), y(n), z(n), px(n), py(n), pz(n), rx(n), ry(n), rz(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = test::sample<double>(-2.0, 2.0);
        y[i] = test::sample<double>(-2.0, 2.0);
        z[i] = test::sample<double>(-2.0, 2.0);
    }
    const auto q = *from_axis_angle(vec3d{1, 2, 3}, 0.7);
    const const_soa3<double> input{x, y, z};
    CHECK(rotate_vectors(q, input, soa3<double>{rx, ry, rz}, 1U));
    CHECK(rotate_vectors(q, input, soa3<double>{px, py, pz}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR((vec3d{px[i], py[i], pz[i]}), (vec3d{rx[i], ry[i], rz[i]}), 1e-11);
    CHECK(normalize_vectors(input, soa3<double>{px, py, pz}, 8U));
    CHECK(normalize_vectors(input, soa3<double>{rx, ry, rz}, 1U));
    for (std::size_t i = 0; i < n; ++i) {
        NEAR((vec3d{px[i], py[i], pz[i]}), (vec3d{rx[i], ry[i], rz[i]}), 1e-11);
        NEAR(length(vec3d{px[i], py[i], pz[i]}), 1.0, 1e-11);
    }
    std::vector<vec3d> packed(n), aos_out(n), aos_ref(n);
    for (std::size_t i = 0; i < n; ++i)
        packed[i] = {x[i], y[i], z[i]};
    CHECK(normalize_vectors(std::span<const vec3d>{packed}, std::span<vec3d>{aos_ref}, 1U));
    CHECK(normalize_vectors(std::span<const vec3d>{packed}, std::span<vec3d>{aos_out}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR(aos_out[i], aos_ref[i], 1e-11);
}

CH_TEST(parallel_batch_cross_and_dot_match_serial) {
    const std::size_t n = parallel_min_chunk * 2U + 5U;
    std::vector<float> ax(n), ay(n), az(n), bx(n), by(n), bz(n);
    for (std::size_t i = 0; i < n; ++i) {
        ax[i] = test::sample<float>(-3.0F, 3.0F);
        ay[i] = test::sample<float>(-3.0F, 3.0F);
        az[i] = test::sample<float>(-3.0F, 3.0F);
        bx[i] = test::sample<float>(-3.0F, 3.0F);
        by[i] = test::sample<float>(-3.0F, 3.0F);
        bz[i] = test::sample<float>(-3.0F, 3.0F);
    }
    const const_soa3<float> a{ax, ay, az}, b{bx, by, bz};
    std::vector<float> dots_ref(n), dots_out(n);
    CHECK(dot_batch(a, b, std::span<float>{dots_ref}, 1U));
    CHECK(dot_batch(a, b, std::span<float>{dots_out}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR(dots_out[i], dots_ref[i], 1e-5F);
    std::vector<float> cx(n), cy(n), cz(n), dx(n), dy(n), dz(n);
    CHECK(cross_batch(a, b, soa3<float>{cx, cy, cz}, 1U));
    CHECK(cross_batch(a, b, soa3<float>{dx, dy, dz}, 8U));
    for (std::size_t i = 0; i < n; ++i)
        NEAR((vec3f{dx[i], dy[i], dz[i]}), (vec3f{cx[i], cy[i], cz[i]}), 1e-5F);
}

CH_TEST(parallel_batch_validation_precedes_writes) {
    constexpr float sentinel = -9876.5F;
    const std::size_t n = parallel_min_chunk * 2U;
    const auto a = translation(vec3f{1, 2, 3});
    std::vector<float> x(n), y(n), z(n), ox(n, sentinel), oy(n, sentinel), oz(n, sentinel);
    std::vector<float> short_output(n - 1U, sentinel);
    const const_soa3<float> input{x, y, z};
    // Mismatched channel lengths are rejected before any worker starts.
    CHECK(!transform_points(a, input, soa3<float>{short_output, oy, oz}, 8U));
    CHECK(!transform_points(a, input, soa3<float>{ox, short_output, oz}, 8U));
    CHECK(!dot_batch(input, input, std::span<float>{short_output}, 8U));
    // Partially overlapping output channels are rejected as well.
    CHECK(!transform_points(a, input, soa3<float>{ox, std::span<float>{ox}.subspan(1U), oz}, 8U));
    // A valid channel layout that partially overlaps the input is rejected too.
    std::vector<float> wide(n + 1U, sentinel);
    const const_soa3<float> overlapping_input{std::span<float>{wide}.subspan(0U, n), y, z};
    const soa3<float> overlapping_output{std::span<float>{wide}.subspan(1U, n), oy, oz};
    CHECK(!transform_points(a, overlapping_input, overlapping_output, 8U));
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(ox[i] == sentinel);
        CHECK(oy[i] == sentinel);
        CHECK(oz[i] == sentinel);
        if (i + 1U < n)
            CHECK(short_output[i] == sentinel);
    }
    // Queries share the same validate-first contract.
    const auto prepared = *prepared_ray<float>::create(ray3f{{0, 0, 4}, {0, 0, -1}});
    std::vector<aabb3f> boxes(n);
    std::vector<std::optional<ray_interval<float>>> short_hits(n - 1U);
    CHECK(!intersect_many(prepared, std::span<const aabb3f>{boxes},
                          std::span<std::optional<ray_interval<float>>>{short_hits}, 0.0F,
                          std::numeric_limits<float>::infinity(), 8U));
    const auto frustum_view = *extract_frustum(*perspective(1.0F, 1.5F, 0.1F, 100.0F));
    std::vector<containment> short_visible(n - 1U);
    CHECK(!classify_batch(frustum_view, std::span<const aabb3f>{boxes},
                          std::span<containment>{short_visible}, 8U));
}

CH_TEST(parallel_queries_match_serial) {
    const std::size_t n = parallel_min_chunk * 2U + 9U;
    const ray3d ray{{-3, -2, 3}, {0.125, 0.25, -1}};
    const auto prepared = *prepared_ray<double>::create(ray);
    std::vector<aabb3d> boxes(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double lo0 = static_cast<double>(i % 17U) - 8.0;
        const double lo1 = static_cast<double>(i % 11U) - 5.0;
        const double lo2 = static_cast<double>(i % 7U) - 3.0;
        boxes[i] = {vec3d{lo0, lo1, lo2}, vec3d{lo0 + 1.0, lo1 + 1.0, lo2 + 1.0}};
    }
    std::vector<std::optional<ray_interval<double>>> reference(n), parallel(n);
    CHECK(intersect_many(prepared, std::span<const aabb3d>{boxes},
                         std::span<std::optional<ray_interval<double>>>{reference}, 0.0,
                         std::numeric_limits<double>::infinity(), 1U));
    CHECK(intersect_many(prepared, std::span<const aabb3d>{boxes},
                         std::span<std::optional<ray_interval<double>>>{parallel}, 0.0,
                         std::numeric_limits<double>::infinity(), 8U));
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(bool(parallel[i]) == bool(reference[i]));
        if (reference[i]) {
            NEAR(parallel[i]->enter, reference[i]->enter, 1e-11);
            NEAR(parallel[i]->exit, reference[i]->exit, 1e-11);
        }
    }
    const auto frustum_view = *extract_frustum(*perspective(1.0, 1.5, 0.1, 100.0));
    std::vector<containment> classified(n), classified_parallel(n);
    CHECK(classify_batch(frustum_view, std::span<const aabb3d>{boxes},
                         std::span<containment>{classified}, 1U));
    CHECK(classify_batch(frustum_view, std::span<const aabb3d>{boxes},
                         std::span<containment>{classified_parallel}, 8U));
    CHECK(classified == classified_parallel);
}

CH_TEST(parallel_rotation_rejects_nonunit_quaternion_before_workers) {
    const std::size_t n = parallel_min_chunk * 2U;
    std::vector<float> x(n), y(n), z(n), ox(n, 7.0F), oy(n, 7.0F), oz(n, 7.0F);
    const const_soa3<float> input{x, y, z};
    const soa3<float> output{ox, oy, oz};
    CHECK(!rotate_vectors(quatf{0.0F, 0.0F, 0.0F, 2.0F}, input, output, 8U));
    CHECK(!rotate_vectors(quatf{std::numeric_limits<float>::quiet_NaN(), 0.0F, 0.0F, 1.0F}, input,
                          output, 8U));
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(ox[i] == 7.0F);
        CHECK(oy[i] == 7.0F);
        CHECK(oz[i] == 7.0F);
    }
}

CH_TEST(safety_parallel_batch_guarded_boundaries) {
    const std::size_t n = parallel_min_chunk * 2U + 3U;
    const auto a =
        compose(trs<double>{{1, 2, 3}, from_euler_xyz(vec3d{0.3, 0.5, 0.7}), {2, 3, 4}});
    std::vector<double> bx(n), by(n), bz(n), expected_x(n), expected_y(n), expected_z(n);
    for (std::size_t i = 0; i < n; ++i) {
        bx[i] = static_cast<double>(i % 29U) * 0.25 - 3.0;
        by[i] = static_cast<double>(i % 13U) * -0.5 + 1.0;
        bz[i] = static_cast<double>(i % 7U) - 2.0;
    }
    const const_soa3<double> plain{bx, by, bz};
    CHECK(transform_points(a, plain, soa3<double>{expected_x, expected_y, expected_z}, 1U));
    for (bool at_end : {false, true}) {
        safety::guarded_buffer<double> x(n, at_end), y(n, at_end), z(n, at_end), ox(n, at_end),
            oy(n, at_end), oz(n, at_end);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = bx[i];
            y[i] = by[i];
            z[i] = bz[i];
        }
        x.read_only();
        y.read_only();
        z.read_only();
        CHECK(transform_points(a,
                               const_soa3<double>{x.const_span(), y.const_span(), z.const_span()},
                               soa3<double>{ox.span(), oy.span(), oz.span()}, 8U));
        for (std::size_t i = 0; i < n; ++i) {
            NEAR(ox[i], expected_x[i], 1e-12);
            NEAR(oy[i], expected_y[i], 1e-12);
            NEAR(oz[i], expected_z[i], 1e-12);
        }
    }
    // The query batch clones the same guard-page contract.
    std::vector<aabb3d> boxes_plain(n);
    for (std::size_t i = 0; i < n; ++i)
        boxes_plain[i] = {vec3d{bx[i] - 1.0, by[i] - 1.0, bz[i] - 1.0},
                          vec3d{bx[i] + 1.0, by[i] + 1.0, bz[i] + 1.0}};
    const ray3d ray{{0, 0, 30}, {0.1, 0.2, -1}};
    const auto prepared = *prepared_ray<double>::create(ray);
    std::vector<std::optional<ray_interval<double>>> reference(n);
    CHECK(intersect_many(prepared, std::span<const aabb3d>{boxes_plain},
                         std::span<std::optional<ray_interval<double>>>{reference}, 0.0,
                         std::numeric_limits<double>::infinity(), 1U));
    safety::guarded_buffer<aabb3d> guarded_boxes(n);
    safety::guarded_buffer<std::optional<ray_interval<double>>> guarded_hits(n);
    for (std::size_t i = 0; i < n; ++i)
        guarded_boxes[i] = boxes_plain[i];
    CHECK(intersect_many(prepared, std::span<const aabb3d>{guarded_boxes.const_span()},
                         std::span<std::optional<ray_interval<double>>>{guarded_hits.span()}, 0.0,
                         std::numeric_limits<double>::infinity(), 8U));
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(bool(guarded_hits[i]) == bool(reference[i]));
        if (reference[i]) {
            NEAR(guarded_hits[i]->enter, reference[i]->enter, 1e-11);
            NEAR(guarded_hits[i]->exit, reference[i]->exit, 1e-11);
        }
    }
}
