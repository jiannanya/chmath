#include "boundary_helpers.hpp"
#include "guarded_buffer.hpp"
#include "test.hpp"
#include <bit>
#include <cstring>
#include <vector>
using namespace chm;

namespace safety {
inline affine3f transform() {
    affine3f a;
    a.matrix = mat<float, 3, 4>::from_rows(
        {vec4f{2, -1, 0.5f, 7}, vec4f{-3, 4, 2, -1}, vec4f{1, 2, -2, 3}});
    return a;
}
inline void soa_guard(std::size_t n, bool tail = true, bool readonly = false, int inplace = 0) {
    guarded_buffer<float> x(n, tail), y(n, tail), z(n, tail), ox(n, tail), oy(n, tail), oz(n, tail);
    std::vector<vec3f> reference(n);
    const auto a = transform();
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = float(i) + 0.25f;
        y[i] = -float(i) - 2;
        z[i] = float(i) * 2 + 3;
        reference[i] = transform_point(a, vec3f{x[i], y[i], z[i]});
    }
    if (readonly) {
        x.read_only();
        y.read_only();
        z.read_only();
    }
    const_soa3<float> input{x.const_span(), y.const_span(), z.const_span()};
    soa3<float> output{ox.span(), oy.span(), oz.span()};
    if (inplace == 1)
        output = {x.span(), y.span(), z.span()};
    if (inplace == 2)
        output = {y.span(), z.span(), x.span()};
    CHECK(transform_points(a, input, output));
    for (std::size_t i = 0; i < n; ++i)
        NEAR((vec3f{output.x[i], output.y[i], output.z[i]}), reference[i], 1e-6f);
}
inline void aos_guard(std::size_t n, bool tail) {
    guarded_buffer<vec3f> input(n, tail), output(n, tail);
    for (std::size_t i = 0; i < n; ++i)
        input[i] = {float(i), 2, -3};
    input.read_only();
    const auto a = transform();
    CHECK(transform_points(a, input.const_span(), output.span()));
    for (std::size_t i = 0; i < n; ++i)
        NEAR(output[i], transform_point(a, input[i]), 1e-6f);
}
inline void dot_guard(std::size_t n, bool tail, bool inplace = false) {
    guarded_buffer<float> x(n, tail), y(n, tail), z(n, tail), ox(n, tail);
    std::vector<float> expected(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = float(i) + 1;
        y[i] = -2;
        z[i] = 3;
        expected[i] = x[i] * x[i] + 13;
    }
    const_soa3<float> input{x.const_span(), y.const_span(), z.const_span()};
    const auto out = inplace ? x.span() : ox.span();
    CHECK(dot_batch(input, input, out));
    for (std::size_t i = 0; i < n; ++i)
        CHECK(out[i] == expected[i]);
}
} // namespace safety
CH_TEST(safety_soa_head_one) {
    safety::soa_guard(1, false);
}
CH_TEST(safety_soa_tail_one) {
    safety::soa_guard(1);
}
CH_TEST(safety_soa_tail_two) {
    safety::soa_guard(2);
}
CH_TEST(safety_soa_tail_three) {
    safety::soa_guard(3);
}
CH_TEST(safety_soa_tail_four) {
    safety::soa_guard(4);
}
CH_TEST(safety_soa_tail_five) {
    safety::soa_guard(5);
}
CH_TEST(safety_soa_tail_seven) {
    safety::soa_guard(7);
}
CH_TEST(safety_soa_tail_eight) {
    safety::soa_guard(8);
}
CH_TEST(safety_soa_tail_seventeen) {
    safety::soa_guard(17);
}
CH_TEST(safety_soa_head_seventeen) {
    safety::soa_guard(17, false);
}
CH_TEST(safety_aos_tail_one) {
    safety::aos_guard(1, true);
}
CH_TEST(safety_aos_tail_five) {
    safety::aos_guard(5, true);
}
CH_TEST(safety_aos_head_five) {
    safety::aos_guard(5, false);
}
CH_TEST(safety_dot_tail_three) {
    safety::dot_guard(3, true);
}
CH_TEST(safety_dot_tail_five) {
    safety::dot_guard(5, true);
}
CH_TEST(safety_dot_head_nine) {
    safety::dot_guard(9, false);
}
CH_TEST(safety_soa_read_only_inputs) {
    safety::soa_guard(9, true, true);
}
CH_TEST(safety_soa_inplace_tail_seven) {
    safety::soa_guard(7, true, false, 1);
}
CH_TEST(safety_soa_cyclic_channels_tail_five) {
    safety::soa_guard(5, true, false, 2);
}
CH_TEST(safety_dot_inplace_tail_seven) {
    safety::dot_guard(7, true, true);
}
CH_TEST(safety_aos_forward_overlap_transactional) {
    std::array<vec3f, 9> buffer{};
    for (std::size_t i = 0; i < buffer.size(); ++i)
        buffer[i] = vec3f(float(i));
    const auto before = buffer;
    CHECK(!transform_points(safety::transform(), std::span<const vec3f>{buffer}.first(8),
                            std::span<vec3f>{buffer}.subspan(1)));
    CHECK(buffer == before);
}
CH_TEST(safety_aos_backward_overlap_transactional) {
    std::array<vec3f, 9> buffer{};
    for (std::size_t i = 0; i < buffer.size(); ++i)
        buffer[i] = vec3f(float(i));
    const auto before = buffer;
    CHECK(!transform_points(safety::transform(), std::span<const vec3f>{buffer}.subspan(1),
                            std::span<vec3f>{buffer}.first(8)));
    CHECK(buffer == before);
}
CH_TEST(safety_aos_size_mismatch_transactional) {
    const std::array<vec3f, 5> input{};
    std::array<vec3f, 4> out{};
    out.fill(vec3f{99});
    CHECK(!transform_points(safety::transform(), std::span<const vec3f>{input},
                            std::span<vec3f>{out}));
    for (auto v : out)
        CHECK(v == vec3f{99});
}
CH_TEST(safety_soa_output_channel_alias_rejected) {
    std::array<float, 24> input{}, out{};
    out.fill(77);
    const auto before = out;
    const_soa3<float> in{std::span<const float>{input}.subspan(0, 8),
                         std::span<const float>{input}.subspan(8, 8),
                         std::span<const float>{input}.subspan(16, 8)};
    for (std::size_t a = 0; a < 3; ++a)
        for (std::size_t b = a + 1; b < 3; ++b) {
            std::array<std::span<float>, 3> channels{std::span<float>{out}.subspan(0, 8),
                                                     std::span<float>{out}.subspan(8, 8),
                                                     std::span<float>{out}.subspan(16, 8)};
            channels[b] = channels[a];
            CHECK(!transform_points(safety::transform(), in,
                                    soa3<float>{channels[0], channels[1], channels[2]}));
            CHECK(out == before);
        }
}
CH_TEST(safety_soa_partial_overlap_all_channel_pairs) {
    std::array<float, 96> storage{};
    storage.fill(77);
    const auto before = storage;
    for (std::size_t src = 0; src < 3; ++src)
        for (std::size_t dst = 0; dst < 3; ++dst)
            for (int direction : {-1, 1}) {
                std::array<std::span<const float>, 3> in{
                    std::span<const float>{storage}.subspan(1, 8),
                    std::span<const float>{storage}.subspan(17, 8),
                    std::span<const float>{storage}.subspan(33, 8)};
                std::array<std::span<float>, 3> out{std::span<float>{storage}.subspan(48, 8),
                                                    std::span<float>{storage}.subspan(64, 8),
                                                    std::span<float>{storage}.subspan(80, 8)};
                out[dst] = {storage.data() + static_cast<std::ptrdiff_t>(1 + src * 16) + direction,
                            8};
                CHECK(!transform_points(safety::transform(), const_soa3<float>{in[0], in[1], in[2]},
                                        soa3<float>{out[0], out[1], out[2]}));
                CHECK(storage == before);
            }
}
CH_TEST(safety_soa_length_mismatch_each_channel) {
    std::array<float, 48> storage{};
    storage.fill(77);
    const auto before = storage;
    for (std::size_t bad = 0; bad < 6; ++bad) {
        std::array<std::span<float>, 6> c{};
        for (std::size_t i = 0; i < 6; ++i)
            c[i] = std::span<float>{storage}.subspan(i * 8, 8);
        c[bad] = c[bad].first(7);
        CHECK(!transform_points(safety::transform(), const_soa3<float>{c[0], c[1], c[2]},
                                soa3<float>{c[3], c[4], c[5]}));
        CHECK(storage == before);
    }
}
CH_TEST(safety_dot_partial_overlap_each_input) {
    std::array<float, 96> storage{};
    storage.fill(77);
    const auto before = storage;
    std::array<std::span<const float>, 6> c{};
    for (std::size_t i = 0; i < 6; ++i)
        c[i] = std::span<const float>{storage}.subspan(i * 16 + 1, 8);
    for (std::size_t i = 0; i < 6; ++i)
        for (int direction : {-1, 1}) {
            const std::span<float> out{
                storage.data() + static_cast<std::ptrdiff_t>(i * 16 + 1) + direction, 8};
            CHECK(!dot_batch(const_soa3<float>{c[0], c[1], c[2]},
                             const_soa3<float>{c[3], c[4], c[5]}, out));
            CHECK(storage == before);
        }
}
CH_TEST(safety_dot_length_mismatch_each_argument) {
    std::array<float, 56> storage{};
    storage.fill(77);
    const auto before = storage;
    for (std::size_t bad = 0; bad < 7; ++bad) {
        std::array<std::span<float>, 7> c{};
        for (std::size_t i = 0; i < 7; ++i)
            c[i] = std::span<float>{storage}.subspan(i * 8, 8);
        c[bad] = c[bad].first(7);
        CHECK(!dot_batch(const_soa3<float>{c[0], c[1], c[2]}, const_soa3<float>{c[3], c[4], c[5]},
                         c[6]));
        CHECK(storage == before);
    }
}
CH_TEST(safety_empty_null_buffers) {
    CHECK(transform_points(affine3f{}, std::span<const vec3f>{}, std::span<vec3f>{}));
    CHECK(transform_points(affine3f{}, const_soa3<float>{}, soa3<float>{}));
    CHECK(dot_batch(const_soa3<float>{}, const_soa3<float>{}, std::span<float>{}));
}
CH_TEST(safety_adjacent_spans_do_not_overlap) {
    std::array<vec3f, 8> buffer{};
    for (std::size_t i = 0; i < 4; ++i)
        buffer[i] = vec3f(float(i));
    CHECK(transform_points(safety::transform(), std::span<const vec3f>{buffer}.first(4),
                           std::span<vec3f>{buffer}.subspan(4)));
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(buffer[i + 4] == transform_point(safety::transform(), buffer[i]));
}
CH_TEST(safety_packed_layout_and_array_stride) {
    static_assert(std::is_standard_layout_v<vec3f> && std::is_trivially_copyable_v<mat4f>);
    CHECK(sizeof(vec3f) == 12);
    CHECK(sizeof(quatf) == 16);
    CHECK(sizeof(mat4f) == 64);
    CHECK(sizeof(affine3f) == 48);
    const std::array<vec3f, 2> p{};
    CHECK(reinterpret_cast<const std::byte *>(&p[1]) - reinterpret_cast<const std::byte *>(&p[0]) ==
          12);
    CHECK(alignof(vec3f) == alignof(float));
    CHECK(alignof(mat4f) == alignof(float));
}
CH_TEST(safety_matrix_natural_alignment_without_simd_alignment) {
    constexpr auto identity_squared = mat4f::identity() * mat4f::identity();
    static_assert(identity_squared == mat4f::identity());
    struct alignas(16) packed {
        float prefix{};
        mat4f matrix{};
        std::array<float, 3> suffix{};
    };
    packed a;
    a.matrix = mat4f::identity();
    a.matrix(0, 3) = 4;
    CHECK(reinterpret_cast<std::uintptr_t>(&a.matrix) % 16 == 4);
    const auto m = a.matrix * a.matrix;
    CHECK(m(0, 3) == 8);
    CHECK(m(3, 3) == 1);
}
CH_TEST(safety_vector_representation_roundtrip) {
    const vec4f v{std::numeric_limits<float>::infinity(), -0.0f, 1.25f,
                  std::numeric_limits<float>::denorm_min()};
    const auto bytes = std::bit_cast<std::array<std::byte, sizeof(vec4f)>>(v);
    const auto r = std::bit_cast<vec4f>(bytes);
    CHECK(std::memcmp(&v, &r, sizeof(v)) == 0);
    CHECK(std::signbit(r[1]));
}
CH_TEST(safety_bspline_max_degree_guarded_arrays) {
    safety::guarded_buffer<vec2d> c(33);
    safety::guarded_buffer<double> k(66);
    for (std::size_t i = 0; i < 33; ++i) {
        c[i] = {double(i), 0};
        k[i] = 0;
        k[i + 33] = 1;
    }
    c.read_only();
    k.read_only();
    const auto s = bspline_view<double, 2, 32>::create(c.const_span(), k.const_span());
    CHECK(s);
    for (double t : {0.0, 0.5, 1.0}) {
        const auto v = s->evaluate(t), d = s->derivative(t);
        CHECK(v && d);
        NEAR(*v, vec2d{32 * t, 0}, 1e-12);
        NEAR(*d, vec2d{32, 0}, 1e-12);
    }
}
CH_TEST(safety_nurbs_degree_zero_guarded_arrays) {
    safety::guarded_buffer<vec2d> c(1);
    safety::guarded_buffer<double> k(2), w(1);
    c[0] = {2, 3};
    k[0] = 0;
    k[1] = 1;
    w[0] = 1e300;
    c.read_only();
    k.read_only();
    w.read_only();
    const auto s = nurbs_view<double, 2, 0>::create(c.const_span(), w.const_span(), k.const_span());
    CHECK(s);
    for (double t : {0.0, 0.5, 1.0})
        CHECK(s->evaluate(t) == vec2d{2, 3});
}
