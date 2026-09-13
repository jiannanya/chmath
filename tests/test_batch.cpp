#include "test.hpp"
#include <chmath/simd/batch.hpp>
#include <vector>

using namespace chm;

template <class T> void batch_properties() {
    const auto a = compose(trs<T>{{T(1), T(2), T(3)},
                                  from_euler_xyz(vec<T, 3>{T(0.3), T(0.7), T(-0.9)}),
                                  {T(2), T(-3), T(4)}});
    // Every short tail, zero count, deliberately unaligned offsets, and a large batch.
    for (std::size_t n : {0U, 1U, 2U, 3U, 4U, 5U, 7U, 8U, 15U, 16U, 17U, 255U, 1025U}) {
        std::vector<T> x(n + 2, T(-99)), y = x, z = x, ox = x, oy = x, oz = x, dots = x;
        std::vector<vec<T, 3>> packed(n), reference(n), output(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i + 1] = test::sample<T>();
            y[i + 1] = test::sample<T>();
            z[i + 1] = test::sample<T>();
            packed[i] = {x[i + 1], y[i + 1], z[i + 1]};
            reference[i] = transform_point(a, packed[i]);
        }
        const const_soa3<T> input{std::span<const T>{x}.subspan(1, n),
                                  std::span<const T>{y}.subspan(1, n),
                                  std::span<const T>{z}.subspan(1, n)};
        const soa3<T> out{std::span<T>{ox}.subspan(1, n), std::span<T>{oy}.subspan(1, n),
                          std::span<T>{oz}.subspan(1, n)};
        CHECK(transform_points(a, input, out));
        CHECK(dot_batch(input, input, std::span<T>{dots}.subspan(1, n)));
        CHECK(transform_points<T>(a, std::span<const vec<T, 3>>{packed},
                                  std::span<vec<T, 3>>{output}));
        for (std::size_t i = 0; i < n; ++i) {
            NEAR(vec<T, 3>{ox[i + 1], oy[i + 1], oz[i + 1]}, reference[i], T(1e-5));
            NEAR(output[i], reference[i], T(1e-5));
            NEAR(dots[i + 1], dot(packed[i], packed[i]), T(1e-5));
        }
        CHECK(ox.front() == T(-99) && ox.back() == T(-99) && oy.front() == T(-99) &&
              oy.back() == T(-99) && oz.front() == T(-99) && oz.back() == T(-99));
        CHECK(dots.front() == T(-99) && dots.back() == T(-99));
        const soa3<T> inplace{std::span<T>{x}.subspan(1, n), std::span<T>{y}.subspan(1, n),
                              std::span<T>{z}.subspan(1, n)};
        CHECK(transform_points(a, inplace.as_const(), inplace));
        CHECK(transform_points<T>(a, std::span<const vec<T, 3>>{packed},
                                  std::span<vec<T, 3>>{packed}));
        for (std::size_t i = 0; i < n; ++i) {
            NEAR(vec<T, 3>{x[i + 1], y[i + 1], z[i + 1]}, reference[i], T(1e-5));
            NEAR(packed[i], reference[i], T(1e-5));
        }
    }
}
CH_TEST(batch_scalar_equivalence_and_memory_boundaries) {
    batch_properties<float>();
    batch_properties<double>();
}
CH_TEST(batch_rejects_invalid_without_writing) {
    std::array<float, 16> x{}, y{}, z{}, out{};
    out.fill(42);
    const const_soa3<float> bad{x, std::span<const float>{y}.first(15), z};
    CHECK(!dot_batch(bad, bad, std::span<float>{out}));
    CHECK(out[0] == 42);
    CHECK(!transform_points(affine3f{}, bad, soa3<float>{x, y, z}));
    const const_soa3<float> input{x, y, z};
    CHECK(!dot_batch(input, input, std::span<float>{out}.first(15)));
    CHECK(!transform_points(affine3f{}, input, soa3<float>{out, out, out}));
    CHECK(out[0] == 42);
    const const_soa3<float> source{std::span<const float>{x}.first(15),
                                   std::span<const float>{y}.first(15),
                                   std::span<const float>{z}.first(15)};
    CHECK(!dot_batch(source, source, std::span<float>{x}.subspan(1)));
    CHECK(
        !transform_points(affine3f{}, source,
                          soa3<float>{std::span<float>{x}.subspan(1), std::span<float>{y}.first(15),
                                      std::span<float>{z}.first(15)}));
    std::array<vec3f, 8> packed{};
    CHECK(!transform_points<float>(affine3f{}, std::span<const vec3f>{packed}.first(7),
                                   std::span<vec3f>{packed}.subspan(1)));
    CHECK(!transform_points<float>(affine3f{}, std::span<const vec3f>{packed},
                                   std::span<vec3f>{packed}.first(7)));
    CHECK(!transform_points(affine3f{}, input, soa3<float>{std::span<float>{out}.first(15), y, z}));
}
CH_TEST(batch_exact_channel_permutation) {
    std::array<float, 5> x{1, 2, 3, 4, 5}, y{6, 7, 8, 9, 10}, z{11, 12, 13, 14, 15};
    const auto original_x = x, original_y = y, original_z = z;
    CHECK(transform_points(affine3f{}, const_soa3<float>{x, y, z}, soa3<float>{y, z, x}));
    CHECK(y == original_x && z == original_y && x == original_z);
    const auto a = x, b = y, c = z;
    CHECK(dot_batch(const_soa3<float>{x, y, z}, const_soa3<float>{x, y, z}, std::span<float>{x}));
    for (std::size_t i = 0; i < x.size(); ++i)
        CHECK(x[i] == a[i] * a[i] + b[i] * b[i] + c[i] * c[i]);
    std::printf("batch backend: %s\n", batch_backend());
}
