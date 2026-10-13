#include "guarded_buffer.hpp"
#include "test.hpp"
#include <chmath/simd/batch.hpp>
#include <vector>

using namespace chm;
namespace batch_hotpaths {
template <floating T> void normalize_guarded(bool at_end) {
    using V = vec<T, 3>;
    const T tol = std::same_as<T, float> ? T(3e-6) : T(3e-15);
    const std::array<V, 10> values{V{T(1), T(2), T(-3)},
                                   V{},
                                   V{std::numeric_limits<T>::max()},
                                   V{std::numeric_limits<T>::denorm_min()},
                                   V{std::numeric_limits<T>::min()},
                                   V{std::numeric_limits<T>::infinity()},
                                   V{std::numeric_limits<T>::quiet_NaN()},
                                   V{T(-0.0), T(1), T(-0.0)},
                                   V{T(-4), T(3), T(2)},
                                   V{T(0.1), T(-0.2), T(0.3)}};
    for (std::size_t n = 1; n <= 19; ++n) {
        safety::guarded_buffer<V> input(n, at_end), output(n, at_end), inplace(n, at_end);
        for (std::size_t offset = 0; offset < values.size(); ++offset) {
            for (std::size_t i = 0; i < n; ++i)
                input[i] = inplace[i] = values[(i + offset) % values.size()];
            CHECK(normalize_vectors(input.const_span(), output.span()));
            CHECK(normalize_vectors(inplace.const_span(), inplace.span()));
            for (std::size_t i = 0; i < n; ++i) {
                const auto expected = normalize(input[i]);
                NEAR(output[i], expected, tol);
                CHECK(output[i] == inplace[i]);
                if (is_finite(input[i]) && input[i][1] == T(1)) {
                    CHECK(std::signbit(output[i][0]));
                    CHECK(std::signbit(output[i][2]));
                }
            }
        }
        input.read_only();
        CHECK(normalize_vectors(input.const_span(), output.span()));
    }
}
} // namespace batch_hotpaths
CH_TEST(safety_aos_normalize_float_head) {
    batch_hotpaths::normalize_guarded<float>(false);
}
CH_TEST(safety_aos_normalize_float_tail) {
    batch_hotpaths::normalize_guarded<float>(true);
}
CH_TEST(safety_aos_normalize_double_head) {
    batch_hotpaths::normalize_guarded<double>(false);
}
CH_TEST(safety_aos_normalize_double_tail) {
    batch_hotpaths::normalize_guarded<double>(true);
}
CH_TEST(safety_aos_normalize_overlap_no_write) {
    std::array<vec3d, 13> values;
    values.fill(vec3d{1, 2, 3});
    const auto original = values;
    const std::span<const vec3d> input{values};
    const std::span<vec3d> output{values};
    CHECK(!normalize_vectors(input.first(12), output.subspan(1)));
    CHECK(values == original);
    CHECK(!normalize_vectors(input.subspan(1), output.first(12)));
    CHECK(values == original);
    CHECK(!normalize_vectors(input, output.first(12)));
    CHECK(values == original);
    CHECK(normalize_vectors(std::span<const vec3d>{}, std::span<vec3d>{}));
}
