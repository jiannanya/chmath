#pragma once
#include <chmath/core/span.hpp>
#include <chmath/simd/detail/packet.hpp>
#include <chmath/transform/transform.hpp>

namespace chm {
template <floating T> struct const_soa3 {
    std::span<const T> x, y, z;
    [[nodiscard]] bool valid() const noexcept {
        return x.size() == y.size() && x.size() == z.size();
    }
    [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
};
template <floating T> struct soa3 {
    std::span<T> x, y, z;
    [[nodiscard]] bool valid() const noexcept {
        return x.size() == y.size() && x.size() == z.size();
    }
    [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
    [[nodiscard]] const_soa3<T> as_const() const noexcept { return {x, y, z}; }
};
namespace detail {
template <floating T>
[[nodiscard]] inline bool valid_batch(const_soa3<T> in, soa3<T> out) noexcept {
    if (!in.valid() || !out.valid() || in.size() != out.size() || spans_overlap(out.x, out.y) ||
        spans_overlap(out.x, out.z) || spans_overlap(out.y, out.z))
        return false;
    for (auto dst : {out.x, out.y, out.z})
        for (auto src : {in.x, in.y, in.z})
            if (unsafe_overlap(src, dst))
                return false;
    return true;
}
template <bool Translate, floating T>
[[nodiscard]] inline bool transform_soa(const affine3<T> &a, const_soa3<T> input,
                                        soa3<T> output) noexcept {
    if (!valid_batch(input, output))
        return false;
    std::size_t i = 0;
    using P = packet<T>;
    if constexpr (P::lanes > 1) {
        typename P::value m[3][4];
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t col = 0; col < (Translate ? 4U : 3U); ++col)
                m[row][col] = P::splat(a.matrix(row, col));
        for (; input.size() - i >= P::lanes; i += P::lanes) {
            const auto x = P::load(input.x.data() + i), y = P::load(input.y.data() + i),
                       z = P::load(input.z.data() + i);
            const auto eval = [&](std::size_t row) {
                auto v = P::add(P::add(P::mul(m[row][0], x), P::mul(m[row][1], y)),
                                P::mul(m[row][2], z));
                if constexpr (Translate)
                    v = P::add(v, m[row][3]);
                return v;
            };
            P::store(output.x.data() + i, eval(0));
            P::store(output.y.data() + i, eval(1));
            P::store(output.z.data() + i, eval(2));
        }
    }
    for (; i < input.size(); ++i) {
        const vec<T, 3> v{input.x[i], input.y[i], input.z[i]};
        const auto r = Translate ? transform_point(a, v) : transform_vector(a, v);
        output.x[i] = r[0];
        output.y[i] = r[1];
        output.z[i] = r[2];
    }
    return true;
}
} // namespace detail
[[nodiscard]] inline const char *batch_backend() noexcept {
    return detail::packet<float>::name;
}
template <floating T> [[nodiscard]] inline const char *batch_backend() noexcept {
    return detail::packet<T>::name;
}

// Packed AoS, exact in-place supported. Invalid sizes/partial overlaps never write.
template <floating T>
[[nodiscard]] inline bool transform_points(const affine3<T> &a, std::span<const vec<T, 3>> input,
                                           std::span<vec<T, 3>> output) noexcept {
    if (input.size() != output.size() || detail::unsafe_overlap(input, output))
        return false;
    for (std::size_t i = 0; i < input.size(); ++i)
        output[i] = transform_point(a, input[i]);
    return true;
}
template <floating T>
[[nodiscard]] inline bool transform_vectors(const affine3<T> &a, std::span<const vec<T, 3>> input,
                                            std::span<vec<T, 3>> output) noexcept {
    if (input.size() != output.size() || detail::unsafe_overlap(input, output))
        return false;
    for (std::size_t i = 0; i < input.size(); ++i)
        output[i] = transform_vector(a, input[i]);
    return true;
}
template <floating T>
[[nodiscard]] inline bool transform_points(const affine3<T> &a, const_soa3<T> input,
                                           soa3<T> output) noexcept {
    return detail::transform_soa<true>(a, input, output);
}
template <floating T>
[[nodiscard]] inline bool transform_vectors(const affine3<T> &a, const_soa3<T> input,
                                            soa3<T> output) noexcept {
    return detail::transform_soa<false>(a, input, output);
}

// A unit quaternion is checked once; its matrix is shared by the entire batch.
template <floating T, class Input, class Output>
[[nodiscard]] inline bool rotate_vectors(const quat<T> &q, Input input, Output output) noexcept {
    if (!is_finite(q) || !almost_equal(dot(q, q), T(1), T(4) * epsilon<T>, T(4) * epsilon<T>))
        return false;
    return transform_vectors(rotation(q), input, output);
}
template <floating T>
[[nodiscard]] inline bool dot_batch(const_soa3<T> a, const_soa3<T> b,
                                    std::span<T> output) noexcept {
    if (!a.valid() || !b.valid() || a.size() != b.size() || a.size() != output.size())
        return false;
    for (auto src : {a.x, a.y, a.z, b.x, b.y, b.z})
        if (detail::unsafe_overlap(src, output))
            return false;
    std::size_t i = 0;
    using P = detail::packet<T>;
    // Keep double as a contiguous loop: compiler vectorization/unrolling is
    // faster than forcing two-lane packets in the measured Clang workload.
    if constexpr (P::lanes > 1 && std::same_as<T, float>)
        for (; a.size() - i >= P::lanes; i += P::lanes) {
            const auto x = P::mul(P::load(a.x.data() + i), P::load(b.x.data() + i));
            const auto y = P::mul(P::load(a.y.data() + i), P::load(b.y.data() + i));
            const auto z = P::mul(P::load(a.z.data() + i), P::load(b.z.data() + i));
            P::store(output.data() + i, P::add(P::add(x, y), z));
        }
    for (; i < a.size(); ++i)
        output[i] = a.x[i] * b.x[i] + a.y[i] * b.y[i] + a.z[i] * b.z[i];
    return true;
}
template <floating T>
[[nodiscard]] inline bool cross_batch(const_soa3<T> a, const_soa3<T> b, soa3<T> output) noexcept {
    if (!detail::valid_batch(a, output) || !detail::valid_batch(b, output))
        return false;
    std::size_t i = 0;
    using P = detail::packet<T>;
    if constexpr (P::lanes > 1)
        for (; a.size() - i >= P::lanes; i += P::lanes) {
            const auto ax = P::load(a.x.data() + i), ay = P::load(a.y.data() + i),
                       az = P::load(a.z.data() + i);
            const auto bx = P::load(b.x.data() + i), by = P::load(b.y.data() + i),
                       bz = P::load(b.z.data() + i);
            P::store(output.x.data() + i, P::sub(P::mul(ay, bz), P::mul(az, by)));
            P::store(output.y.data() + i, P::sub(P::mul(az, bx), P::mul(ax, bz)));
            P::store(output.z.data() + i, P::sub(P::mul(ax, by), P::mul(ay, bx)));
        }
    for (; i < a.size(); ++i) {
        const auto r = cross(vec<T, 3>{a.x[i], a.y[i], a.z[i]}, vec<T, 3>{b.x[i], b.y[i], b.z[i]});
        output.x[i] = r[0];
        output.y[i] = r[1];
        output.z[i] = r[2];
    }
    return true;
}
// Like normalize(vec): zero/nonfinite vectors produce zero. No approximate rsqrt.
template <floating T>
[[nodiscard]] inline bool normalize_vectors(const_soa3<T> input, soa3<T> output) noexcept {
    if (!detail::valid_batch(input, output))
        return false;
    std::size_t i = 0;
    using P = detail::packet<T>;
    if constexpr (P::lanes > 1 && P::vector_sqrt)
        for (; input.size() - i >= P::lanes; i += P::lanes) {
            const auto x = P::load(input.x.data() + i), y = P::load(input.y.data() + i),
                       z = P::load(input.z.data() + i);
            const auto squared = P::add(P::add(P::mul(x, x), P::mul(y, y)), P::mul(z, z));
            if (P::normal_range(squared)) {
                const auto reciprocal = P::div(P::splat(T(1)), P::sqrt(squared));
                P::store(output.x.data() + i, P::mul(x, reciprocal));
                P::store(output.y.data() + i, P::mul(y, reciprocal));
                P::store(output.z.data() + i, P::mul(z, reciprocal));
            } else
                for (std::size_t lane = 0; lane < P::lanes; ++lane) {
                    const auto j = i + lane;
                    const auto r = normalize(vec<T, 3>{input.x[j], input.y[j], input.z[j]});
                    output.x[j] = r[0];
                    output.y[j] = r[1];
                    output.z[j] = r[2];
                }
        }
    for (; i < input.size(); ++i) {
        const auto r = normalize(vec<T, 3>{input.x[i], input.y[i], input.z[i]});
        output.x[i] = r[0];
        output.y[i] = r[1];
        output.z[i] = r[2];
    }
    return true;
}
template <floating T>
[[nodiscard]] inline bool normalize_vectors(std::span<const vec<T, 3>> input,
                                            std::span<vec<T, 3>> output) noexcept {
    if (input.size() != output.size() || detail::unsafe_overlap(input, output))
        return false;
    for (std::size_t i = 0; i < input.size(); ++i)
        output[i] = normalize(input[i]);
    return true;
}
} // namespace chm
