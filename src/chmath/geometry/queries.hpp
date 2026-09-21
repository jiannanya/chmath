#pragma once
#include <chmath/core/span.hpp>
#include <chmath/geometry/frustum.hpp>
#include <chmath/geometry/intersection.hpp>
#include <chmath/parallel/parallel.hpp>

namespace chm {
// Owns its ray and cached reciprocals. Intended for many boxes sharing one ray;
// it does not require the original ray to remain alive or immutable.
template <floating T, std::size_t N = 3> class prepared_ray {
    ray<T, N> ray_;
    vec<T, N> reciprocal_{};
    explicit prepared_ray(const ray<T, N> &r) noexcept : ray_(r) {
        for (std::size_t i = 0; i < N; ++i)
            // Parallel axes never use the reciprocal; keep the divisor nonzero
            // even under aggressive constant propagation.
            reciprocal_[i] = T(1) / (r.direction[i] == T(0) ? T(1) : r.direction[i]);
    }

  public:
    [[nodiscard]] static std::optional<prepared_ray> create(const ray<T, N> &r) noexcept {
        return r.valid() ? std::optional<prepared_ray>{prepared_ray(r)} : std::nullopt;
    }
    [[nodiscard]] const ray<T, N> &source() const noexcept { return ray_; }
    [[nodiscard]] std::optional<ray_interval<T>>
    intersect(const aabb<T, N> &box, T lo = T(0),
              T hi = std::numeric_limits<T>::infinity()) const noexcept {
        if (!box.valid() || !detail::valid_interval(lo, hi))
            return std::nullopt;
        const T original_lo = lo, original_hi = hi;
        for (std::size_t i = 0; i < N; ++i) {
            if (ray_.direction[i] == T(0)) {
                if (ray_.origin[i] < box.lower[i] || ray_.origin[i] > box.upper[i])
                    return std::nullopt;
                continue;
            }
            const T dl = box.lower[i] - ray_.origin[i], dh = box.upper[i] - ray_.origin[i];
            T a =
                is_finite(dl) && is_finite(reciprocal_[i])
                    ? dl * reciprocal_[i]
                    : detail::difference_quotient(box.lower[i], ray_.origin[i], ray_.direction[i]);
            T b =
                is_finite(dh) && is_finite(reciprocal_[i])
                    ? dh * reciprocal_[i]
                    : detail::difference_quotient(box.upper[i], ray_.origin[i], ray_.direction[i]);
            if (a > b)
                std::swap(a, b);
            lo = std::max(lo, a);
            hi = std::min(hi, b);
            if (lo > hi) {
                // Resolve borderline separation using division, as in the
                // original query, instead of losing tangencies to a reciprocal.
                const T scale = std::max(std::abs(lo), std::abs(hi));
                return is_finite(lo) && is_finite(hi) &&
                               lo - hi <= T(8) * std::numeric_limits<T>::epsilon() * scale
                           ? chm::intersect(ray_, box, original_lo, original_hi)
                           : std::nullopt;
            }
        }
        if (!is_finite(lo))
            return std::nullopt;
        if (is_finite(hi) && hi - lo <= T(8) * std::numeric_limits<T>::epsilon() *
                                            std::max(std::abs(lo), std::abs(hi)))
            return chm::intersect(ray_, box, original_lo, original_hi);
        return ray_interval<T>{lo, hi};
    }
};
template <floating T, std::size_t N>
[[nodiscard]] inline std::optional<ray_interval<T>>
intersect(const prepared_ray<T, N> &r, const aabb<T, N> &b, T lo = T(0),
          T hi = std::numeric_limits<T>::infinity()) noexcept {
    return r.intersect(b, lo, hi);
}
template <floating T, std::size_t N>
[[nodiscard]] inline bool
intersect_many(const prepared_ray<T, N> &r, std::span<const aabb<T, N>> boxes,
               std::span<std::optional<ray_interval<T>>> output, T lo = T(0),
               T hi = std::numeric_limits<T>::infinity()) noexcept {
    if (boxes.size() != output.size() || !detail::valid_interval(lo, hi) ||
        detail::spans_overlap(boxes, output) ||
        detail::spans_overlap(std::span<const prepared_ray<T, N>>{&r, 1}, output))
        return false;
    for (std::size_t i = 0; i < boxes.size(); ++i)
        output[i] = r.intersect(boxes[i], lo, hi);
    return true;
}
template <floating T, class Shape>
    requires requires(const frustum<T> &f, const Shape &shape) {
        { classify(f, shape) } -> std::same_as<containment>;
    }
[[nodiscard]] inline bool classify_batch(const frustum<T> &f, std::span<const Shape> input,
                                         std::span<containment> output) noexcept {
    if (input.size() != output.size() || detail::spans_overlap(input, output) ||
        detail::spans_overlap(std::span<const frustum<T>>{&f, 1}, output))
        return false;
    for (std::size_t i = 0; i < input.size(); ++i)
        output[i] = classify(f, input[i]);
    return true;
}
// Opt-in parallel overloads. The whole-batch contract is checked first, then
// the query is split into independent cache-line-aligned ranges; see the note
// in chmath/simd/batch.hpp for the exact guarantees. The interval endpoints are
// explicit here so that the worker count is never confused with a t value.
template <floating T, std::size_t N>
[[nodiscard]] inline bool
intersect_many(const prepared_ray<T, N> &r, std::span<const aabb<T, N>> boxes,
               std::span<std::optional<ray_interval<T>>> output, T lo, T hi, unsigned workers) {
    if (boxes.size() != output.size() || !detail::valid_interval(lo, hi) ||
        detail::spans_overlap(boxes, output) ||
        detail::spans_overlap(std::span<const prepared_ray<T, N>>{&r, 1}, output))
        return false;
    parallel_for(output.size(), workers, sizeof(std::optional<ray_interval<T>>),
                 [&](std::size_t begin, std::size_t end) {
                     for (std::size_t i = begin; i < end; ++i)
                         output[i] = r.intersect(boxes[i], lo, hi);
                 });
    return true;
}
template <floating T, class Shape>
    requires requires(const frustum<T> &f, const Shape &shape) {
        { classify(f, shape) } -> std::same_as<containment>;
    }
[[nodiscard]] inline bool classify_batch(const frustum<T> &f, std::span<const Shape> input,
                                         std::span<containment> output, unsigned workers) {
    if (input.size() != output.size() || detail::spans_overlap(input, output) ||
        detail::spans_overlap(std::span<const frustum<T>>{&f, 1}, output))
        return false;
    parallel_for(output.size(), workers, sizeof(containment), [&](std::size_t begin,
                                                                  std::size_t end) {
        for (std::size_t i = begin; i < end; ++i)
            output[i] = classify(f, input[i]);
    });
    return true;
}
} // namespace chm
