#pragma once

#include <chmath/vector/vector.hpp>
#include <span>
#include <utility>

namespace chm {

namespace detail {
template <floating T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> lerp_finite(const vec<T, N> &a, const vec<T, N> &b,
                                              T t) noexcept {
    vec<T, N> r;
    for (std::size_t i = 0; i < N; ++i)
        r[i] = lerp_finite(a[i], b[i], t);
    return r;
}
} // namespace detail

template <floating T, std::size_t N> struct curve_sample {
    vec<T, N> position{}, derivative{};
};

template <floating T, std::size_t N, std::size_t Degree> struct bezier {
    std::array<vec<T, N>, Degree + 1> control{};
    [[nodiscard]] constexpr curve_sample<T, N> evaluate_with_derivative(T t) const noexcept {
        if constexpr (Degree == 0) {
            return {control[0], {}};
        } else {
            auto work = control;
            for (std::size_t k = Degree; k > 1; --k)
                for (std::size_t i = 0; i < k; ++i)
                    work[i] = detail::lerp_finite(work[i], work[i + 1], t);
            return {detail::lerp_finite(work[0], work[1], t), (work[1] - work[0]) * T(Degree)};
        }
    }
    [[nodiscard]] constexpr vec<T, N> evaluate(T t) const noexcept {
        auto work = control;
        for (std::size_t k = Degree; k > 0; --k)
            for (std::size_t i = 0; i < k; ++i)
                work[i] = detail::lerp_finite(work[i], work[i + 1], t);
        return work[0];
    }
    [[nodiscard]] constexpr auto derivative() const noexcept
        requires(Degree > 0)
    {
        bezier<T, N, Degree - 1> result;
        for (std::size_t i = 0; i < Degree; ++i)
            result.control[i] = (control[i + 1] - control[i]) * T(Degree);
        return result;
    }
    [[nodiscard]] constexpr std::pair<bezier, bezier> split(T t) const noexcept {
        auto work = control;
        bezier left, right;
        left.control[0] = work[0];
        right.control[Degree] = work[Degree];
        for (std::size_t k = Degree; k > 0; --k) {
            for (std::size_t i = 0; i < k; ++i)
                work[i] = detail::lerp_finite(work[i], work[i + 1], t);
            left.control[Degree - k + 1] = work[0];
            right.control[k - 1] = work[k - 1];
        }
        return {left, right};
    }
};
template <floating T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> hermite(const vec<T, N> &a, const vec<T, N> &tangent_a,
                                          const vec<T, N> &b, const vec<T, N> &tangent_b,
                                          T t) noexcept {
    const T t2 = t * t, t3 = t2 * t;
    return a * (T(2) * t3 - T(3) * t2 + T(1)) + tangent_a * (t3 - T(2) * t2 + t) +
           b * (-T(2) * t3 + T(3) * t2) + tangent_b * (t3 - t2);
}
// Uniform Catmull-Rom, tension=0 gives the standard interpolating cubic p1 -> p2.
template <floating T, std::size_t N>
[[nodiscard]] constexpr vec<T, N> catmull_rom(const vec<T, N> &p0, const vec<T, N> &p1,
                                              const vec<T, N> &p2, const vec<T, N> &p3, T t,
                                              T tension = T(0)) noexcept {
    const T s = (T(1) - tension) / T(2);
    return hermite(p1, (p2 - p0) * s, p2, (p3 - p1) * s, t);
}
template <floating T, std::size_t N, std::size_t UDegree, std::size_t VDegree> struct bezier_patch {
    std::array<std::array<vec<T, N>, UDegree + 1>, VDegree + 1> control{}; // [v][u]
    [[nodiscard]] constexpr vec<T, N> evaluate(T u, T v) const noexcept {
        bezier<T, N, VDegree> along_v;
        for (std::size_t j = 0; j <= VDegree; ++j)
            along_v.control[j] = bezier<T, N, UDegree>{control[j]}.evaluate(u);
        return along_v.evaluate(v);
    }
    [[nodiscard]] constexpr vec<T, N> derivative_u(T u, T v) const noexcept
        requires(UDegree > 0)
    {
        bezier<T, N, VDegree> along_v;
        for (std::size_t j = 0; j <= VDegree; ++j)
            along_v.control[j] = bezier<T, N, UDegree>{control[j]}.derivative().evaluate(u);
        return along_v.evaluate(v);
    }
    [[nodiscard]] constexpr vec<T, N> derivative_v(T u, T v) const noexcept
        requires(VDegree > 0)
    {
        bezier<T, N, VDegree> along_v;
        for (std::size_t j = 0; j <= VDegree; ++j)
            along_v.control[j] = bezier<T, N, UDegree>{control[j]}.evaluate(u);
        return along_v.derivative().evaluate(v);
    }
};

namespace detail {
template <floating T, std::size_t Degree>
[[nodiscard]] bool valid_knots(std::span<const T> knots, std::size_t count) noexcept {
    if (count < Degree + 1 || count > std::numeric_limits<std::size_t>::max() - Degree - 1 ||
        knots.size() != count + Degree + 1)
        return false;
    std::size_t multiplicity = 0;
    T previous{};
    for (std::size_t i = 0; i < knots.size(); ++i) {
        const T k = knots[i];
        if (!is_finite(k) || (i > 0 && k < previous))
            return false;
        multiplicity = i > 0 && k == previous ? multiplicity + 1 : 1;
        if (multiplicity > Degree + 1)
            return false;
        previous = k;
    }
    return knots[Degree] < knots[count];
}
template <floating T, std::size_t Degree>
[[nodiscard]] std::size_t knot_span(std::span<const T> knots, std::size_t count, T u) noexcept {
    if (u == knots[count])
        return static_cast<std::size_t>(
                   std::lower_bound(knots.begin() + Degree + 1,
                                    knots.begin() + static_cast<std::ptrdiff_t>(count) + 1, u) -
                   knots.begin()) -
               1;
    return static_cast<std::size_t>(
               std::upper_bound(knots.begin() + Degree,
                                knots.begin() + static_cast<std::ptrdiff_t>(count) + 1, u) -
               knots.begin()) -
           1;
}
template <floating T, std::size_t N, std::size_t Degree>
[[nodiscard]] vec<T, N> de_boor(std::array<vec<T, N>, Degree + 1> &work, std::span<const T> knots,
                                std::size_t span, T u) noexcept {
    for (std::size_t r = 1; r <= Degree; ++r)
        for (std::size_t j = Degree; j >= r; --j) {
            const std::size_t i = span - Degree + j;
            const T hi = knots[i + Degree - r + 1], lo = knots[i];
            const T denominator = hi - lo;
            const T alpha =
                denominator == T(0)
                    ? T(0)
                    : (is_finite(denominator) ? (u - lo) / denominator
                                              : (u / T(2) - lo / T(2)) / (hi / T(2) - lo / T(2)));
            work[j] = lerp_finite(work[j - 1], work[j], alpha);
        }
    return work[Degree];
}
} // namespace detail

// Non-owning validated view. Control and knot arrays must remain alive and unchanged.
// Degree is compile-time bounded; evaluation uses O(Degree*N) stack space, O(Degree^2*N) work.
template <floating T, std::size_t N, std::size_t Degree>
    requires(Degree <= 32)
class bspline_view {
    std::span<const vec<T, N>> control_;
    std::span<const T> knots_;
    bspline_view(std::span<const vec<T, N>> control, std::span<const T> knots) noexcept
        : control_(control), knots_(knots) {}

  public:
    [[nodiscard]] static std::optional<bspline_view> create(std::span<const vec<T, N>> control,
                                                            std::span<const T> knots) noexcept {
        if (!detail::valid_knots<T, Degree>(knots, control.size()))
            return std::nullopt;
        for (const auto &p : control)
            if (!is_finite(p))
                return std::nullopt;
        return bspline_view(control, knots);
    }
    [[nodiscard]] T domain_min() const noexcept { return knots_[Degree]; }
    [[nodiscard]] T domain_max() const noexcept { return knots_[control_.size()]; }
    [[nodiscard]] std::optional<vec<T, N>> evaluate(T u) const noexcept {
        if (!is_finite(u) || u < domain_min() || u > domain_max())
            return std::nullopt;
        const auto span = detail::knot_span<T, Degree>(knots_, control_.size(), u);
        std::array<vec<T, N>, Degree + 1> work;
        for (std::size_t j = 0; j <= Degree; ++j)
            work[j] = control_[span - Degree + j];
        const auto result = detail::de_boor<T, N, Degree>(work, knots_, span, u);
        if (!is_finite(result))
            return std::nullopt;
        return result;
    }
    [[nodiscard]] std::optional<vec<T, N>> derivative(T u) const noexcept
        requires(Degree > 0)
    {
        if (!is_finite(u) || u < domain_min() || u > domain_max())
            return std::nullopt;
        const auto span = detail::knot_span<T, Degree>(knots_, control_.size(), u);
        std::array<vec<T, N>, Degree> work;
        for (std::size_t j = 0; j < Degree; ++j) {
            const std::size_t i = span - Degree + j;
            const T d = knots_[i + Degree + 1] - knots_[i + 1];
            if (d == T(0)) {
                work[j] = {};
                continue;
            }
            const T factor =
                is_finite(d)
                    ? T(Degree) / d
                    : (T(Degree) / (knots_[i + Degree + 1] / T(2) - knots_[i + 1] / T(2))) / T(2);
            work[j] = (control_[i + 1] - control_[i]) * factor;
            if (!is_finite(work[j])) {
                for (std::size_t k = 0; k < N; ++k) {
                    const T delta = control_[i + 1][k] - control_[i][k];
                    if (!is_finite(d)) {
                        const T half_width = knots_[i + Degree + 1] / T(2) - knots_[i + 1] / T(2);
                        const T half_delta = control_[i + 1][k] / T(2) - control_[i][k] / T(2);
                        work[j][k] = multiply_divide(half_delta, T(Degree), half_width);
                    } else if (!is_finite(delta)) {
                        const T half_delta = control_[i + 1][k] / T(2) - control_[i][k] / T(2);
                        work[j][k] = multiply_divide(half_delta, T(2 * Degree), d);
                    } else {
                        work[j][k] = multiply_divide(delta, T(Degree), d);
                    }
                }
            }
        }
        const auto result = detail::de_boor<T, N, Degree - 1>(
            work, knots_.subspan(1, knots_.size() - 2), span - 1, u);
        if (!is_finite(result))
            return std::nullopt;
        return result;
    }
};
template <floating T, std::size_t N, std::size_t Degree>
    requires(Degree <= 32)
class nurbs_view {
    std::span<const vec<T, N>> control_;
    std::span<const T> weights_, knots_;
    nurbs_view(std::span<const vec<T, N>> control, std::span<const T> weights,
               std::span<const T> knots) noexcept
        : control_(control), weights_(weights), knots_(knots) {}

  public:
    [[nodiscard]] static std::optional<nurbs_view> create(std::span<const vec<T, N>> control,
                                                          std::span<const T> weights,
                                                          std::span<const T> knots) noexcept {
        if (weights.size() != control.size() ||
            !detail::valid_knots<T, Degree>(knots, control.size()))
            return std::nullopt;
        for (std::size_t i = 0; i < control.size(); ++i)
            if (!is_finite(control[i]) || !is_finite(weights[i]) || weights[i] <= T(0))
                return std::nullopt;
        return nurbs_view(control, weights, knots);
    }
    [[nodiscard]] T domain_min() const noexcept { return knots_[Degree]; }
    [[nodiscard]] T domain_max() const noexcept { return knots_[control_.size()]; }
    [[nodiscard]] std::optional<vec<T, N>> evaluate(T u) const noexcept {
        if (!is_finite(u) || u < domain_min() || u > domain_max())
            return std::nullopt;
        const auto span = detail::knot_span<T, Degree>(knots_, control_.size(), u);
        std::array<vec<T, N + 1>, Degree + 1> work;
        T scale{};
        for (std::size_t j = 0; j <= Degree; ++j)
            scale = std::max(scale, weights_[span - Degree + j]);
        for (std::size_t j = 0; j <= Degree; ++j) {
            const std::size_t i = span - Degree + j;
            const T w = weights_[i] / scale;
            for (std::size_t k = 0; k < N; ++k)
                work[j][k] = control_[i][k] * w;
            work[j][N] = w;
        }
        const auto h = detail::de_boor<T, N + 1, Degree>(work, knots_, span, u);
        if (h[N] <= T(0))
            return std::nullopt;
        vec<T, N> result;
        for (std::size_t k = 0; k < N; ++k)
            result[k] = h[k] / h[N];
        if (!is_finite(result))
            return std::nullopt;
        return result;
    }
};

} // namespace chm
