#pragma once

#include <array>
#include <chmath/core/scalar.hpp>
#include <cstdint>
#include <numeric>
#include <optional>
#include <span>

namespace chm {

template <floating T> struct quadratic_roots {
    std::array<T, 2> values{}; // Ascending, distinct real roots.
    std::size_t count{};
    bool all_reals = false;
};
template <floating T>
[[nodiscard]] inline std::optional<quadratic_roots<T>> solve_quadratic(T a, T b, T c) noexcept {
    if (!is_finite(a) || !is_finite(b) || !is_finite(c))
        return std::nullopt;
    quadratic_roots<T> r;
    if (a == T(0)) {
        if (b == T(0)) {
            r.all_reals = c == T(0);
            return r;
        }
        r.values[0] = -c / b;
        r.count = 1;
        if (!is_finite(r.values[0]))
            return std::nullopt;
        return r;
    }
    const T scale = std::max({std::abs(a), std::abs(b), std::abs(c)});
    a /= scale;
    b /= scale;
    c /= scale;
    // Scaling can erase a coefficient outside T's dynamic range; do not divide by zero.
    if (a == T(0))
        return std::nullopt;
    const T bb = b * b;
    const T d = std::fma(-T(4) * a, c, bb) + std::fma(b, b, -bb);
    if (d < T(0))
        return r;
    if (d == T(0)) {
        r.values[0] = -b / (T(2) * a);
        r.count = 1;
    } else {
        const T q = -T(0.5) * (b + std::copysign(std::sqrt(d), b));
        r.values = {q / a, c / q};
        r.count = 2;
        if (r.values[0] > r.values[1])
            std::swap(r.values[0], r.values[1]);
    }
    for (std::size_t i = 0; i < r.count; ++i)
        if (!is_finite(r.values[i]))
            return std::nullopt;
    return r;
}
// Coefficients are ordered from constant term to highest degree; empty polynomial = 0.
template <floating T>
[[nodiscard]] constexpr T evaluate_polynomial(std::span<const T> coefficients, T x) noexcept {
    T r{};
    for (std::size_t i = coefficients.size(); i-- > 0;)
        r = r * x + coefficients[i];
    return r;
}
template <floating T> struct root_result {
    T value{}, residual{};
    std::size_t iterations{};
    bool converged = false;
};
template <floating T, class F>
[[nodiscard]] std::optional<root_result<T>> bisect(F &&f, T lo, T hi, T x_tolerance = epsilon<T>,
                                                   T f_tolerance = epsilon<T>,
                                                   std::size_t max_iterations = 128) {
    if (!is_finite(lo) || !is_finite(hi) || lo > hi || !is_finite(x_tolerance) ||
        !is_finite(f_tolerance) || x_tolerance < T(0) || f_tolerance < T(0))
        return std::nullopt;
    T fl = f(lo), fh = f(hi);
    if (!is_finite(fl) || !is_finite(fh))
        return std::nullopt;
    if (std::abs(fl) <= f_tolerance)
        return root_result<T>{lo, fl, 0, true};
    if (std::abs(fh) <= f_tolerance)
        return root_result<T>{hi, fh, 0, true};
    if (std::signbit(fl) == std::signbit(fh))
        return std::nullopt;
    root_result<T> r{lo, fl, 0, false};
    for (std::size_t i = 0; i < max_iterations; ++i) {
        const T mid = std::midpoint(lo, hi), fm = f(mid);
        if (!is_finite(fm))
            return std::nullopt;
        r = {mid, fm, i + 1, false};
        if (std::abs(fm) <= f_tolerance || std::abs(hi / T(2) - lo / T(2)) <= x_tolerance ||
            mid == lo || mid == hi) {
            r.converged = true;
            return r;
        }
        if (std::signbit(fl) == std::signbit(fm)) {
            lo = mid;
            fl = fm;
        } else {
            hi = mid;
        }
    }
    return r;
}

template <floating T> struct integration_result {
    T value{}, estimated_error{};
    std::size_t evaluations{};
    bool converged = true;
};
namespace detail {
// Keep the usual short expression on ordinary ranges; scale exponents only when
// an intermediate overflows or erases a nonzero integral through underflow.
template <floating T> T simpson_panel(T a, T b, T fa, T fm, T fb) noexcept {
    const T width = b - a, sum = fa + T(4) * fm + fb;
    const T ordinary = (width / T(6)) * sum;
    if (is_finite(ordinary) && (ordinary != T(0) || sum == T(0)))
        return ordinary;
    const T scale = std::max({std::abs(fa), std::abs(fm), std::abs(fb)});
    if (scale == T(0))
        return T(0);
    const T average = ((fa / scale) + T(4) * (fm / scale) + (fb / scale)) / T(6);
    const bool wide = !is_finite(width);
    int ew{}, es{}, ef{};
    const T mw = std::frexp(wide ? b / T(2) - a / T(2) : width, &ew);
    const T ms = std::frexp(scale, &es), mf = std::frexp(average, &ef);
    return std::ldexp(mw * ms * mf, ew + es + ef + (wide ? 1 : 0));
}
template <floating T, class F>
bool simpson_recurse(F &f, T a, T b, T fa, T fm, T fb, T whole, T tol, unsigned depth,
                     std::size_t budget, integration_result<T> &r) {
    const T m = std::midpoint(a, b), l = std::midpoint(a, m), h = std::midpoint(m, b);
    if (depth == 0 || r.evaluations > budget - 2 || l == a || h == b) {
        r.value += whole;
        r.estimated_error = std::numeric_limits<T>::infinity();
        r.converged = false;
        return true;
    }
    const T fl = f(l), fh = f(h);
    r.evaluations += 2;
    if (!is_finite(fl) || !is_finite(fh))
        return false;
    const T left = simpson_panel(a, m, fa, fl, fm), right = simpson_panel(m, b, fm, fh, fb),
            delta = left + right - whole;
    if (!is_finite(left) || !is_finite(right) || !is_finite(delta))
        return false;
    if (std::abs(delta) <= T(15) * tol) {
        r.value += left + right + delta / T(15);
        r.estimated_error += std::abs(delta) / T(15);
        return true;
    }
    return simpson_recurse(f, a, m, fa, fl, fm, left, tol / T(2), depth - 1, budget, r) &&
           simpson_recurse(f, m, b, fm, fh, fb, right, tol / T(2), depth - 1, budget, r);
}
} // namespace detail
// Bounded stack and evaluation budget; error is an estimate, not a rigorous interval bound.
template <floating T, class F>
[[nodiscard]] std::optional<integration_result<T>>
integrate(F &&f, T a, T b, T absolute_tolerance = epsilon<T>, unsigned max_depth = 20,
          std::size_t max_evaluations = 100001) {
    if (!is_finite(a) || !is_finite(b) || !is_finite(absolute_tolerance) ||
        absolute_tolerance <= T(0) || max_depth > 32 || max_evaluations < 3)
        return std::nullopt;
    integration_result<T> r;
    if (a == b)
        return r;
    const bool reversed = a > b;
    if (reversed)
        std::swap(a, b);
    const T fa = f(a), fm = f(std::midpoint(a, b)), fb = f(b);
    if (!is_finite(fa) || !is_finite(fm) || !is_finite(fb))
        return std::nullopt;
    r.evaluations = 3;
    const T whole = detail::simpson_panel(a, b, fa, fm, fb);
    if (!is_finite(whole) ||
        !detail::simpson_recurse(f, a, b, fa, fm, fb, whole, absolute_tolerance, max_depth,
                                 max_evaluations, r) ||
        !is_finite(r.value))
        return std::nullopt;
    if (reversed)
        r.value = -r.value;
    return r;
}

// Compensated summation for CAD measurements and reductions.
template <floating T> class compensated_sum {
    T sum_{}, correction_{};

  public:
    constexpr void add(T x) noexcept {
        const T next = sum_ + x;
        correction_ += std::abs(sum_) >= std::abs(x) ? (sum_ - next) + x : (x - next) + sum_;
        sum_ = next;
    }
    [[nodiscard]] constexpr T value() const noexcept { return sum_ + correction_; }
    constexpr void reset() noexcept { sum_ = correction_ = T(0); }
};

} // namespace chm
