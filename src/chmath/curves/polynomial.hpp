#pragma once
#include <chmath/curves/curves.hpp>

namespace chm {
// Reusable power-basis cubic. Factory rejects unrepresentable coefficients;
// evaluation is the fast Horner form, requiring representable intermediates.
template <floating T, std::size_t N> struct cubic_polynomial {
    std::array<vec<T, N>, 4> coefficients{}; // constant through cubic
    [[nodiscard]] static constexpr std::optional<cubic_polynomial>
    from_bezier(const bezier<T, N, 3> &b) noexcept {
        for (const auto &p : b.control)
            if (!is_finite(p))
                return std::nullopt;
        cubic_polynomial c;
        c.coefficients[0] = b.control[0];
        c.coefficients[1] = (b.control[1] - b.control[0]) * T(3);
        c.coefficients[2] = (b.control[2] - b.control[1]) * T(3) - c.coefficients[1];
        c.coefficients[3] = b.control[3] - b.control[0] - c.coefficients[1] - c.coefficients[2];
        for (const auto &v : c.coefficients)
            if (!is_finite(v))
                return std::nullopt;
        return c;
    }
    [[nodiscard]] constexpr vec<T, N> evaluate(T t) const noexcept {
        return ((coefficients[3] * t + coefficients[2]) * t + coefficients[1]) * t +
               coefficients[0];
    }
    [[nodiscard]] constexpr vec<T, N> derivative(T t) const noexcept {
        return (coefficients[3] * (T(3) * t) + coefficients[2] * T(2)) * t + coefficients[1];
    }
    [[nodiscard]] constexpr vec<T, N> second_derivative(T t) const noexcept {
        return coefficients[3] * (T(6) * t) + coefficients[2] * T(2);
    }
    [[nodiscard]] constexpr curve_sample<T, N> evaluate_with_derivative(T t) const noexcept {
        const auto q = coefficients[3] * t + coefficients[2];
        const auto r = q * t + coefficients[1];
        return {r * t + coefficients[0], r + (q + coefficients[3] * t) * t};
    }
};
} // namespace chm
