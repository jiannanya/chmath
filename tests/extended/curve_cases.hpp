#pragma once
#include "common.hpp"
namespace extended_curve {
using namespace chm;
template <class T, std::size_t D> bezier<T, 3, D> curve() {
    bezier<T, 3, D> c;
    for (auto &p : c.control)
        p = extended::random_vector<T, 3>();
    return c;
}
template <class T, std::size_t D> vec<T, 3> bernstein(const bezier<T, 3, D> &c, T t) {
    std::array<long double, 3> sum{};
    long double choose = 1;
    for (std::size_t i = 0; i <= D; ++i) {
        const long double weight =
            choose * std::pow(static_cast<long double>(t), static_cast<int>(i)) *
            std::pow(1 - static_cast<long double>(t), static_cast<int>(D - i));
        for (std::size_t k = 0; k < 3; ++k)
            sum[k] += weight * c.control[i][k];
        if (i < D)
            choose = choose * static_cast<long double>(D - i) / static_cast<long double>(i + 1);
    }
    return {T(sum[0]), T(sum[1]), T(sum[2])};
}
template <class T, std::size_t D> void polynomial_oracle() {
    const auto c = curve<T, D>();
    for (T t : {T(0), T(0.01), T(0.25), T(0.75), T(0.99), T(1)})
        NEAR(c.evaluate(t), bernstein(c, t), extended::tolerance<T>);
}
template <class T, std::size_t D> void joint_derivative() {
    const auto c = curve<T, D>();
    const auto d = c.derivative();
    for (T t : {T(0), T(0.1), T(0.5), T(0.9), T(1)}) {
        const auto sample = c.evaluate_with_derivative(t);
        NEAR(sample.position, bernstein(c, t), extended::tolerance<T>);
        NEAR(sample.derivative, bernstein(d, t), T(4) * extended::tolerance<T>);
    }
}
template <class T, std::size_t D> void split_reparameterization() {
    const auto c = curve<T, D>();
    for (T cut : {T(0), T(0.25), T(1)}) {
        const auto [a, b] = c.split(cut);
        for (T t : {T(0), T(0.33), T(0.75), T(1)}) {
            NEAR(a.evaluate(t), c.evaluate(cut * t), extended::tolerance<T>);
            NEAR(b.evaluate(t), c.evaluate(cut + (T(1) - cut) * t), extended::tolerance<T>);
        }
    }
}
template <class T, std::size_t D> void spline_bezier_equivalence() {
    const auto c = curve<T, D>();
    std::array<T, 2 * (D + 1)> knots{};
    std::fill(knots.begin() + D + 1, knots.end(), T(1));
    const auto s = bspline_view<T, 3, D>::create(c.control, knots);
    CHECK(s);
    for (T t : {T(0), T(0.01), T(0.25), T(0.75), T(0.99), T(1)}) {
        const auto p = s->evaluate(t), d = s->derivative(t);
        CHECK(p && d);
        NEAR(*p, bernstein(c, t), extended::tolerance<T>);
        NEAR(*d, bernstein(c.derivative(), t), T(4) * extended::tolerance<T>);
    }
}
template <class T, std::size_t D> void rational_reference() {
    const auto c = curve<T, D>();
    std::array<T, D + 1> weights{}, small_weights{};
    std::array<T, 2 * (D + 1)> knots{};
    std::fill(knots.begin() + D + 1, knots.end(), T(1));
    const T scale = std::same_as<T, float> ? T(1e-30L) : T(1e-250L);
    for (std::size_t i = 0; i <= D; ++i) {
        weights[i] = T(1 + i % 3);
        small_weights[i] = weights[i] * scale;
    }
    const auto s = nurbs_view<T, 3, D>::create(c.control, weights, knots),
               small = nurbs_view<T, 3, D>::create(c.control, small_weights, knots);
    CHECK(s && small);
    for (T t : {T(0), T(0.1), T(0.3), T(0.8), T(1)}) {
        long double denominator = 0, choose = 1;
        std::array<long double, 3> numerator{};
        for (std::size_t i = 0; i <= D; ++i) {
            const long double w = weights[i] * choose *
                                  std::pow(static_cast<long double>(t), int(i)) *
                                  std::pow(1 - static_cast<long double>(t), int(D - i));
            denominator += w;
            for (std::size_t k = 0; k < 3; ++k)
                numerator[k] += w * c.control[i][k];
            if (i < D)
                choose = choose * static_cast<long double>(D - i) / static_cast<long double>(i + 1);
        }
        const vec<T, 3> expected{T(numerator[0] / denominator), T(numerator[1] / denominator),
                                 T(numerator[2] / denominator)};
        const auto p = s->evaluate(t), q = small->evaluate(t);
        CHECK(p && q);
        NEAR(*p, expected, extended::tolerance<T>);
        NEAR(*q, expected, extended::tolerance<T>);
    }
}
template <class T, std::size_t D> void nonclamped_derivative() {
    std::array<vec<T, 3>, D + 3> controls{};
    std::array<T, 2 * D + 4> knots{};
    for (auto &p : controls)
        p = extended::random_vector<T, 3>();
    for (std::size_t i = 0; i < knots.size(); ++i)
        knots[i] = T(i);
    const auto s = bspline_view<T, 3, D>::create(controls, knots);
    CHECK(s);
    const T h = std::same_as<T, float> ? T(0.01) : T(1e-5);
    for (T t : {T(D) + T(0.35), T(D) + T(1.35), T(D) + T(2.35)}) {
        const auto d = s->derivative(t), left = s->evaluate(t - h), right = s->evaluate(t + h);
        CHECK(d && left && right);
        NEAR(*d, (*right - *left) / (T(2) * h), std::same_as<T, float> ? T(2e-3) : T(2e-7));
    }
}
} // namespace extended_curve
