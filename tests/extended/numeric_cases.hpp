#pragma once
#include "common.hpp"
namespace extended_numeric {
using namespace chm;
template <class T, int Case> void quadratic_roots_known() {
    const T lo = T(-2 - Case), hi = T(1 + Case),
            scale = Case % 2                 ? T(1)
                    : std::same_as<T, float> ? T(1e-20L)
                                             : T(1e-150L);
    const auto r = solve_quadratic(scale, -(lo + hi) * scale, lo * hi * scale);
    CHECK(r && r->count == 2);
    NEAR(r->values[0], lo, extended::tolerance<T>);
    NEAR(r->values[1], hi, extended::tolerance<T>);
}
template <class T, int Case> void horner_reference() {
    constexpr std::array<std::size_t, 4> degrees{1, 3, 8, 16};
    std::array<T, degrees[Case] + 1> c{};
    for (T &v : c)
        v = test::sample<T>(T(-1), T(1));
    for (T x : {T(-0.3), T(0), T(0.5), T(1.3)}) {
        long double sum = 0, power = 1;
        for (T v : c) {
            sum += static_cast<long double>(v) * power;
            power *= x;
        }
        NEAR(evaluate_polynomial<T>(c, x), T(sum), T(4) * extended::tolerance<T>);
    }
}
template <class T, int Case> void bisect_irrational() {
    constexpr std::array<int, 4> values{2, 3, 5, 7};
    const T value = T(values[Case]);
    const auto r = bisect([value](T x) { return x * x - value; }, T(0), T(4),
                          extended::tolerance<T> / T(10), extended::tolerance<T> / T(10));
    CHECK(r && r->converged);
    NEAR(r->value, std::sqrt(value), extended::tolerance<T>);
    CHECK(r->iterations <= 128);
}
template <class T, int Case> void integrate_polynomial() {
    constexpr int degree = Case + 2;
    const T hi = T(Case + 1) / T(2);
    const T tol = std::same_as<T, float> ? T(1e-5) : T(1e-10);
    const auto r = integrate(
        [](T x) {
            T y = T(1);
            for (int i = 0; i < degree; ++i)
                y *= x;
            return y;
        },
        T(0), hi, tol, 24, 8193);
    CHECK(r && r->converged);
    NEAR(r->value, std::pow(hi, T(degree + 1)) / T(degree + 1), T(8) * extended::tolerance<T>);
    CHECK(r->evaluations <= 8193);
}
template <class T, int Case> void compensated_cancellation() {
    constexpr std::array<int, 4> counts{1, 8, 256, 4096};
    const T large = std::same_as<T, float> ? T(1e8) : T(1e20);
    compensated_sum<T> sum;
    for (int i = 0; i < counts[Case]; ++i) {
        sum.add(large);
        sum.add(T(0.25));
        sum.add(-large);
    }
    CHECK(sum.value() == T(counts[Case]) * T(0.25));
    sum.reset();
    CHECK(sum.value() == T(0));
}
template <class T, int Case> void bounded_root_work() {
    constexpr std::array<std::size_t, 4> budgets{0, 1, 4, 8};
    std::size_t calls{};
    const auto r = bisect(
        [&](T x) {
            ++calls;
            return x * x * x - T(2);
        },
        T(0), T(2), T(0), T(0), budgets[Case]);
    CHECK(r && !r->converged);
    CHECK(r->iterations == budgets[Case]);
    CHECK(calls == budgets[Case] + 2);
    CHECK(r->value >= T(0) && r->value <= T(2));
}
} // namespace extended_numeric
