#pragma once
#include "common.hpp"
namespace extended_stress {
using namespace chm;
template <class T, int Work, int Tier> void run() {
    constexpr std::size_t n = 512U << (3 * Tier);
    using V = vec<T, 3>;
    using M = mat<T, 3, 3>;
    std::size_t completed{};
    double checksum{};
    if constexpr (Work == 0) {
        const V direction = normalize(V{T(1), T(-2), T(3)});
        for (std::size_t i = 0; i < n; ++i) {
            const T scale = std::ldexp(T(1), i % 2 ? (std::same_as<T, float> ? 100 : 800)
                                                   : (std::same_as<T, float> ? -100 : -800));
            const auto v = normalize(direction * scale);
            NEAR(v, direction, extended::tolerance<T>);
            checksum += v[0];
            ++completed;
        }
    } else if constexpr (Work == 1) {
        const T step = T(0.001);
        const auto dq = *from_axis_angle(V{T(0), T(1), T(0)}, step);
        quat<T> q;
        for (std::size_t i = 0; i < n; ++i) {
            q = dq * q;
            if ((i % 32) == 31)
                q = normalize(q);
            CHECK(is_finite(q));
            ++completed;
        }
        const auto expected = *from_axis_angle(V{T(0), T(1), T(0)}, T(n) * step);
        NEAR(rotate(q, V{T(1), T(0), T(0)}), rotate(expected, V{T(1), T(0), T(0)}),
             T(3) * extended::tolerance<T>);
        checksum = q.w();
    } else if constexpr (Work == 2) {
        const V delta{T(0.001), T(-0.002), T(0.0005)};
        const auto step = translation(delta);
        affine3<T> a;
        for (std::size_t i = 0; i < n; ++i) {
            a = step * a;
            CHECK(is_finite(a.matrix));
            ++completed;
        }
        NEAR(transform_point(a, V{}), delta * T(n), std::same_as<T, float> ? T(8e-4) : T(3e-11));
        checksum = a.matrix(0, 3);
    } else if constexpr (Work == 3) {
        const auto a = extended::positive_matrix<T, 3>();
        const auto lu = factor_lu(a);
        CHECK(lu);
        for (std::size_t i = 0; i < n; ++i) {
            const V x{T(i % 13), T(-2), T(0.25)};
            const auto solved = lu->solve(a * x);
            NEAR(solved, x, extended::tolerance<T>);
            checksum += solved[0];
            ++completed;
        }
    } else if constexpr (Work == 4) {
        const ray<T, 3> r{V{T(-3), T(-2), T(3)}, V{T(0.125), T(0.25), T(-1)}};
        const auto p = *prepared_ray<T>::create(r);
        for (std::size_t i = 0; i < n; ++i) {
            const V lo{T(int(i % 17) - 8) * T(0.5), T(int(i % 11) - 5), T(int(i % 7) - 3)};
            const aabb<T, 3> b{lo, lo + V{T(1), T(1), T(1)}};
            const auto a = intersect(r, b), c = p.intersect(b);
            CHECK(bool(a) == bool(c));
            if (a) {
                NEAR(a->enter, c->enter, extended::tolerance<T>);
                NEAR(a->exit, c->exit, extended::tolerance<T>);
                checksum += c->enter;
            }
            ++completed;
        }
    } else if constexpr (Work == 5) {
        std::array<V, 9> controls{};
        std::array<T, 18> knots{};
        std::fill(knots.begin() + 9, knots.end(), T(1));
        for (std::size_t j = 0; j < 9; ++j)
            controls[j] = V{T(j) / T(8), T(j) / T(4), T(1)};
        const auto s = bspline_view<T, 3, 8>::create(controls, knots);
        CHECK(s);
        for (std::size_t i = 0; i < n; ++i) {
            const T t = T(i % 1025) / T(1024);
            const auto p = s->evaluate(t);
            CHECK(p);
            NEAR(*p, (V{t, T(2) * t, T(1)}), extended::tolerance<T>);
            checksum += (*p)[0];
            ++completed;
        }
    } else if constexpr (Work == 6) {
        std::vector<T> x(n), y(n), z(n), ox(n), oy(n), oz(n);
        const auto q = *from_axis_angle(V{T(1), T(2), T(3)}, T(0.4));
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = T(i % 17);
            y[i] = T(-2);
            z[i] = T(i % 7);
        }
        const soa3<T> out{ox, oy, oz};
        CHECK(rotate_vectors(q, const_soa3<T>{x, y, z}, out));
        CHECK(normalize_vectors(out.as_const(), out));
        for (std::size_t i = 0; i < n; ++i) {
            NEAR((V{ox[i], oy[i], oz[i]}), normalize(rotate(q, V{x[i], y[i], z[i]})),
                 extended::tolerance<T>);
            CHECK(x[i] == T(i % 17) && y[i] == T(-2) && z[i] == T(i % 7));
            checksum += ox[i];
            ++completed;
        }
    } else if constexpr (Work == 7) {
        for (std::size_t i = 0; i < n; ++i) {
            const auto r = to_matrix(*from_axis_angle(V{T(1), T(2), T(3)}, T(i % 31) / T(17)));
            auto d = M::identity();
            d(1, 1) = std::same_as<T, float> ? T(1e-3) : T(1e-8);
            const auto a = r * d;
            const auto inv = inverse(a);
            CHECK(inv);
            NEAR(*inv * a, M::identity(), std::same_as<T, float> ? T(2e-3) : T(2e-6));
            checksum += (*inv)(0, 0);
            ++completed;
        }
    } else if constexpr (Work == 8) {
        for (std::size_t i = 0; i < n / 128; ++i) {
            const T w = T(3 + 4 * (i % 7));
            const auto result = integrate([w](T x) { return std::sin(w * x) + x * x; }, T(0), pi<T>,
                                          std::same_as<T, float> ? T(1e-4) : T(1e-9), 24, 8193);
            CHECK(result && result->converged);
            NEAR(result->value, T(2) / w + pi<T> * pi<T> * pi<T> / T(3),
                 std::same_as<T, float> ? T(2e-5) : T(2e-10));
            CHECK(result->evaluations <= 8193);
            checksum += result->value;
            completed += 128;
        }
    } else {
        constexpr std::size_t workers = 2U << Tier;
        const auto a = extended::positive_matrix<T, 3>();
        const auto lu = *factor_lu(a);
        const bezier<T, 3, 3> b{
            {V{T(0), T(1), T(2)}, V{T(1), T(2), T(3)}, V{T(2), T(3), T(4)}, V{T(3), T(4), T(5)}}};
        const auto poly = *cubic_polynomial<T, 3>::from_bezier(b);
        std::array<std::size_t, workers> failures{}, done{};
        std::array<double, workers> sums{};
        std::vector<std::jthread> threads;
        threads.reserve(workers);
        for (std::size_t w = 0; w < workers; ++w)
            threads.emplace_back([&, w] {
                for (std::size_t i = w; i < n; i += workers) {
                    const T t = T(i % 1025) / T(1024);
                    const V x{t, T(2), T(-1)};
                    const auto v = lu.solve(a * x), p = poly.evaluate(t);
                    if (!almost_equal(v, x, extended::tolerance<T>, extended::tolerance<T>) ||
                        !almost_equal(p, V{T(3) * t, T(3) * t + T(1), T(3) * t + T(2)},
                                      extended::tolerance<T>, extended::tolerance<T>))
                        ++failures[w];
                    sums[w] += v[0] + p[0];
                    ++done[w];
                }
            });
        for (auto &thread : threads)
            thread.join();
        for (std::size_t w = 0; w < workers; ++w) {
            CHECK(failures[w] == 0);
            completed += done[w];
            checksum += sums[w];
        }
    }
    CHECK(completed == n);
    CHECK(std::isfinite(checksum));
}
} // namespace extended_stress
