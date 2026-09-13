#pragma once
#include "../performance_helpers.hpp"
#include "common.hpp"
namespace extended_perf {
using namespace chm;
template <class F> CH_TEST_NOINLINE double invoke(F &f) {
    return f();
}
template <class F> double time(const char *name, std::size_t n, F work) {
    return perf::time(name, n, [&] { return invoke(work); });
}
template <class T, int Op> void run() {
    constexpr std::size_t n = 4096;
    using V = vec<T, 3>;
    if constexpr (Op < 4) {
        std::vector<T> x(n), y(n), z(n), ox(n), oy(n), oz(n), rx(n), ry(n), rz(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = T(i % 17) * T(0.1);
            y[i] = T(-2);
            z[i] = T(i % 13) * T(0.1);
        }
        const auto q = *from_axis_angle(V{T(1), T(2), T(3)}, T(0.7));
        const auto a = rotation(q);
        const const_soa3<T> input{x, y, z};
        const soa3<T> output{ox, oy, oz};
        const double reference = time("scalar vector reference", n, [&] {
            for (std::size_t i = 0; i < n; ++i) {
                const V v{x[i], y[i], z[i]};
                V r;
                if constexpr (Op == 0)
                    r = transform_vector(a, v);
                else if constexpr (Op == 1)
                    r = rotate(q, v);
                else if constexpr (Op == 2)
                    r = normalize(v);
                else
                    r = cross(v, V{z[i], x[i], y[i]});
                rx[i] = r[0];
                ry[i] = r[1];
                rz[i] = r[2];
            }
            return double(rx[n / 2]);
        });
        const double actual = time("SIMD vector batch", n, [&] {
            bool ok;
            if constexpr (Op == 0)
                ok = transform_vectors(a, input, output);
            else if constexpr (Op == 1)
                ok = rotate_vectors(q, input, output);
            else if constexpr (Op == 2)
                ok = normalize_vectors(input, output);
            else
                ok = cross_batch(input, const_soa3<T>{z, x, y}, output);
            return ok ? double(ox[n / 2]) : std::numeric_limits<double>::quiet_NaN();
        });
        for (std::size_t i = 0; i < n; ++i)
            NEAR((V{ox[i], oy[i], oz[i]}), (V{rx[i], ry[i], rz[i]}), extended::tolerance<T>);
        perf::relative(actual, reference, 6);
    } else if constexpr (Op == 4) {
        using B = mat<T, 4, 8>;
        const auto a = extended::positive_matrix<T, 4>();
        const auto lu = *factor_lu(a);
        std::vector<B> input(n / 16), out(n / 16), ref(n / 16);
        for (auto &b : input)
            for (auto &x : b.elements)
                x = test::sample<T>();
        const double reference = time("LU separate RHS", input.size(), [&] {
            for (std::size_t i = 0; i < input.size(); ++i)
                for (std::size_t j = 0; j < 8; ++j)
                    ref[i].set_column(j, lu.solve(input[i].column(j)));
            return double(ref[0](0, 0));
        });
        const double actual = time("LU shared RHS", input.size(), [&] {
            for (std::size_t i = 0; i < input.size(); ++i)
                out[i] = lu.solve(input[i]);
            return double(out[0](0, 0));
        });
        for (std::size_t i = 0; i < input.size(); ++i) {
            NEAR(out[i], ref[i], extended::tolerance<T>);
            NEAR(a * out[i], input[i], extended::tolerance<T>);
        }
        perf::relative(actual, reference, 4);
    } else if constexpr (Op == 5 || Op == 6) {
        const bezier<T, 3, 3> b{{V{T(1), T(-2), T(3)}, V{T(2), T(4), T(-1)}, V{T(-3), T(2), T(1)},
                                 V{T(4), T(3), T(2)}}};
        const auto p = *cubic_polynomial<T, 3>::from_bezier(b);
        std::vector<curve_sample<T, 3>> out(n), ref(n);
        std::vector<T> parameters(n);
        for (std::size_t i = 0; i < n; ++i)
            parameters[i] = T(i) / T(n - 1);
        const double reference = time("Bezier separate value slope", n, [&] {
            for (std::size_t i = 0; i < n; ++i)
                ref[i] = {b.evaluate(parameters[i]), b.derivative().evaluate(parameters[i])};
            return double(ref[n / 2].position[0]);
        });
        const double actual =
            time(Op == 5 ? "Bezier joint value slope" : "cubic Horner value slope", n, [&] {
                for (std::size_t i = 0; i < n; ++i) {
                    if constexpr (Op == 5)
                        out[i] = b.evaluate_with_derivative(parameters[i]);
                    else
                        out[i] = p.evaluate_with_derivative(parameters[i]);
                }
                return double(out[n / 2].position[0]);
            });
        for (std::size_t i = 0; i < n; ++i) {
            NEAR(out[i].position, ref[i].position, extended::tolerance<T>);
            NEAR(out[i].derivative, ref[i].derivative, extended::tolerance<T>);
        }
        perf::relative(actual, reference, 4);
    } else if constexpr (Op == 7) {
        const ray<T, 3> r{V{T(-3), T(-2), T(3)}, V{T(0.125), T(0.25), T(-1)}};
        const auto p = *prepared_ray<T>::create(r);
        std::vector<aabb<T, 3>> boxes(n);
        std::vector<std::optional<ray_interval<T>>> out(n), ref(n);
        for (std::size_t i = 0; i < n; ++i) {
            const V lo{T(int(i % 17) - 8), T(int(i % 11) - 5), T(int(i % 7) - 3)};
            boxes[i] = {lo, lo + V{T(2), T(2), T(2)}};
        }
        const double reference = time("ray per box", n, [&] {
            double hits{};
            for (std::size_t i = 0; i < n; ++i) {
                ref[i] = intersect(r, boxes[i]);
                hits += ref[i].has_value();
            }
            return hits;
        });
        const double actual = time("cached ray batch", n, [&] {
            if (!intersect_many(p, std::span<const aabb<T, 3>>{boxes},
                                std::span<std::optional<ray_interval<T>>>{out}))
                return -1.0;
            double hits{};
            for (const auto &v : out)
                hits += v.has_value();
            return hits;
        });
        for (std::size_t i = 0; i < n; ++i) {
            CHECK(bool(out[i]) == bool(ref[i]));
            if (out[i]) {
                NEAR(out[i]->enter, ref[i]->enter, extended::tolerance<T>);
                NEAR(out[i]->exit, ref[i]->exit, extended::tolerance<T>);
            }
        }
        perf::relative(actual, reference, 4);
    } else if constexpr (Op == 8) {
        const auto f = *extract_frustum(*perspective(T(1), T(1.5), T(0.1), T(100)));
        std::vector<aabb<T, 3>> boxes(n);
        std::vector<containment> out(n), ref(n);
        for (std::size_t i = 0; i < n; ++i) {
            const V lo{T(int(i % 17) - 8), T(int(i % 11) - 5), T(-int(i % 103))};
            boxes[i] = {lo, lo + V{T(1), T(1), T(1)}};
        }
        const double reference = time("frustum per box", n, [&] {
            for (std::size_t i = 0; i < n; ++i)
                ref[i] = classify(f, boxes[i]);
            return double(int(ref[n / 2]));
        });
        const double actual = time("frustum batch", n, [&] {
            if (!classify_batch(f, std::span<const aabb<T, 3>>{boxes}, std::span<containment>{out}))
                return -1.0;
            return double(int(out[n / 2]));
        });
        CHECK(out == ref);
        perf::relative(actual, reference, 4);
    } else {
        const mat<T, 3, 4> a = mat<T, 3, 4>::from_rows({vec<T, 4>{T(1), T(2), T(3), T(4)},
                                                        vec<T, 4>{T(2), T(1), T(-1), T(3)},
                                                        vec<T, 4>{T(3), T(2), T(1), T(-1)}});
        std::vector<mat<T, 4, 3>> input(n / 4);
        std::vector<mat<T, 3, 3>> out(n / 4), ref(n / 4);
        for (auto &b : input)
            for (auto &x : b.elements)
                x = test::sample<T>();
        const double reference = time("rectangular triple loop", input.size(), [&] {
            for (std::size_t t = 0; t < input.size(); ++t)
                for (std::size_t i = 0; i < 3; ++i)
                    for (std::size_t j = 0; j < 3; ++j) {
                        T sum{};
                        for (std::size_t k = 0; k < 4; ++k)
                            sum += a(i, k) * input[t](k, j);
                        ref[t](i, j) = sum;
                    }
            return double(ref[0](0, 0));
        });
        const double actual = time("rectangular product", input.size(), [&] {
            for (std::size_t i = 0; i < input.size(); ++i)
                out[i] = a * input[i];
            return double(out[0](0, 0));
        });
        for (std::size_t i = 0; i < input.size(); ++i)
            NEAR(out[i], ref[i], extended::tolerance<T>);
        perf::relative(actual, reference, 4);
    }
}
} // namespace extended_perf
