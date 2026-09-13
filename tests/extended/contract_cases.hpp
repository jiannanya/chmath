#pragma once
#include "common.hpp"
namespace extended_contract {
using namespace chm;
template <class T, int Case> void run() {
    using V = vec<T, 3>;
    const T nan = std::numeric_limits<T>::quiet_NaN(), inf = std::numeric_limits<T>::infinity();
    if constexpr (Case == 0) {
        std::array<V, 3> in{V{T(1)}, V{T(2)}, V{T(3)}}, out = in;
        CHECK(!rotate_vectors(quat<T>{T(0), T(0), T(0), T(2)}, std::span<const V>{in},
                              std::span<V>{out}));
        CHECK(out == in);
        const auto unit = *from_axis_angle(V{T(1), T(2), T(3)}, T(0.5));
        CHECK(rotate_vectors(unit, std::span<const V>{in}, std::span<V>{out}));
        for (std::size_t i = 0; i < in.size(); ++i)
            NEAR(out[i], rotate(unit, in[i]), extended::tolerance<T>);
        CHECK(!rotate_vectors(quat<T>{nan, T(0), T(0), T(1)}, std::span<const V>{in},
                              std::span<V>{out}));
        for (std::size_t i = 0; i < in.size(); ++i)
            NEAR(out[i], rotate(unit, in[i]), extended::tolerance<T>);
    } else if constexpr (Case == 1) {
        std::array<T, 9> x{T(0), nan, inf, T(3), T(1e-30), T(1e30), T(1), -inf, T(-4)}, y{}, z{},
            ox{}, oy{}, oz{};
        y[3] = T(4);
        y[8] = T(3);
        CHECK(normalize_vectors(const_soa3<T>{x, y, z}, soa3<T>{ox, oy, oz}));
        for (std::size_t i = 0; i < x.size(); ++i)
            NEAR((V{ox[i], oy[i], oz[i]}), extended::normalization_reference(V{x[i], y[i], z[i]}),
                 extended::tolerance<T>);
    } else if constexpr (Case == 2) {
        std::array<T, 12> a{};
        a.fill(T(7));
        const auto saved = a;
        CHECK(!cross_batch(const_soa3<T>{{a.data(), 4}, {a.data() + 4, 4}, {a.data() + 8, 4}},
                           const_soa3<T>{{a.data(), 4}, {a.data() + 4, 4}, {a.data() + 8, 4}},
                           soa3<T>{{a.data() + 1, 4}, {a.data() + 5, 4}, {a.data() + 8, 4}}));
        CHECK(a == saved);
    } else if constexpr (Case == 3) {
        std::array<V, 5> a{};
        for (auto &v : a)
            v = V{T(1), T(2), T(3)};
        const auto saved = a;
        CHECK(!normalize_vectors(std::span<const V>{a.data(), 4}, std::span<V>{a.data() + 1, 4}));
        CHECK(a == saved);
        CHECK(normalize_vectors(std::span<const V>{a}, std::span<V>{a}));
        for (const auto &v : a)
            NEAR(length(v), T(1), extended::tolerance<T>);
        const auto normalized = a;
        const auto model = compose(trs<T>{V{T(5), T(6), T(7)}, quat<T>{}, V{T(2), T(3), T(4)}});
        CHECK(!transform_vectors(model, std::span<const V>{a.data(), 4},
                                 std::span<V>{a.data() + 1, 4}));
        CHECK(a == normalized);
        CHECK(!transform_vectors(model, std::span<const V>{a}, std::span<V>{a.data(), 4}));
        CHECK(a == normalized);
        CHECK(transform_vectors(model, std::span<const V>{a}, std::span<V>{a}));
        for (std::size_t i = 0; i < a.size(); ++i)
            NEAR(a[i], hadamard(normalized[i], V{T(2), T(3), T(4)}), extended::tolerance<T>);
        CHECK(!normalize_vectors(std::span<const V>{a}, std::span<V>{a.data(), 4}));
    } else if constexpr (Case == 4) {
        const auto a = mat<T, 3, 3>::identity();
        mat<T, 3, 2> b;
        b(1, 1) = nan;
        CHECK(!solve(a, b));
        b(1, 1) = inf;
        CHECK(!solve(a, b));
        b(1, 1) = T(1);
        CHECK(!solve(mat<T, 3, 3>{}, b));
        CHECK(!solve(a * std::numeric_limits<T>::denorm_min(), b));
        CHECK(b(1, 1) == T(1));
    } else if constexpr (Case == 5) {
        // Constant evaluation must not form an overflowing cached reciprocal.
        constexpr auto compile_time = [] {
            const auto m = mat<T, 2, 2>::identity() * (std::numeric_limits<T>::denorm_min() * T(8));
            const auto result = solve(m, m);
            return result && *result == mat<T, 2, 2>::identity();
        }();
        static_assert(compile_time);
        CHECK(compile_time);
        const T tiny = std::numeric_limits<T>::denorm_min() * T(8);
        const auto a = mat<T, 3, 3>::identity() * tiny;
        mat<T, 3, 2> x;
        x(0, 0) = T(2);
        x(1, 1) = T(-3);
        x(2, 0) = T(1);
        const auto result = solve(a, a * x);
        CHECK(result);
        NEAR(*result, x, extended::tolerance<T>);
    } else if constexpr (Case == 6) {
        CHECK(!prepared_ray<T>::create(ray<T, 3>{V{}, V{}}));
        CHECK(!prepared_ray<T>::create(ray<T, 3>{V{nan, T(0), T(0)}, V{T(1), T(0), T(0)}}));
        CHECK(!prepared_ray<T>::create(ray<T, 3>{V{}, V{inf, T(0), T(0)}}));
        ray<T, 3> original{V{}, V{T(1), T(0), T(0)}};
        const auto owned = prepared_ray<T>::create(original);
        CHECK(owned);
        original.origin = V{T(99), T(99), T(99)};
        original.direction = V{};
        CHECK(owned->source().origin == V{});
        CHECK((owned->source().direction == V{T(1), T(0), T(0)}));
    } else if constexpr (Case == 7) {
        const T large = std::numeric_limits<T>::max() / T(2);
        for (T scale : {large, std::numeric_limits<T>::denorm_min() * T(8)}) {
            const ray<T, 3> r{V{-scale, T(0), T(0)}, V{scale, T(0), T(0)}};
            const aabb<T, 3> b{V{scale / T(2), T(-1), T(-1)}, V{scale, T(1), T(1)}};
            const auto cached = prepared_ray<T>::create(r);
            CHECK(cached);
            const auto ref = intersect(r, b), v = intersect(*cached, b);
            CHECK(ref && v);
            NEAR(v->enter, ref->enter, extended::tolerance<T>);
            NEAR(v->exit, ref->exit, extended::tolerance<T>);
        }
    } else if constexpr (Case == 8) {
        const ray<T, 3> r{V{}, V{T(3), T(7), -T(0)}};
        const auto p = prepared_ray<T>::create(r);
        CHECK(p);
        const aabb<T, 3> b{V{T(3), T(7), T(0)}, V{T(6), T(7), T(0)}};
        const auto v = p->intersect(b);
        CHECK(v);
        NEAR(v->enter, T(1), T(0));
        NEAR(v->exit, T(1), T(0));
        CHECK(!p->intersect(b, T(2), T(1)));
        CHECK(!p->intersect(b, nan, inf));
    } else if constexpr (Case == 9) {
        const auto p = *prepared_ray<T>::create(ray<T, 3>{V{}, V{T(1), T(0), T(0)}});
        std::array<aabb<T, 3>, 2> b{};
        std::array<std::optional<ray_interval<T>>, 1> out{ray_interval<T>{T(11), T(12)}};
        CHECK(!intersect_many(p, std::span<const aabb<T, 3>>{b},
                              std::span<std::optional<ray_interval<T>>>{out}));
        CHECK(out[0]);
        NEAR(out[0]->enter, T(11), T(0));
        CHECK(intersect_many(p, std::span<const aabb<T, 3>>{},
                             std::span<std::optional<ray_interval<T>>>{}));
    } else if constexpr (Case == 10) {
        const frustum<T> f{};
        std::array<aabb<T, 3>, 2> b{};
        std::array<containment, 1> out{containment::inside};
        CHECK(!classify_batch(f, std::span<const aabb<T, 3>>{b}, std::span{out}));
        CHECK(out[0] == containment::inside);
        CHECK(classify_batch(f, std::span<const aabb<T, 3>>{}, std::span<containment>{}));
        const auto actual = *extract_frustum(*perspective(T(1), T(1.5), T(0.1), T(100)));
        const std::array<V, 3> points{V{T(0), T(0), T(-5)}, V{T(100), T(0), T(-1)},
                                      V{nan, T(0), T(0)}};
        const std::array<sphere<T>, 3> spheres{
            sphere<T>{points[0], T(1)}, sphere<T>{points[1], T(0.5)}, sphere<T>{points[2], T(1)}};
        std::array<containment, 3> classified{};
        CHECK(
            classify_batch(actual, std::span<const V>{points}, std::span<containment>{classified}));
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(classified[i] == classify(actual, points[i]));
        CHECK(classify_batch(actual, std::span<const sphere<T>>{spheres},
                             std::span<containment>{classified}));
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(classified[i] == classify(actual, spheres[i]));
    } else if constexpr (Case == 11) {
        bezier<T, 3, 3> b;
        b.control[2][1] = nan;
        CHECK(!cubic_polynomial<T, 3>::from_bezier(b));
        b.control[2][1] = inf;
        CHECK(!cubic_polynomial<T, 3>::from_bezier(b));
        b.control[2][1] = std::numeric_limits<T>::max();
        CHECK(!cubic_polynomial<T, 3>::from_bezier(b));
    } else if constexpr (Case == 12) {
        bezier<T, 3, 3> b{{V{T(1), T(-2), T(3)}, V{T(2), T(4), T(-1)}, V{T(-3), T(2), T(1)},
                           V{T(4), T(3), T(2)}}};
        const auto p = cubic_polynomial<T, 3>::from_bezier(b);
        CHECK(p);
        for (int i = 0; i <= 32; ++i) {
            const T t = T(i) / T(32);
            const auto s = p->evaluate_with_derivative(t);
            NEAR(s.position, b.evaluate(t), extended::tolerance<T>);
            NEAR(s.derivative, b.derivative().evaluate(t), extended::tolerance<T>);
            NEAR(p->derivative(t), s.derivative, extended::tolerance<T>);
            NEAR(p->second_derivative(t),
                 (b.control[2] - b.control[1] * T(2) + b.control[0]) * (T(6) * (T(1) - t)) +
                     (b.control[3] - b.control[2] * T(2) + b.control[1]) * (T(6) * t),
                 extended::tolerance<T>);
        }
    } else {
        const bezier<T, 3, 0> b{{V{T(2), T(-3), T(4)}}};
        for (T t : {T(-3), T(0), T(1), T(2)}) {
            const auto s = b.evaluate_with_derivative(t);
            CHECK(s.position == b.control[0]);
            CHECK(s.derivative == V{});
        }
    }
}
template <class T, int Mode> void projection() {
    const auto hand = (Mode & 1) ? handedness::left : handedness::right;
    const auto range = (Mode & 2) ? depth_range::minus_one_to_one : depth_range::zero_to_one;
    const auto dir = (Mode & 4) ? depth_direction::reverse : depth_direction::forward;
    const T sign = (Mode & 1) ? T(1) : T(-1);
    const viewport<T> vp{T(17), T(-23), T(1920), T(1080), T(0.2), T(0.9)};
    const auto m = perspective(T(1.1), T(16) / T(9), T(0.5), T(500), hand, range, dir);
    CHECK(m);
    const auto inv = inverse(*m);
    CHECK(inv);
    for (T z : {T(0.5), T(3), T(70), T(500)}) {
        const vec<T, 3> p{T(0.1), T(-0.2), sign * z};
        const auto screen = project(p, *m, vp, range);
        CHECK(screen);
        const auto back = unproject(*screen, *inv, vp, range);
        CHECK(back);
        NEAR(*back, p, (std::same_as<T, float> ? T(8e-4) : T(3e-11)));
        if (z == T(0.5) || z == T(500)) {
            const bool near = z == T(0.5);
            const T expected = (near != bool(Mode & 4)) ? vp.min_depth : vp.max_depth;
            NEAR((*screen)[2], expected, extended::tolerance<T>);
        }
    }
}
} // namespace extended_contract
