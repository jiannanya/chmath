#pragma once
#include "../guarded_buffer.hpp"
#include "common.hpp"
namespace extended_batch {
using namespace chm;
template <class T, int Op>
bool invoke(const affine3<T> &a, const quat<T> &q, const_soa3<T> x, const_soa3<T> y, soa3<T> out) {
    if constexpr (Op == 0)
        return transform_points(a, x, out);
    else if constexpr (Op == 1)
        return transform_vectors(a, x, out);
    else if constexpr (Op == 2)
        return rotate_vectors(q, x, out);
    else if constexpr (Op == 3)
        return dot_batch(x, y, out.x);
    else if constexpr (Op == 4)
        return cross_batch(x, y, out);
    else
        return normalize_vectors(x, out);
}
template <class T, int Op>
vec<T, 3> oracle(const affine3<T> &a, const quat<T> &q, const vec<T, 3> &x, const vec<T, 3> &y) {
    if constexpr (Op == 0)
        return transform_point(a, x);
    else if constexpr (Op == 1)
        return transform_vector(a, x);
    else if constexpr (Op == 2)
        return rotate(q, x);
    else if constexpr (Op == 3)
        return {dot(x, y), T(0), T(0)};
    else if constexpr (Op == 4)
        return cross(x, y);
    else
        return extended::normalization_reference(x);
}
template <class T, int Op, std::size_t N> void guarded_contract() {
    const auto q = *from_axis_angle(vec<T, 3>{T(1), T(2), T(3)}, T(0.7));
    auto a = rotation(q);
    a.matrix(0, 3) = T(3);
    a.matrix(1, 3) = T(-2);
    if constexpr (N == 0) {
        CHECK((invoke<T, Op>(a, q, const_soa3<T>{}, const_soa3<T>{}, soa3<T>{})));
    } else {
        for (bool end : {false, true})
            for (int alias : {0, 1, 2}) {
                safety::guarded_buffer<T> x(N, end), y(N, end), z(N, end), bx(N, end), by(N, end),
                    bz(N, end), ox(N, end), oy(N, end), oz(N, end);
                std::array<vec<T, 3>, N> expected{};
                for (std::size_t i = 0; i < N; ++i) {
                    x[i] = T(i + 1) * T(0.25);
                    y[i] = T(-2);
                    z[i] = T(i % 3) * T(0.5);
                    bx[i] = T(-1);
                    by[i] = T(i % 5) * T(0.2);
                    bz[i] = T(3);
                    expected[i] = oracle<T, Op>(a, q, {x[i], y[i], z[i]}, {bx[i], by[i], bz[i]});
                    ox[i] = oy[i] = oz[i] = T(-99);
                }
                bx.read_only();
                by.read_only();
                bz.read_only();
                if (alias == 0) {
                    x.read_only();
                    y.read_only();
                    z.read_only();
                }
                soa3<T> out = alias == 0   ? soa3<T>{ox.span(), oy.span(), oz.span()}
                              : alias == 1 ? soa3<T>{x.span(), y.span(), z.span()}
                                           : soa3<T>{y.span(), z.span(), x.span()};
                CHECK((invoke<T, Op>(a, q, {x.const_span(), y.const_span(), z.const_span()},
                                     {bx.const_span(), by.const_span(), bz.const_span()}, out)));
                for (std::size_t i = 0; i < N; ++i) {
                    NEAR(out.x[i], expected[i][0], extended::tolerance<T>);
                    if constexpr (Op != 3) {
                        NEAR(out.y[i], expected[i][1], extended::tolerance<T>);
                        NEAR(out.z[i], expected[i][2], extended::tolerance<T>);
                    }
                    NEAR(bx[i], T(-1), T(0));
                    if (alias == 0) {
                        NEAR(x[i], T(i + 1) * T(0.25), T(0));
                        NEAR(y[i], T(-2), T(0));
                    }
                }
            }
    }
}
} // namespace extended_batch
