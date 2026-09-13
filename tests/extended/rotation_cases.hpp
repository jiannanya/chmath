#pragma once
#include "common.hpp"
namespace extended_rotation {
using namespace chm;
template <class T, int Axis> vec<T, 3> axis() {
    if constexpr (Axis == 3)
        return normalize(vec<T, 3>{T(1), T(2), T(-3)});
    else {
        vec<T, 3> v;
        v[Axis] = T(1);
        return v;
    }
}
template <class T, int Angle> T angle() {
    if constexpr (Angle == 0)
        return T(0);
    else if constexpr (Angle == 1)
        return T(1e-5);
    else if constexpr (Angle == 2)
        return pi<T> / T(2);
    else
        return pi<T> - T(1e-3);
}
template <class T, int Axis, int Angle> void rodrigues_reference() {
    const auto a = axis<T, Axis>();
    const T t = angle<T, Angle>();
    const auto q = from_axis_angle(a, t);
    CHECK(q);
    const vec<T, 3> v{T(0.3), T(-0.2), T(1)};
    const auto oracle =
        v * std::cos(t) + cross(a, v) * std::sin(t) + a * (dot(a, v) * (T(1) - std::cos(t)));
    NEAR(rotate(*q, v), oracle, extended::tolerance<T>);
    const auto matrix = to_matrix(*q);
    NEAR(matrix * v, oracle, extended::tolerance<T>);
    const auto restored = from_matrix(matrix);
    CHECK(restored);
    NEAR(rotate(*restored, v), oracle, extended::tolerance<T>);
    NEAR(transpose(matrix) * matrix, (mat<T, 3, 3>::identity()), extended::tolerance<T>);
}
template <class T, int Axis, int Angle> void between_directions() {
    const auto q = from_axis_angle(axis<T, Axis>(), angle<T, Angle>());
    CHECK(q);
    const vec<T, 3> v{T(0.2), T(-0.3), T(0.7)};
    const auto target = rotate(*q, v);
    const auto between = rotation_between(v, target);
    CHECK(between);
    NEAR(rotate(*between, normalize(v)), normalize(target), extended::tolerance<T>);
    const auto inv = inverse(*between);
    CHECK(inv);
    NEAR(rotate(*inv, normalize(target)), normalize(v), extended::tolerance<T>);
}
template <class T, int Axis, int Angle> void interpolation_arc() {
    const auto a = axis<T, Axis>();
    const T theta = angle<T, Angle>();
    const auto q = from_axis_angle(a, theta);
    CHECK(q);
    for (T t : {T(0), T(0.125), T(0.5), T(1)}) {
        const auto value = slerp(quat<T>{}, *q, t);
        const auto expected = from_axis_angle(a, theta * t);
        CHECK(expected);
        NEAR(to_matrix(value), to_matrix(*expected), extended::tolerance<T>);
        NEAR(dot(value, value), T(1), extended::tolerance<T>);
        NEAR(to_matrix(slerp(quat<T>{}, -*q, t)), to_matrix(value), extended::tolerance<T>);
    }
}
} // namespace extended_rotation
