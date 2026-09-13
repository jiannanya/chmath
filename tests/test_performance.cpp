#include "performance_helpers.hpp"
using namespace chm;
namespace perf {
template <std::size_t N>
CH_TEST_NOINLINE double inverse_kernel(std::span<const mat<double, N, N>> in,
                                       std::span<mat<double, N, N>> out) {
    for (std::size_t i = 0; i < in.size(); ++i) {
        const auto r = inverse(in[i]);
        if (!r)
            return boundary::nan;
        out[i] = *r;
    }
    return out[in.size() / 2](0, 0);
}
template <std::size_t N> void inverse_case() {
    constexpr std::size_t n = 128;
    std::vector<mat<double, N, N>> in(n), out(n);
    for (auto &m : in) {
        m = boundary::well_conditioned<double, N>();
        for (std::size_t j = 0; j < N; ++j)
            m(j, j) += test::sample<double>(0, 1);
    }
    time("LU inverse", n, [&] { return inverse_kernel<N>(in, out); });
    for (std::size_t i = 0; i < n; ++i)
        NEAR(in[i] * out[i], (mat<double, N, N>::identity()), 1e-12);
}
CH_TEST_NOINLINE inline double chain_kernel(std::span<const mat4f> in, mat4f &out) {
    out = mat4f::identity();
    for (const auto &m : in)
        out = out * m;
    return out(0, 0);
}
CH_TEST_NOINLINE inline double reused_solve(const lu_factorization<double, 4> &lu,
                                            std::span<const vec4d> rhs, std::span<vec4d> out) {
    for (std::size_t i = 0; i < rhs.size(); ++i)
        out[i] = lu.solve(rhs[i]);
    return out[0][0];
}
CH_TEST_NOINLINE inline double repeated_solve(const mat4d &a, std::span<const vec4d> rhs,
                                              std::span<vec4d> out) {
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const auto x = solve(a, rhs[i]);
        if (!x)
            return boundary::nan;
        out[i] = *x;
    }
    return out[0][0];
}
CH_TEST_NOINLINE inline double rotate_kernel(const quatf &q, std::span<const vec3f> in,
                                             std::span<vec3f> out) {
    for (std::size_t i = 0; i < in.size(); ++i)
        out[i] = rotate(q, in[i]);
    return out[in.size() / 2][0];
}
template <class Curve>
CH_TEST_NOINLINE double evaluate_curve(const Curve &c, std::span<const double> parameters,
                                       std::span<vec3d> out) {
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        if constexpr (requires { c.evaluate(parameters[i]).has_value(); }) {
            const auto p = c.evaluate(parameters[i]);
            if (!p)
                return boundary::nan;
            out[i] = *p;
        } else
            out[i] = c.evaluate(parameters[i]);
    }
    return out[parameters.size() / 2][0];
}
CH_TEST_NOINLINE inline double eigen_kernel(const mat<double, 6, 6> &a,
                                            symmetric_eigen_result<double, 6> &out) {
    const auto e = symmetric_eigen(a, 1e-12, 64);
    if (!e)
        return boundary::nan;
    out = *e;
    return out.values[0];
}
CH_TEST_NOINLINE inline double integration_kernel(integration_result<double> &out) {
    const auto r =
        integrate([](double x) { return std::sin(x); }, 0.0, pi<double>, 1e-10, 20, 4097);
    if (!r)
        return boundary::nan;
    out = *r;
    return out.value;
}
} // namespace perf
CH_TEST(perf_matrix_two_float) {
    perf::matrix_case<float, 2>();
}
CH_TEST(perf_matrix_three_double) {
    perf::matrix_case<double, 3>();
}
CH_TEST(perf_matrix_four_float) {
    perf::matrix_case<float, 4>();
}
CH_TEST(perf_matrix_eight_double) {
    perf::matrix_case<double, 8>();
}
CH_TEST(perf_matrix_transform_chain) {
    const mat4f step = to_matrix(translation(vec3f{0.01f, 0.02f, 0.03f}));
    const std::vector<mat4f> in(128, step);
    mat4f out;
    perf::time("mat4 chain", in.size(), [&] { return perf::chain_kernel(in, out); });
    NEAR(out, to_matrix(translation(vec3f{1.28f, 2.56f, 3.84f})), 1e-5f);
}
CH_TEST(perf_inverse_three_double) {
    perf::inverse_case<3>();
}
CH_TEST(perf_inverse_four_double) {
    perf::inverse_case<4>();
}
CH_TEST(perf_lu_repeated_right_hand_sides) {
    const auto a = boundary::well_conditioned<double, 4>();
    const auto lu = factor_lu(a);
    CHECK(lu);
    std::vector<vec4d> rhs(128), out(128), reference(128);
    for (auto &v : rhs)
        for (auto &x : v.elements)
            x = test::sample<double>();
    const double rt = perf::time("refactor each RHS", rhs.size(),
                                 [&] { return perf::repeated_solve(a, rhs, reference); });
    const double at = perf::time("reuse LU factors", rhs.size(),
                                 [&] { return perf::reused_solve(*lu, rhs, out); });
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        NEAR(a * out[i], rhs[i], 1e-12);
        NEAR(out[i], reference[i], 1e-12);
    }
    perf::relative(at, rt, 2.0);
}
CH_TEST(perf_normalize_float) {
    perf::normalize_case<float>();
}
CH_TEST(perf_normalize_double) {
    perf::normalize_case<double>();
}
CH_TEST(perf_normalize_extreme_dynamic_range) {
    perf::normalize_case<double>(true);
}
CH_TEST(perf_quaternion_rotate) {
    std::vector<vec3f> in(1024), out(1024);
    for (auto &v : in)
        v = {test::sample<float>(), test::sample<float>(), test::sample<float>()};
    const auto q = from_euler_xyz(vec3f{0.1f, 0.2f, 0.3f});
    const auto matrix = to_matrix(q);
    perf::time("quaternion rotate", in.size(), [&] { return perf::rotate_kernel(q, in, out); });
    for (std::size_t i = 0; i < in.size(); ++i) {
        NEAR(out[i], matrix * in[i], 1e-5f);
        NEAR(length(out[i]), length(in[i]), 1e-5f);
    }
}
CH_TEST(perf_aos_small_batch) {
    perf::aos_case(17);
}
CH_TEST(perf_aos_large_batch) {
    perf::aos_case(65536);
}
CH_TEST(perf_soa_small_batch) {
    perf::soa_case(17);
}
CH_TEST(perf_soa_large_batch) {
    perf::soa_case(65536);
}
CH_TEST(perf_dot_small_batch) {
    perf::dot_case(17);
}
CH_TEST(perf_dot_large_batch) {
    perf::dot_case(65536);
}
CH_TEST(perf_bezier_cubic_evaluation) {
    const bezier<double, 3, 3> c{{vec3d{0, 0, 0}, vec3d{1, 2, 0}, vec3d{2, 2, 1}, vec3d{3, 0, 1}}};
    std::vector<double> t(256);
    std::vector<vec3d> out(256);
    for (auto &v : t)
        v = test::sample<double>(0, 1);
    perf::time("cubic Bezier", t.size(), [&] { return perf::evaluate_curve(c, t, out); });
    for (std::size_t i = 0; i < t.size(); ++i) {
        const double u = t[i], v = 1 - u;
        NEAR(out[i],
             c.control[0] * (v * v * v) + c.control[1] * (3 * u * v * v) +
                 c.control[2] * (3 * u * u * v) + c.control[3] * (u * u * u),
             1e-12);
    }
}
CH_TEST(perf_bspline_degree_eight_evaluation) {
    bezier<double, 3, 8> polynomial;
    for (auto &p : polynomial.control)
        p = {test::sample<double>(), test::sample<double>(), test::sample<double>()};
    std::array<double, 18> knots{};
    std::fill(knots.begin() + 9, knots.end(), 1.0);
    const auto c = bspline_view<double, 3, 8>::create(polynomial.control, knots);
    CHECK(c);
    std::vector<double> t(128);
    std::vector<vec3d> out(128);
    for (auto &v : t)
        v = test::sample<double>(0, 1);
    perf::time("degree8 B-spline", t.size(), [&] { return perf::evaluate_curve(*c, t, out); });
    for (std::size_t i = 0; i < t.size(); ++i)
        NEAR(out[i], polynomial.evaluate(t[i]), 1e-12);
}
CH_TEST(perf_eigen_convergence_budget) {
    auto a = boundary::well_conditioned<double, 6>();
    symmetric_eigen_result<double, 6> e;
    perf::time("Jacobi eigen6", 1, [&] { return perf::eigen_kernel(a, e); });
    CHECK(e.converged);
    CHECK(e.sweeps <= 64);
    for (std::size_t i = 0; i < 6; ++i)
        NEAR(a * e.vectors.column(i), e.vectors.column(i) * e.values[i], 1e-10);
}
CH_TEST(perf_integration_evaluation_budget) {
    integration_result<double> r;
    perf::time("adaptive sine integral", 1, [&] { return perf::integration_kernel(r); });
    CHECK(r.converged);
    CHECK(r.evaluations <= 4097);
    NEAR(r.value, 2.0, 1e-10);
    std::size_t calls = 0;
    const auto limited = integrate(
        [&](double x) {
            ++calls;
            return std::exp(x);
        },
        0.0, 1.0, 1e-15, 32, 33);
    CHECK(limited);
    CHECK(calls == limited->evaluations && calls <= 33);
    CHECK(!limited->converged);
}
