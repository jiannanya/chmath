#pragma once
#include <chmath/core/scalar.hpp>
#include <cstdint>
#ifndef CHMATH_ENABLE_SIMD
#define CHMATH_ENABLE_SIMD 1
#endif
#if CHMATH_ENABLE_SIMD &&                                                                          \
    (defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2))
#define CHMATH_PACKET_SSE2 1
#include <emmintrin.h>
#elif CHMATH_ENABLE_SIMD && (defined(__ARM_NEON) || defined(__ARM_NEON__))
#define CHMATH_PACKET_NEON 1
#include <arm_neon.h>
#endif

namespace chm::detail {
template <floating T> struct packet {
    using value = T;
    static constexpr std::size_t lanes = 1;
    static constexpr bool vector_sqrt = false;
    static constexpr const char *name = "scalar";
    static value load(const T *p) noexcept { return *p; }
    static void store(T *p, value v) noexcept { *p = v; }
    static value splat(T x) noexcept { return x; }
    static value add(value a, value b) noexcept { return a + b; }
    static value sub(value a, value b) noexcept { return a - b; }
    static value mul(value a, value b) noexcept { return a * b; }
    static value div(value a, value b) noexcept { return a / b; }
    static value sqrt(value a) noexcept { return std::sqrt(a); }
    static bool normal_range(value a) noexcept {
        return a >= std::numeric_limits<T>::min() && is_finite(a);
    }
};
#if defined(CHMATH_PACKET_SSE2)
template <> struct packet<float> {
    using value = __m128;
    static constexpr std::size_t lanes = 4;
    static constexpr bool vector_sqrt = true;
    static constexpr const char *name = "SSE2";
    static value load(const float *p) noexcept { return _mm_loadu_ps(p); }
    static void store(float *p, value v) noexcept { _mm_storeu_ps(p, v); }
    static value splat(float x) noexcept { return _mm_set1_ps(x); }
    static value add(value a, value b) noexcept { return _mm_add_ps(a, b); }
    static value sub(value a, value b) noexcept { return _mm_sub_ps(a, b); }
    static value mul(value a, value b) noexcept { return _mm_mul_ps(a, b); }
    static value div(value a, value b) noexcept { return _mm_div_ps(a, b); }
    static value sqrt(value a) noexcept { return _mm_sqrt_ps(a); }
    static bool normal_range(value a) noexcept {
        return _mm_movemask_ps(
                   _mm_and_ps(_mm_cmpge_ps(a, splat(std::numeric_limits<float>::min())),
                              _mm_cmple_ps(a, splat(std::numeric_limits<float>::max())))) == 15;
    }
};
template <> struct packet<double> {
    using value = __m128d;
    static constexpr std::size_t lanes = 2;
    static constexpr bool vector_sqrt = true;
    static constexpr const char *name = "SSE2";
    static value load(const double *p) noexcept { return _mm_loadu_pd(p); }
    static void store(double *p, value v) noexcept { _mm_storeu_pd(p, v); }
    static value splat(double x) noexcept { return _mm_set1_pd(x); }
    static value add(value a, value b) noexcept { return _mm_add_pd(a, b); }
    static value sub(value a, value b) noexcept { return _mm_sub_pd(a, b); }
    static value mul(value a, value b) noexcept { return _mm_mul_pd(a, b); }
    static value div(value a, value b) noexcept { return _mm_div_pd(a, b); }
    static value sqrt(value a) noexcept { return _mm_sqrt_pd(a); }
    static bool normal_range(value a) noexcept {
        return _mm_movemask_pd(
                   _mm_and_pd(_mm_cmpge_pd(a, splat(std::numeric_limits<double>::min())),
                              _mm_cmple_pd(a, splat(std::numeric_limits<double>::max())))) == 3;
    }
};
#elif defined(CHMATH_PACKET_NEON)
template <> struct packet<float> {
    using value = float32x4_t;
    static constexpr std::size_t lanes = 4;
    static constexpr const char *name = "NEON";
#if defined(__aarch64__) || defined(_M_ARM64)
    static constexpr bool vector_sqrt = true;
    static value div(value a, value b) noexcept { return vdivq_f32(a, b); }
    static value sqrt(value a) noexcept { return vsqrtq_f32(a); }
#else
    static constexpr bool vector_sqrt = false;
    static value div(value a, value b) noexcept {
        float x[4], y[4];
        vst1q_f32(x, a);
        vst1q_f32(y, b);
        for (std::size_t i = 0; i < 4; ++i)
            x[i] /= y[i];
        return vld1q_f32(x);
    }
    static value sqrt(value a) noexcept {
        float x[4];
        vst1q_f32(x, a);
        for (float &v : x)
            v = std::sqrt(v);
        return vld1q_f32(x);
    }
#endif
    static value load(const float *p) noexcept { return vld1q_f32(p); }
    static void store(float *p, value v) noexcept { vst1q_f32(p, v); }
    static value splat(float x) noexcept { return vdupq_n_f32(x); }
    static value add(value a, value b) noexcept { return vaddq_f32(a, b); }
    static value sub(value a, value b) noexcept { return vsubq_f32(a, b); }
    static value mul(value a, value b) noexcept { return vmulq_f32(a, b); }
    static bool normal_range(value a) noexcept {
        const auto m = vandq_u32(vcgeq_f32(a, splat(std::numeric_limits<float>::min())),
                                 vcleq_f32(a, splat(std::numeric_limits<float>::max())));
        return (vgetq_lane_u32(m, 0) & vgetq_lane_u32(m, 1) & vgetq_lane_u32(m, 2) &
                vgetq_lane_u32(m, 3)) == ~std::uint32_t(0);
    }
};
#if defined(__aarch64__) || defined(_M_ARM64)
template <> struct packet<double> {
    using value = float64x2_t;
    static constexpr std::size_t lanes = 2;
    static constexpr bool vector_sqrt = true;
    static constexpr const char *name = "NEON";
    static value load(const double *p) noexcept { return vld1q_f64(p); }
    static void store(double *p, value v) noexcept { vst1q_f64(p, v); }
    static value splat(double x) noexcept { return vdupq_n_f64(x); }
    static value add(value a, value b) noexcept { return vaddq_f64(a, b); }
    static value sub(value a, value b) noexcept { return vsubq_f64(a, b); }
    static value mul(value a, value b) noexcept { return vmulq_f64(a, b); }
    static value div(value a, value b) noexcept { return vdivq_f64(a, b); }
    static value sqrt(value a) noexcept { return vsqrtq_f64(a); }
    static bool normal_range(value a) noexcept {
        const auto m = vandq_u64(vcgeq_f64(a, splat(std::numeric_limits<double>::min())),
                                 vcleq_f64(a, splat(std::numeric_limits<double>::max())));
        return (vgetq_lane_u64(m, 0) & vgetq_lane_u64(m, 1)) == ~std::uint64_t(0);
    }
};
#endif
#endif
} // namespace chm::detail
#undef CHMATH_PACKET_SSE2
#undef CHMATH_PACKET_NEON
