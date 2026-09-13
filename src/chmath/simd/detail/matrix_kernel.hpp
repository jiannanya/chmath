#pragma once
#include <cstddef>

#ifndef CHMATH_ENABLE_SIMD
#define CHMATH_ENABLE_SIMD 1
#endif
#if CHMATH_ENABLE_SIMD &&                                                                          \
    (defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2))
#define CHMATH_DETAIL_MATRIX_SSE2 1
#include <emmintrin.h>
#elif CHMATH_ENABLE_SIMD && (defined(__ARM_NEON) || defined(__ARM_NEON__))
#define CHMATH_DETAIL_MATRIX_NEON 1
#include <arm_neon.h>
#endif

namespace chm::detail {
#if defined(CHMATH_DETAIL_MATRIX_SSE2) || defined(CHMATH_DETAIL_MATRIX_NEON)
inline constexpr bool native_mat4_kernel = true;
#else
inline constexpr bool native_mat4_kernel = false;
#endif

// Internal kernel: column-major, naturally aligned float storage. Output must
// not overlap either input. The public value-returning operator ensures this.
inline void multiply_mat4_float(const float *a, const float *b, float *out) noexcept {
#if defined(CHMATH_DETAIL_MATRIX_SSE2)
    const auto a0 = _mm_loadu_ps(a), a1 = _mm_loadu_ps(a + 4), a2 = _mm_loadu_ps(a + 8),
               a3 = _mm_loadu_ps(a + 12);
    for (std::size_t j = 0; j < 4; ++j) {
        const auto column = _mm_loadu_ps(b + 4 * j);
        auto sum =
            _mm_add_ps(_mm_setzero_ps(),
                       _mm_mul_ps(a0, _mm_shuffle_ps(column, column, _MM_SHUFFLE(0, 0, 0, 0))));
        sum = _mm_add_ps(sum,
                         _mm_mul_ps(a1, _mm_shuffle_ps(column, column, _MM_SHUFFLE(1, 1, 1, 1))));
        sum = _mm_add_ps(sum,
                         _mm_mul_ps(a2, _mm_shuffle_ps(column, column, _MM_SHUFFLE(2, 2, 2, 2))));
        sum = _mm_add_ps(sum,
                         _mm_mul_ps(a3, _mm_shuffle_ps(column, column, _MM_SHUFFLE(3, 3, 3, 3))));
        _mm_storeu_ps(out + 4 * j, sum);
    }
#elif defined(CHMATH_DETAIL_MATRIX_NEON)
    const auto a0 = vld1q_f32(a), a1 = vld1q_f32(a + 4), a2 = vld1q_f32(a + 8),
               a3 = vld1q_f32(a + 12);
    for (std::size_t j = 0; j < 4; ++j) {
        auto sum = vaddq_f32(vdupq_n_f32(0), vmulq_f32(a0, vdupq_n_f32(b[4 * j])));
        sum = vaddq_f32(sum, vmulq_f32(a1, vdupq_n_f32(b[4 * j + 1])));
        sum = vaddq_f32(sum, vmulq_f32(a2, vdupq_n_f32(b[4 * j + 2])));
        sum = vaddq_f32(sum, vmulq_f32(a3, vdupq_n_f32(b[4 * j + 3])));
        vst1q_f32(out + 4 * j, sum);
    }
#else
    for (std::size_t j = 0; j < 4; ++j)
        for (std::size_t i = 0; i < 4; ++i) {
            float sum{};
            for (std::size_t k = 0; k < 4; ++k)
                sum += a[4 * k + i] * b[4 * j + k];
            out[4 * j + i] = sum;
        }
#endif
}
} // namespace chm::detail
#undef CHMATH_DETAIL_MATRIX_SSE2
#undef CHMATH_DETAIL_MATRIX_NEON
