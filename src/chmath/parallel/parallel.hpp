#pragma once

// Portable, opt-in parallel execution for the batch entry points.
//
// The library is header only and has no threading dependency by default. Set
// CHMATH_ENABLE_PARALLEL=1 (CMake option CHMATH_ENABLE_PARALLEL) to compile the
// worker spawning path; otherwise every parallel entry point runs serially with
// identical results, so calling code does not need to change. Only standard
// C++ is used - no platform threads API, intrinsics or inline assembly.
//
// The batch kernels are pure functions over caller-provided, disjoint buffers,
// so splitting the index range across workers never introduces a data race.
// Splits are contiguous and start on element boundaries; for elements at most
// one cache line in size they also start on different cache lines, which keeps
// adjacent workers from sharing a line (see detail::aligned_split).
//
// Parallel execution is opt-in per call: the worker-count overloads of the
// batch entry points are the only functions that can start threads. The
// default (single-argument) entry points always stay serial.

#include <algorithm>
#include <cstddef>

#ifndef CHMATH_ENABLE_PARALLEL
#define CHMATH_ENABLE_PARALLEL 0
#endif

#if CHMATH_ENABLE_PARALLEL
#include <exception>
#include <thread>
#include <vector>
#endif

namespace chm {

// Number of hardware threads reported by the platform, at least one.
[[nodiscard]] inline unsigned parallel_hardware_threads() noexcept {
#if CHMATH_ENABLE_PARALLEL
    const unsigned reported = std::thread::hardware_concurrency();
    return reported > 0U ? reported : 1U;
#else
    return 1U;
#endif
}

// Worker count used when a parallel entry point is called without an explicit
// count. Falls back to one on builds without CHMATH_ENABLE_PARALLEL.
[[nodiscard]] inline unsigned parallel_default_workers() noexcept {
    return parallel_hardware_threads();
}

// Smallest batch length that is worth splitting. Below this the request runs on
// the calling thread, which avoids paying thread creation for short batches.
// The value is deliberately conservative: spawning a handful of threads costs
// on the order of a hundred microseconds on desktop hardware, so a split only
// pays off once the batch contains enough items for the per-item work to
// dominate that cost. Callers with expensive per-item work can split smaller
// batches by calling parallel_for directly.
inline constexpr std::size_t parallel_min_chunk = 1U << 15;

namespace detail {

// Cache line size used for split alignment. std::hardware_destructive_
// interference_size is preferred when the standard library provides it, so
// splits follow the target's line size (64 bytes on most x86 and many ARM
// cores, 128 on Apple Silicon). The value only tunes where splits fall; it is
// never part of the public ABI.
#if defined(__cpp_lib_hardware_interference_size)
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winterference-size"
#endif
inline constexpr std::size_t parallel_line_bytes = std::hardware_destructive_interference_size;
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic pop
#endif
#else
inline constexpr std::size_t parallel_line_bytes = 64;
#endif

// Split [0, count) into `workers` disjoint contiguous ranges. Interior split
// points advance by whole cache lines' worth of `item_bytes` sized elements, so
// adjacent chunks start on different lines (and can never share one) for any
// element at most one line in size; larger elements are split on element
// boundaries because an element already spans multiple lines. The ranges still
// cover [0, count) exactly once: every chunk i receives
// [split(i), split(i + 1)).
[[nodiscard]] inline std::size_t aligned_split(std::size_t count, unsigned workers,
                                               std::size_t item_bytes,
                                               unsigned index) noexcept {
    if (index == 0U)
        return 0;
    if (index >= workers)
        return count;
    if (item_bytes == 0 || count == 0)
        return count * index / workers;
    const std::size_t per_line =
        item_bytes >= parallel_line_bytes
            ? 1
            : (parallel_line_bytes + item_bytes - 1) / item_bytes;
    std::size_t split = count * index / workers;
    const std::size_t remainder = split % per_line;
    if (remainder != 0)
        split += per_line - remainder;
    return std::min(split, count);
}

// Clamp a requested worker count to something usable for `count` items.
[[nodiscard]] inline unsigned effective_workers(std::size_t count, unsigned requested) noexcept {
    if (requested <= 1U || count < parallel_min_chunk)
        return 1U;
    const unsigned hardware = parallel_hardware_threads();
    const unsigned limited = std::min(requested, hardware);
    const std::size_t by_items = count / parallel_min_chunk;
    const unsigned usable =
        static_cast<unsigned>(std::min<std::size_t>(limited, std::max<std::size_t>(by_items, 1)));
    return usable > 1U ? usable : 1U;
}

} // namespace detail

// Number of workers a request would actually use. An empty range uses zero,
// matching the return value of parallel_for. Exposed for callers that want to
// size their own buffers or log the decision, and used by tests.
[[nodiscard]] inline unsigned parallel_workers(std::size_t count,
                                               unsigned requested = 0) noexcept {
    if (count == 0)
        return 0;
    return detail::effective_workers(count,
                                     requested == 0U ? parallel_default_workers() : requested);
}

// Run body(begin, end) over a set of disjoint contiguous sub-ranges of
// [0, count) and return the number of workers that ran (1 on a serial call).
//
// `body` is called once per worker with the half-open range the worker owns. It
// is called on the calling thread for the first range and on worker threads for
// the rest; it must not throw in builds that translate to a `noexcept` entry
// point. When exceptions are enabled anyway, an exception escaping a worker is
// captured, every other worker is joined, and the first failure is rethrown on
// the calling thread.
//
// If the platform refuses to create a worker thread, the unclaimed ranges run
// on the calling thread and the function returns the reduced worker count
// instead of failing. item_bytes is the size of one element, used to keep
// adjacent ranges off a shared cache line.
template <class Body>
unsigned parallel_for(std::size_t count, unsigned workers, std::size_t item_bytes, Body &&body) {
    if (count == 0)
        return 0;
    const unsigned active = detail::effective_workers(count, workers);
    if (active <= 1U) {
        body(std::size_t{0}, count);
        return 1;
    }
#if CHMATH_ENABLE_PARALLEL
    std::vector<std::thread> pool;
    std::vector<std::exception_ptr> failures(active);
    pool.reserve(active - 1U);
    unsigned spawned = 1U;
    try {
        for (unsigned i = 1U; i < active; ++i) {
            pool.emplace_back([&, i] {
                const std::size_t begin = detail::aligned_split(count, active, item_bytes, i);
                const std::size_t end = detail::aligned_split(count, active, item_bytes, i + 1U);
                try {
                    body(begin, end);
                } catch (...) {
                    failures[i] = std::current_exception();
                }
            });
            ++spawned;
        }
    } catch (...) {
        // Out of thread resources: keep the workers that started and run the
        // remaining ranges on the calling thread below.
    }
    try {
        body(detail::aligned_split(count, active, item_bytes, 0U),
             detail::aligned_split(count, active, item_bytes, 1U));
    } catch (...) {
        failures[0] = std::current_exception();
    }
    if (spawned < active) {
        try {
            body(detail::aligned_split(count, active, item_bytes, spawned),
                 detail::aligned_split(count, active, item_bytes, active));
        } catch (...) {
            if (!failures[0])
                failures[0] = std::current_exception();
        }
    }
    for (std::thread &worker : pool)
        worker.join();
    for (const std::exception_ptr &failure : failures)
        if (failure)
            std::rethrow_exception(failure);
    return spawned;
#else
    (void)item_bytes;
    body(std::size_t{0}, count);
    return 1;
#endif
}

// Convenience overload for element type T: parallel_for<T>(count, workers, body).
template <class T, class Body>
unsigned parallel_for(std::size_t count, unsigned workers, Body &&body) {
    return parallel_for(count, workers, sizeof(T), static_cast<Body &&>(body));
}

} // namespace chm
