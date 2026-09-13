#pragma once
#include <cstdint>
#include <span>

namespace chm::detail {
template <class A, class B>
[[nodiscard]] inline bool spans_overlap(std::span<A> a, std::span<B> b) noexcept {
    if (a.empty() || b.empty())
        return false;
    const auto x = reinterpret_cast<std::uintptr_t>(a.data());
    const auto y = reinterpret_cast<std::uintptr_t>(b.data());
    return x <= y ? y - x < a.size_bytes() : x - y < b.size_bytes();
}
template <class A, class B>
[[nodiscard]] inline bool unsafe_overlap(std::span<A> a, std::span<B> b) noexcept {
    return spans_overlap(a, b) &&
           !(static_cast<const void *>(a.data()) == static_cast<const void *>(b.data()) &&
             a.size_bytes() == b.size_bytes());
}
} // namespace chm::detail
