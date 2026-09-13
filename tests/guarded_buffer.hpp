#pragma once
#include <cstddef>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#pragma push_macro("NEAR")
#pragma push_macro("near")
#undef NEAR
#undef near
#include <windows.h>
#pragma pop_macro("near")
#pragma pop_macro("NEAR")
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace safety {
// Two inaccessible pages surround a writable region. Place the logical buffer
// against either boundary so even a single SIMD lane crossing it faults.
template <class T> class guarded_buffer {
    static_assert(std::is_trivially_destructible_v<T>);
    void *allocation_{};
    std::size_t allocation_bytes_{}, page_bytes_{}, writable_bytes_{}, count_{};
    T *data_{};

  public:
    explicit guarded_buffer(std::size_t count, bool at_end = true) : count_(count) {
#if defined(_WIN32)
        SYSTEM_INFO info{};
        GetSystemInfo(&info);
        page_bytes_ = info.dwPageSize;
#else
        const long page = sysconf(_SC_PAGESIZE);
        if (page <= 0)
            throw std::runtime_error("sysconf page size failed");
        page_bytes_ = static_cast<std::size_t>(page);
#endif
        if (count == 0 ||
            count > (std::numeric_limits<std::size_t>::max() - 3 * page_bytes_) / sizeof(T))
            throw std::runtime_error("invalid guarded buffer size");
        const std::size_t bytes = count * sizeof(T);
        writable_bytes_ = ((bytes + page_bytes_ - 1) / page_bytes_) * page_bytes_;
        allocation_bytes_ = writable_bytes_ + 2 * page_bytes_;
#if defined(_WIN32)
        allocation_ = VirtualAlloc(nullptr, allocation_bytes_, MEM_RESERVE, PAGE_NOACCESS);
        if (!allocation_)
            throw std::runtime_error("VirtualAlloc reserve failed");
        auto *middle = static_cast<std::byte *>(allocation_) + page_bytes_;
        if (!VirtualAlloc(middle, writable_bytes_, MEM_COMMIT, PAGE_READWRITE)) {
            VirtualFree(allocation_, 0, MEM_RELEASE);
            throw std::runtime_error("VirtualAlloc commit failed");
        }
#else
        allocation_ =
            mmap(nullptr, allocation_bytes_, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (allocation_ == MAP_FAILED)
            throw std::runtime_error("mmap failed");
        auto *middle = static_cast<std::byte *>(allocation_) + page_bytes_;
        if (mprotect(middle, writable_bytes_, PROT_READ | PROT_WRITE) != 0) {
            munmap(allocation_, allocation_bytes_);
            throw std::runtime_error("mprotect failed");
        }
#endif
        data_ = reinterpret_cast<T *>(middle + (at_end ? writable_bytes_ - bytes : 0));
        for (std::size_t i = 0; i < count_; ++i)
            std::construct_at(data_ + i);
    }
    guarded_buffer(const guarded_buffer &) = delete;
    guarded_buffer &operator=(const guarded_buffer &) = delete;
    ~guarded_buffer() {
#if defined(_WIN32)
        VirtualFree(allocation_, 0, MEM_RELEASE);
#else
        munmap(allocation_, allocation_bytes_);
#endif
    }
    void read_only() {
        auto *middle = static_cast<std::byte *>(allocation_) + page_bytes_;
#if defined(_WIN32)
        DWORD previous{};
        if (!VirtualProtect(middle, writable_bytes_, PAGE_READONLY, &previous))
            throw std::runtime_error("VirtualProtect failed");
#else
        if (mprotect(middle, writable_bytes_, PROT_READ) != 0)
            throw std::runtime_error("mprotect read only failed");
#endif
    }
    T &operator[](std::size_t i) noexcept { return data_[i]; }
    const T &operator[](std::size_t i) const noexcept { return data_[i]; }
    std::span<T> span() noexcept { return {data_, count_}; }
    std::span<const T> const_span() const noexcept { return {data_, count_}; }
};
} // namespace safety
