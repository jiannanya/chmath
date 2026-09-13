#pragma once
#include <array>
#include <chmath/core/scalar.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <random>
#include <string_view>

namespace test {
struct failure {};
struct entry {
    const char *name;
    void (*run)();
};
inline std::array<entry, 2048> registry{};
inline std::size_t count{}, assertions{};
struct registration {
    registration(const char *name, void (*run)()) {
        if (count == registry.size()) {
            std::fputs("Test registry capacity exceeded\n", stderr);
            std::abort();
        }
        for (std::size_t i = 0; i < count; ++i) {
            if (std::string_view(registry[i].name) == name) {
                std::fprintf(stderr, "Duplicate test name: %s\n", name);
                std::abort();
            }
        }
        registry[count++] = {name, run};
    }
};
inline void check(bool ok, const char *expr, const char *file, int line) {
    ++assertions;
    if (!ok) {
        std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", file, line, expr);
        throw failure{};
    }
}
inline std::mt19937_64 random{0x43484d415448ULL};
template <class T> T sample(T low = T(-10), T high = T(10)) {
    return std::uniform_real_distribution<T>(low, high)(random);
}
template <class A, class B, class T> bool near(const A &a, const B &b, T tolerance) {
    using chm::almost_equal;
    return almost_equal(a, b, tolerance, tolerance);
}
} // namespace test
#define CH_TEST(name)                                                                              \
    static void name();                                                                            \
    static test::registration reg_##name(#name, name);                                             \
    static void name()
#define CHECK(...) test::check(static_cast<bool>((__VA_ARGS__)), #__VA_ARGS__, __FILE__, __LINE__)
#define NEAR(...) CHECK(test::near(__VA_ARGS__))

int main(int argc, char **argv) {
    std::string_view selected;
    if (argc == 2 && std::string_view(argv[1]) == "--list") {
        for (std::size_t i = 0; i < test::count; ++i)
            std::puts(test::registry[i].name);
        return 0;
    }
    if (argc == 3 && std::string_view(argv[1]) == "--case" && argv[2][0] != '\0')
        selected = argv[2];
    else if (argc != 1) {
        std::fputs("usage: test_executable [--list | --case NAME]\n", stderr);
        return 2;
    }
    std::size_t failed{}, executed{};
    for (std::size_t i = 0; i < test::count; ++i) {
        if (!selected.empty() && selected != test::registry[i].name)
            continue;
        ++executed;
        // Each case is reproducible both alone and within the unified suite.
        std::uint64_t seed = 14695981039346656037ULL;
        for (const char c : std::string_view(test::registry[i].name)) {
            seed = (seed ^ static_cast<unsigned char>(c)) * 1099511628211ULL;
        }
        test::random.seed(seed);
        try {
            test::registry[i].run();
            std::printf("PASS %s\n", test::registry[i].name);
        } catch (const test::failure &) {
            ++failed;
        } catch (const std::exception &e) {
            ++failed;
            std::fprintf(stderr, "%s: %s\n", test::registry[i].name, e.what());
        } catch (...) {
            ++failed;
            std::fprintf(stderr, "%s: unknown exception\n", test::registry[i].name);
        }
    }
    std::printf("%zu cases, %zu checks, %zu failed\n", executed, test::assertions, failed);
    return executed == 0 ? 2 : (failed == 0 ? 0 : 1);
}
