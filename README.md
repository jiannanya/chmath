# chmath

**English** | [简体中文](README.zh-CN.md)

A C++20 mathematics library for 3D rendering, game engines, and CAD applications. The source is organized by module under `src/chmath/`. The library is header-only and has no third-party runtime dependencies.

**Version 2.0 uses the `chm` namespace** and exports the CMake target `chm::chmath`. The package name and include paths remain `chmath`. To migrate, replace `chmath::` with `chm::` in C++ code, update the CMake target, and use `find_package(chmath 2 CONFIG REQUIRED)`. No legacy namespace alias is provided.

The core uses compact, fixed-size value types, stack-allocated temporary storage, and explicit borrowing through `std::span`. Common aliases are provided for `float` and `double`. Vectors and matrices also support other arithmetic types; operations involving lengths, rotations, and decompositions require floating-point types. Batch operations for float/double and 4×4 float matrix multiplication can use SSE2 or NEON, with scalar fallbacks.

The test suite contains **1000 independently named, individually runnable cases**, including **179 memory safety, 42 performance, and 60 stress cases**, plus heap allocation, installation, and linking checks. See the [validation report (Chinese)](docs/VALIDATION.md) for the test inventory and measured results.

## Modules

| Directory | Features |
| --- | --- |
| `core/` | Constants, angle conversion, approximate comparison, interpolation, finite-value checks, and multiply/divide operations that avoid intermediate overflow |
| `vector/` | Generic N-dimensional vectors, dot and cross products, stable normalization, distance, projection, reflection, and refraction |
| `matrix/` | Generic R×C matrices, multiplication, transpose, outer products, determinants, LU, inversion, and linear solves with one or multiple right-hand sides |
| `rotation/` | Quaternions, axis-angle representation, fixed-axis XYZ Euler angles, matrix conversion, rotation interpolation, and rotation between vectors |
| `transform/` | 3×4 affine transforms, TRS composition/decomposition, normal transforms, look-at, perspective/orthographic projection, and screen projection |
| `geometry/` | Rays, segments, AABBs, spheres, planes, triangles, and OBBs; closest points, intersections, SAT, frustum culling, prepared rays, and batch queries |
| `curves/` | Bézier curves of any fixed degree, combined value/derivative evaluation, reusable cubic Horner evaluation, splitting, patches and partial derivatives, Hermite, Catmull–Rom, B-splines, and NURBS |
| `numeric/` | Stable quadratic equations, Horner evaluation, reusable LU solves, Householder QR, least squares, Cholesky, symmetric eigendecomposition, bisection, adaptive integration, and compensated summation |
| `simd/` | AoS/SoA point and direction transforms, rotation, normalization, SoA dot/cross products, and 4×4 float multiplication; unaligned loads, tail handling, in-place operations, and overlap checks |

## Quick start

```cpp
#include <chmath/chmath.hpp>

int main() {
    using namespace chm;
    const auto q = from_axis_angle(vec3f{0, 1, 0}, radians(45.0f));
    if (!q) return 1;

    const auto model = compose(trs<float>{{1, 2, -5}, *q, {2, 2, 2}});
    const vec3f world = transform_point(model, vec3f{1, 0, 0});

    const auto inv = inverse(model);
    if (!inv) return 2;
    const vec3f local = transform_point(*inv, world);
    return almost_equal(local, vec3f{1, 0, 0}) ? 0 : 3;
}
```

You can include individual module headers, such as `chmath/vector/vector.hpp`, as needed. Every public header has an independent compilation check.

To integrate the library as a subproject:

```cmake
add_subdirectory(path/to/chmath)
target_link_libraries(your_target PRIVATE chm::chmath)
```

Or install it first:

```sh
cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --parallel
ctest --test-dir build/release --output-on-failure
cmake --install build/release --prefix /your/install/prefix
```

```cmake
find_package(chmath 2 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE chm::chmath)
```

Without CMake, add the project's `src` directory to your include path and enable C++20. See [examples/scene.cpp](examples/scene.cpp) for a complete scene example.

## Conventions

- **Column-major storage and column vectors**: use `m(row, column)`; the storage index is `column * rows + row`. In `A * B * p`, B is applied first.
- All angles are in **radians**. Quaternions are stored in **x, y, z, w** order and use the Hamilton product.
- Quaternions and affine transforms default to the identity; vectors and general matrices default to zero.
- Access vector components through `v.x()` / `v.y()` / `v.z()` / `v.w()` or `v[i]`. The implementation does not use anonymous unions or pointer arithmetic across separate members.
- Projection defaults to a **right-handed coordinate system, Z ∈ [0,1], and forward depth**. Handedness, `[-1,1]` depth, and reverse Z are explicit options. Cameras look along -Z in right-handed view space and +Z in left-handed view space.
- `from_euler_xyz` applies rotations around the fixed X, Y, and Z axes in that order: `qz * qy * qx`.
- Viewport Y increases upward in `project` / `unproject`. Callers using a top-left window origin must flip Y. `unproject` takes a **precomputed inverse MVP matrix**.
- A ray is `origin + t * direction`; its direction need not be normalized. The parameter t represents distance only when the direction has unit length. Tangency and boundary contact count as hits.
- A plane's normal must have unit length; prefer `make_plane` / `plane_from_points` for construction. `rotate`, `to_matrix(quat)`, and `compose` also require unit quaternions.

## Memory and performance

Static assertions verify these layouts on conventional IEEE 754 platforms:

| Type | Bytes | Purpose |
| --- | ---: | --- |
| `vec3f` | 12 | Compact positions, normals, and directions |
| `vec4f` / `quatf` | 16 | Homogeneous coordinates / rotations |
| `mat3f` | 36 | Linear transforms |
| `affine3f` | 48 | Affine transforms with the fixed last row omitted |
| `mat4f` | 64 | Projection and general homogeneous transforms |

`vec3f` is aligned to `alignof(float)`. These are standard-layout, trivially copyable value types without pointers, virtual tables, or implicit extra components. Double-precision versions occupy twice the space. Compact CPU layouts do not automatically satisfy GPU uniform-buffer alignment requirements; pack uploads according to the graphics API's layout rules.

The implementation uses compile-time dimensions, contiguous column access, small inlineable functions, and SoA / mat4f SIMD kernels, without reference counting or runtime dispatch. Mat4f multiplication uses the generic implementation during constant evaluation and unaligned SIMD loads at runtime, without increasing object size or alignment. Normalization uses a squared-norm fast path for ordinary values and retains scaling protection for very large or small inputs. LU shares row-scale and pivot reciprocals, while multiple-right-hand-side substitution reduces repeated division. Small matrix products use local accumulators. Quaternion interpolation, matrix conversion, and rotation between vectors reuse intermediate results. Cubic polynomial and prepared-ray interfaces let repeated queries share preparation work.

The library does not enable `fast-math`, globally target host-specific CPU instructions, or use approximate reciprocal square roots by default, preserving NaN/infinity and error-checking semantics. Matrix inversion uses pivoted LU to balance generality and numerical stability. Reuse a `factor_lu` result when solving repeatedly with the same matrix.

SIMD paths preserve the ABI of public value types. Set `CHMATH_ENABLE_SIMD=OFF` in CMake to disable explicit SIMD; the compiler may still vectorize ordinary loops. Kernels are selected automatically on x86 with SSE2 and ARM toolchains supporting `__ARM_NEON`; other platforms use scalar implementations. Use consistent macro settings across translation units. Mixing different settings within one process is unsupported.

Batch operations use caller-provided memory and do not allocate output buffers. Exact in-place operations are supported; partial overlaps are rejected. SoA output channels must not overlap each other. Transform parameters themselves must not overlap output buffers. Invalid buffer sizes or overlaps return `false` before any output is written.

See the [validation report (Chinese)](docs/VALIDATION.md) for measurements, methodology, and complete logs. Performance depends on the CPU, compiler, data size, and cache behavior; a single benchmark does not predict every application.

## Error handling and numerical limits

Fallible construction, inversion, decomposition, and intersection functions primarily return `std::optional`, without allocating exception objects. `try_normalize` rejects zero vectors and nonfinite values. The convenience function `normalize(vec)` returns a zero vector on failure, while `normalize(quat)` returns the identity quaternion.

The default `epsilon<T>` is `32 * numeric_limits<T>::epsilon()`. `almost_equal` accepts finite, nonnegative relative and absolute tolerances and avoids incorrect comparisons caused by intermediate overflow. LU checks pivots relative to row scales, with a caller-adjustable tolerance. For CAD, prefer `double` and choose tolerances appropriate to the model's scale. Index checks use `assert`; valid indices remain a caller requirement in Release builds. Ordinary arithmetic follows C++ rules, so callers must avoid integer overflow, integer division by zero, and invalid indices.

B-splines and NURBS are **non-owning views**. Control points, knots, and weights must remain alive and unchanged while a view is in use. Construction validates the inputs; evaluation uses a stack-allocated de Boor workspace. The degree is fixed at compile time and limited to 0–32. NURBS accepts only strictly positive, finite weights. Evaluation outside the domain fails. Interior knots use the segment on the right; the domain's right endpoint uses the left-hand limit. Derivatives at discontinuous knots follow the same one-sided convention.

`bisect` requires a continuous function and a bracket with opposite signs. `integrate` uses adaptive Simpson integration with depth and function-evaluation budgets; its error is an estimate. Callers must check `converged`. Exceptions from callbacks propagate normally, and allocations made by callbacks are outside the library's zero-allocation guarantee.

This library provides general mathematical foundations for games, rendering, and CAD. It does not include solid Boolean operations, B-rep topology, exact geometric predicates, general mesh collision acceleration structures, sparse/dynamic large matrices, SVD, or a complete surface-modeling kernel. Floating-point intersections do not guarantee exact topology in every near-degenerate case or watertight triangle boundaries. `barycentric` exposes a relative degeneracy tolerance. Extreme condition numbers, unrepresentable results, or insufficient precision can still cause failure; applications requiring exact geometry should use an appropriate geometry kernel.

## Building, testing, and benchmarking

Requirements: a C++20 compiler and CMake 3.20 or later; presets require CMake 3.21. Presets and validation scripts use Ninja, while ordinary CMake builds can use other generators. On Windows, run the commands in a developer terminal with the MSVC and Windows SDK environment loaded.

```sh
cmake --preset release
cmake --build --preset release --parallel
ctest --preset release
./build/release/chmath_bench
```

On Windows, the benchmark executable is `build/release/chmath_bench.exe`. You can pass an item count, for example `chmath_bench 262144`. It reports the median and range of nine samples after warmup. Its ordinary SoA baseline is available for compiler optimization; automatic vectorization is not forcibly disabled.

Validate Debug, Release, and scalar fallback builds with:

```powershell
./scripts/validate.ps1 -Compiler cl
./scripts/validate.ps1 -Compiler clang++
```

```sh
CXX=g++ sh scripts/validate.sh
CXX=clang++ sh scripts/validate.sh
```

Each named case is registered independently with CTest. CMake rejects duplicate names or a full suite containing fewer than 1000 cases. Checks remain enabled in Release builds. Parameterized cases are named by actual scenarios such as type, dimension, angle, or batch length. Loop iterations and repeated execution under different compilers do not count as additional cases. Run individual categories or cases with:

```sh
ctest --test-dir build/release -L safety --output-on-failure
ctest --test-dir build/release -L stress --output-on-failure
ctest --test-dir build/release -L performance -V
./build/release/test_matrix_boundaries --case matrix_determinant_balanced_underflow_diagonal
python scripts/test_inventory.py --build-dir build/release --output build/test-inventory.json
python scripts/generate_extended_tests.py --check
```

The 179 memory safety cases check buffer contracts and use real guard pages through Windows `VirtualAlloc/VirtualProtect` or POSIX `mmap/mprotect`. They cover SIMD tails, read-only inputs, addresses without SIMD alignment, in-place operations, and rejection before writes. The 144 added cases cover float/double, six operations, and 12 lengths. Violations of documented caller preconditions, such as out-of-range indexing, are not exercised as valid inputs.

The 60 stress cases cover repeated normalization, accumulated transforms, decomposition reuse, ray grids, spline sampling, batch pipelines, ill-conditioned matrices, oscillatory integration, and shared read-only data with 2/4/8 threads. Each type/workload has three tiers: 512, 4096, and 32768 work units. Integration performs one complete evaluation-budgeted integral per 128 units.

The 42 performance cases include correctness comparisons, iteration/evaluation budgets, and timing after warmup. They use fixed random seeds, preallocated buffers, and consumed outputs. Timing is recorded by default. On a stable Release runner, `CHMATH_PERFORMANCE_CHECKS=ON` enables relative runtime ceilings for 30 comparable workloads. Ordinary CI with varying machine load does not enable timing thresholds. On Windows, pin benchmark processes to a logical processor for reproduction:

```powershell
./scripts/benchmark.ps1 -Executable build/release/chmath_bench.exe -Processor 0 -Repeats 3
./scripts/benchmark.ps1 -Executable build/release/chmath_bench_extended.exe -Processor 0 -Repeats 3
```

ASan + UBSan on Linux / macOS:

```sh
CXX=clang++ cmake --preset sanitize
cmake --build --preset sanitize --parallel
ctest --preset sanitize
```

| CMake option | Default | Description |
| --- | --- | --- |
| `CHMATH_BUILD_TESTS` | ON for top-level projects | CTest, standalone headers, multiple translation units, and installation consumer tests |
| `CHMATH_BUILD_EXAMPLES` | ON for top-level projects | Runnable 3D scene example |
| `CHMATH_BUILD_BENCHMARKS` | OFF | Benchmarks without external dependencies |
| `CHMATH_ENABLE_SIMD` | ON | Explicit SIMD batch operations and mat4f multiplication on supported platforms |
| `CHMATH_SANITIZERS` | OFF | ASan and UBSan with GCC / Clang on Unix-like systems |
| `CHMATH_COVERAGE` | OFF | Clang source-based coverage |
| `CHMATH_PERFORMANCE_CHECKS` | OFF | Relative performance thresholds on stable Release runners |
| `CHMATH_WARNINGS_AS_ERRORS` | ON | Strict warnings for this project's tests, examples, and benchmarks only |

Tests include deterministic randomized properties, analytic references, boundary and degenerate inputs, extreme numerical values, layout checks, SIMD/scalar differential checks, guard pages and sentinel buffers, and heap allocation counting. Source coverage does not prove correctness for every template instantiation or input; see the report for the actual validation scope. GitHub Actions configuration is in `.github/workflows/ci.yml` and covers Windows, Linux, macOS, and sanitizer builds.

See the [API reference (Chinese)](docs/API.md) for the detailed interface index and examples.
