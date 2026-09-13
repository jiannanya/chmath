#!/usr/bin/env sh
set -eu
project_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
compiler=${CXX:-c++}
jobs=${CHMATH_JOBS:-8}
for config in Debug Release Scalar; do
    build_type=$config
    simd=ON
    if [ "$config" = Scalar ]; then build_type=Release; simd=OFF; fi
    build_dir="$project_root/build/validate-$(basename "$compiler")-$config"
    cmake -S "$project_root" -B "$build_dir" -G Ninja "-DCMAKE_CXX_COMPILER=$compiler" "-DCMAKE_BUILD_TYPE=$build_type" "-DCHMATH_ENABLE_SIMD=$simd" -DCHMATH_BUILD_TESTS=ON -DCHMATH_BUILD_EXAMPLES=ON -DCHMATH_BUILD_BENCHMARKS=ON
    cmake --build "$build_dir" --parallel "$jobs"
    ctest --test-dir "$build_dir" --output-on-failure --parallel "$jobs"
done
