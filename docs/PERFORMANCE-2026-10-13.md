# 2026-10-13 性能优化与跨平台验证

本轮以 `da183c4` 为基线，保留 C++20、header-only、公共类型布局、现有接口、SSE2 / NEON / 标量后端、可选 std::thread 并行设计。没有引入依赖、全局线程池、近似倒平方根、fast-math 或本机专用编译选项。

## 实现

- **AoS 批量归一化**：复用 packet 后端，一次处理多个 vec3。通过每个对象的实际元素收集/写回，避免把 12 字节 vec3 当作 16 字节加载，也不跨对象做标量指针运算。普通值使用精确平方根；含零、非有限数、溢出或下溢平方和的分组回退到原有稳定归一化，尾部逐点处理。保留完整原地支持、部分重叠拒绝、校验失败前不写入。
- **并行归一化**：分块边界对齐到全局 packet 边界，保证相同向量不会因线程数量变化而切换所在分组。新增混合普通值/极值、不同 workers、原地操作的精确相等测试。SIMD 与逐向量标量计算本身可能有末位舍入差异。
- **最小二乘**：与 `factor_qr` 共用 Householder 消元，但直接将反射作用于右端项，不再构造完整 Q 后计算 Qᵀb。对于 R 行 C 列、R≥C，工作空间由 O(R²+R*C) 降到 O(R*C+R)，运算量由 O(R²*C) 降到 O(R*C²)。保留满列秩、容差和有限值检查；反射投影或更新中间溢出时使用缩放回退。该算法与显式 Q 的舍入不必逐位相同，通过解、残差、正交性及极值测试验证。
- **完整 QR**：R>4 时按连续列累积/更新 Q，保持各行累加顺序，便于编译器对独立行向量化；仍返回原有完整 Q、R。
- **较大 LU**：N≥16 时先计算主元列乘子，再按列连续更新子矩阵，保留选主元、更新顺序和极小主元除法回退。8×8 的列遍历原型曾测得回退，因此小尺寸继续使用原来的行局部循环。
- **可复现测量**：新增 `chmath_bench_hotspots` 和跨平台 Python 交错运行脚本 `scripts/benchmark_ab.py`；Linux/macOS/Windows CI 增加新基准的冒烟执行。

## 原生 ARM 性能

环境为 Apple M4 Max、macOS 26.1、Apple Clang 17.0.0，Release `-O3 -DNDEBUG`，显式 SIMD 开启。基线与优化版使用同一份新增基准源码和相同编译参数。每个进程预热 3 次，再测 9 组、每组 12 轮；前后进程交替顺序运行 7 对，表中为进程中位数的中位数。批量长度 65,536，矩阵批量 256；单位为 ns/项，倍数为前/后耗时。

| 工作负载 | 优化前 | 优化后 | 加速比 |
| --- | ---: | ---: | ---: |
| AoS float transform | 0.344 | 0.344 | 1.00× |
| AoS float normalize | 1.300 | 0.622 | 2.09× |
| AoS float normalize mixed | 1.302 | 0.784 | 1.66× |
| AoS double transform | 0.690 | 0.695 | 0.99× |
| AoS double normalize | 1.054 | 0.810 | 1.30× |
| AoS double normalize mixed | 1.055 | 0.923 | 1.14× |
| least squares double 4x2 | 90.169 | 64.494 | 1.40× |
| least squares float 8x3 | 239.963 | 150.051 | 1.60× |
| least squares double 32x4 | 3916.396 | 686.347 | 5.71× |
| least squares double 64x8 | 26836.575 | 3185.588 | 8.42× |
| QR double 8x3 | 258.816 | 222.900 | 1.16× |
| QR double 32x4 | 3599.515 | 2652.194 | 1.36× |
| LU double 8x8 | 123.386 | 126.410 | 0.98× |
| LU double 16x16 | 759.806 | 564.562 | 1.35× |
| LU double 32x32 | 4688.327 | 2664.388 | 1.76× |

原有核心基准 8 项、扩展基准 15 项的前后中位数变化均在约 ±2% 内。覆盖小矩阵乘法/求逆、逐点归一化、四元数、SoA 变换/点积、Bezier 与光线包围盒查询；这些结果不计为本轮加速收益。

关闭显式 SIMD 后，AoS 归一化保持原标量路径，变化约 1% 以内；32×4 / 64×8 最小二乘仍约快 5.66× / 8.30×，32×32 LU 约快 1.70×。数值分解的收益来自算法和访问顺序，不依赖显式 SIMD。

x86_64/SSE2 二进制也完成交错测量；它们在 Rosetta 上执行，数据只作为第二后端的补充，不能代替原生 Intel/AMD 硬件结论。该环境中 AoS float/double 归一化约快 1.40× / 2.24×，64×8 最小二乘约快 8.73×；16×16 LU 约慢 2.2%，32×32 LU 约快 1.36×，并非每个尺寸都在所有平台加速。

## 并行测量的限制

1,048,576 项、16 workers 的 AoS double 归一化，两轮各 7 对交错测试分别为 0.267→0.249 ns/项、0.212→0.197 ns/项，约快 1.07×。单 worker 分别约快 1.35×。端到端并行收益受到线程创建、调度和内存带宽限制，不能直接套用内核加速倍数。

未改动的 frustum 16-worker 工作负载首轮出现 1.408→1.869 ns/项，中间样本在 1.3–5.0 ns/项间波动；复测为 1.319→1.318 ns/项。SoA transform 的 16-worker 行两轮也有相反方向的小幅变化。因此多线程数据作为带噪声的观察记录，不据此声称稳定收益或全场景无回退；两轮原始样本均保留。

所有测量的 7 次前后样本、中位数和比率见 [机器可读数据](PERFORMANCE-2026-10-13.json)。这些是本机观测，不是跨设备性能保证。

## 正确性和平台验证

| 本机配置 | CTest | 结果 |
| --- | ---: | --- |
| ARM64 Release / NEON | 1036 | 全部通过 |
| ARM64 Debug | 1036 | 全部通过 |
| ARM64 Release / SIMD OFF | 1036 | 全部通过 |
| ARM64 Release / PARALLEL ON | 1036 | 全部通过 |
| ARM64 Debug / ASan + UBSan | 1036 | 全部通过 |
| x86_64 Release / SSE2，Rosetta 执行 | 1036 | 全部通过 |

1036 项包含 1032 个独立命名用例，以及堆分配、ODR、安装消费和示例检查。新增 13 个命名用例覆盖保护页首尾、只读输入、非 SIMD 对齐长度、尾部、零/NaN/Inf/极大/极小值、原地和重叠拒绝，矩形/方形/单元素最小二乘、缩放和秩亏、QR 重建与正交性、较大 LU 选主元/重建/求逆/常量求值/极小主元，以及串并行精确一致性。

独立头文件检查通过；核心操作堆分配计数仍为 0，并已将最小二乘加入该计数路径。测试清单与生成用例同步检查通过。项目目标均启用严格警告与 warnings-as-errors。

本机没有执行 Linux GCC 或 Windows MSVC，也没有原生 x86 硬件；现有 CI 保留这些编译器的构建与测试，新用例自动进入 CI。x86 配置的主库测试为 x86_64，header-only 安装消费子工程使用宿主默认 ARM64 编译器设置。ARM32 等未在本机执行的目标继续使用现有编译期选择和标量回退。

## 复现

先保存旧版构建及头文件，再对相同源码、编译器和配置构建新版本。新增热点基准也必须用同一份 `benchmarks/hotspots.cpp` 分别针对旧/新头文件编译。本次旧版头文件快照和二进制保存在忽略的 `build/performance-20261013/before/`，完整测试及测量日志保存在 `build/performance-20261013/`。

```sh
cmake --preset release
cmake --build --preset release --parallel
ctest --preset release --output-on-failure
python3 scripts/benchmark_ab.py build/performance-20261013/before/chmath_bench_hotspots build/release/chmath_bench_hotspots --repeats 7 --output-dir build/benchmark-ab
python3 scripts/test_inventory.py --build-dir build/release
python3 scripts/generate_extended_tests.py --check
```

Windows 使用对应 `.exe` 路径及 `python`，交错脚本无需修改。其他本机配置使用 `debug`、`scalar`、`parallel`、`sanitize` 预设；x86 验证使用独立构建目录并设置 `-DCMAKE_OSX_ARCHITECTURES=x86_64`。
