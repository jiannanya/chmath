# chmath 2.0 验证与性能报告

验证日期：2026-09-12。目录：`E:/cc/AI/tokmon/chmath`。C++ 命名空间为 **`chm`**，CMake 导出目标为 **`chm::chmath`**，包名及头文件路径继续使用 `chmath`。不提供旧命名空间别名；安装消费程序使用 `find_package(chmath 2 CONFIG REQUIRED)`。

> 2026-09-21/22 的后续优化（并行批处理、热点收敛与临时对象削减）见文末 [2026-09-22 优化续作](#2026-09-22-优化续作)。该轮把套件扩展到 1019 条独立用例，并新增独立并行基准；本文其余内容仍为 2026-09-12 的验证记录。

测试由 300 条增至 **1000 条独立命名用例**。统一执行完成 **1,250,926 次运行时检查，0 失败**。[1000 条清单及源码位置](results/2026-09-12/test-inventory.json)、[统一运行日志](results/2026-09-12/unified-tests.log)。循环次数、压力迭代次数、编译期断言和跨编译器重复执行均不计入独立用例数量。

## 实际验证结果

| 平台 / 编译器 | 配置 | 结果 |
| --- | --- | --- |
| Windows x64 / MSVC 19.38.33130 | Debug、Release、Release + SIMD OFF | 每种 **1004/1004** 通过 |
| Windows x64 / Clang 17.0.6 / MSVC 标准库 | Debug、Release、Release + SIMD OFF | 每种 **1004/1004** 通过 |
| Linux x86_64 / WSL2 / GCC 13.3.0 | Debug、Release、Release + SIMD OFF | 每种 **1004/1004** 通过 |
| Linux x86_64 / Clang 18.1.3 | Debug + ASan + UBSan + 泄漏检测 + coverage | **1005/1005** 通过 |
| Linux x86_64 / Clang 18.1.3 | 独立 ThreadSanitizer 压力程序 | **60/60**，827,652 次检查，未报告数据竞争 |
| Windows x64 / MSVC、Clang | Release + 固定核心 + 相对性能阈值 ON | 各 **42/42** 性能用例通过 |
| Windows x64 / Clang 17 | `-fno-exceptions -fno-rtti` 场景示例 | 编译、执行通过 |

1004 项由 1000 条用例及堆分配、ODR、安装消费程序、场景示例四项组成；覆盖率配置另加统一程序。全部 **18 个头文件**（含内部实现头）独立编译。项目验证目标使用 C++20、严格告警且告警视为错误。新多右端求解另有 float/double 次正规矩阵的 C++20 编译期断言。安装脚本也验证了未指定构建类型的单配置工程。

本机实际运行范围是 Windows/Linux x86_64、SSE2 和显式 SIMD 关闭模式。macOS、ARM/NEON 未在本机运行；提供 CI 与回退实现，不把代码路径存在当作硬件验证完成。

- [MSVC 三种配置](results/2026-09-12/validation-windows-msvc.log)、[Windows Clang](results/2026-09-12/validation-windows-clang.log)、[Linux GCC](results/2026-09-12/validation-linux-gcc.log)。
- [ASan/UBSan](results/2026-09-12/sanitizers.log)、[ThreadSanitizer](results/2026-09-12/thread-sanitizer.log)、[空构建类型安装](results/2026-09-12/install-empty-config.log)、[无异常/RTTI 示例](results/2026-09-12/no-exceptions.log)。
- [编译后清单核对](results/2026-09-12/inventory-audit.txt)、[验证源码 SHA256 清单](results/2026-09-12/source-manifest.json)。

## 测试组成

| 分组 | 独立用例 | 验证内容 |
| --- | ---: | --- |
| 原有完整套件 | 300 | 标量、向量、矩阵、旋转、变换、几何、曲线、数值、35 条安全与 22 条性能用例 |
| 新向量功能 | 120 | float/double × 6 种维度 × 10 种代数、独立范数参考、极端尺度与紧凑布局契约 |
| 新矩阵功能 | 120 | float/double × 5 种尺寸 × 12 种乘法、求解、逆、LU/QR/Cholesky/特征分解契约 |
| 新旋转功能 | 96 | 两种精度、四种轴、四种角度；Rodrigues 独立参考、方向映射、插值弧 |
| 新曲线功能 | 48 | 两种精度与 1/2/3/8 次曲线；Bernstein 参考、联合导数、分割、样条及有理曲线 |
| 新数值功能 | 48 | 已知根、多项式参考、非有理根、解析积分、补偿累加、调用次数预算 |
| 新批处理内存安全 | 144 | 两种精度 × 六种操作 × 12 种长度，保护页、只读输入、原地与通道置换 |
| 新边界与接口契约 | 44 | 非有限值、次正规数、失败时输出不变、AoS、缓存所有权、点/球/盒裁剪、投影约定 |
| 新压力测试 | 60 | 两种精度 × 10 种负载 × 3 级工作量，含 2/4/8 线程共享只读分解和曲线 |
| 新性能测试 | 20 | 两种精度 × 10 种批处理、共享求解、曲线与缓存查询负载 |
| **合计** | **1000** | 其中 **179 条安全、42 条性能、60 条压力** |

新增用例通过共享断言函数和显式登记文件组织，参数分别表示真实类型、尺寸、角度、批量长度或工作量；不会只更换随机种子来凑数。`scripts/generate_extended_tests.py --check` 检查登记表与场景定义一致，正常 C++ 构建不需要 Python。CMake、运行器和清单审计均检查重复名称，完整套件少于 1000 条时配置失败。Release 中检查仍启用。

批处理安全测试用 Windows `VirtualAlloc/VirtualProtect` 或 POSIX `mmap/mprotect` 设置真实不可访问保护页，缓冲区分别紧贴头尾边界。新增长度为 0/1/2/3/4/5/7/8/9/15/16/17；非空场景验证分离输出、完整原地、通道循环置换。span 仍须指向存活且合法的内存；不把悬空指针、虚构容量或非法下标当作受支持输入。

压力工作量为 512、4096、32768 单位，覆盖归一化、旋转/仿射累积、LU 复用、射线网格、8 次样条、批处理流水线、条件数变化的矩阵与振荡积分。积分每 128 单位执行一次完整查询，最多 8193 次函数求值。线程使用局部计数和互不重叠的结果数组，断言在 join 后执行。

## 优化和新增接口

- 向量范数与归一化对普通数值直接计算，极大/极小及非有限输入使用原有稳定回退；不采用近似倒平方根。
- 小矩阵乘法使用局部累加；LU 合并验证/行缩放扫描，共享行尺度和主元倒数。新增矩阵右端 `solve(A,B)` / `lu.solve(B)`，求逆共享多右端回代。常量求值使用直接除法，避免缓存倒数在次正规尺度溢出。
- 四元数插值复用已归一化输入，消除最终归一化可抵消的公共系数；旋转矩阵检查以三重积判断方向；向量间旋转复用轴长度和半角。
- 新增 AoS/SoA 方向变换、批量旋转/归一化、SoA 叉积。float/double 共用 SSE2/NEON 包装，保持非对齐加载、尾部和重叠约定。double 点积保留可自动向量化的连续循环：实测强制双通道展开在 Clang 下反而较慢。
- 新增 Bézier 联合值/导数、三次 `cubic_polynomial` Horner 复用、拥有射线副本的 `prepared_ray`、批量 AABB 查询和视锥分类。缓存接口保留明确的浮点范围和数据生命周期约束，见 [API](API.md)。

公共布局保持 `vec3f=12`、`quatf=16`、`mat3f=36`、`affine3f=48`、`mat4f=64` 字节，double 版本为两倍。核心无可变全局缓存、虚表或隐式堆分配。全局普通及对齐 `new/new[]` 计数覆盖旧功能和全部新增功能，结果 **0 次堆分配**；用户缓冲区、测试框架和回调自行分配不在保证范围内。[分配日志](results/2026-09-12/allocation.log)

## 覆盖率

Clang 18 统一测试程序使用全新 profile 目录，仅统计 `src/chmath`。

| 指标 | 命中 / 总数 | 比例 |
| --- | ---: | ---: |
| 行 | 2114 / 2175 | 97.20% |
| 区域 | 2012 / 2084 | 96.55% |
| 分支 | 1183 / 1392 | 84.99% |
| 已生成函数组 | 280 / 280 | 100.00% |

[原始覆盖率](results/2026-09-12/coverage.txt)。本地逐行报告：`build/linux-sanitize/coverage/html/index.html`；同目录包含完整 JSON 和未覆盖区域索引。新增代码改变了统计分母；覆盖率不证明所有模板实例、所有输入、ARM 分支或近退化几何均正确。编译期分支由静态断言验证，无法全部体现为运行时行命中。

## 性能测量

Intel Core i9-12900K / Windows x64，固定逻辑处理器 0。所有编译/压力任务完成后顺序测量，CMake Release；不启用 fast-math、LTO 或本机专用指令开关。修改前二进制在改动前保存，日志包含 before/after 的 SHA256。

每进程 3 次预热、9 组计时、每组 12 次完整工作；每版本 3 个独立进程，表格取三个进程中位数的中位数。扩展基准点数 65536、矩阵数 2048；原核心基准矩阵数 4096。预分配缓冲区、不可内联调用边界及被消费的输出防止移除工作。

| 运算 | MSVC ns/元素：前 → 后 | 速度比 | Clang ns/元素：前 → 后 | 速度比 |
| --- | ---: | ---: | ---: | ---: |
| mat3f multiply | 7.926 → 7.064 | 1.12× | 3.654 → 3.642 | 1.00× |
| mat3d inverse | 62.211 → 51.168 | 1.22× | 68.335 → 46.847 | 1.46× |
| mat4d inverse | 105.436 → 98.071 | 1.08× | 106.421 → 68.526 | 1.55× |
| normalize vec3d | 8.499 → 6.854 | 1.24× | 4.044 → 2.378 | 1.70× |
| quaternion nlerp | 88.944 → 46.818 | 1.90× | 26.119 → 15.324 | 1.70× |
| quaternion slerp near | 148.478 → 55.098 | 2.69× | 47.465 → 15.177 | 3.13× |
| quaternion slerp wide | 125.627 → 75.346 | 1.67× | 56.075 → 30.103 | 1.86× |
| quaternion from matrix | 117.765 → 53.459 | 2.20× | 99.955 → 60.693 | 1.65× |
| rotation between | 84.029 → 65.938 | 1.27× | 60.978 → 41.606 | 1.47× |
| SoA double transform | 1.601 → 1.146 | 1.40× | 1.191 → 1.160 | 1.03× |
| SoA double dot | 0.634 → 0.570 | 1.11× | 0.568 → 0.585 | 0.97× |

完整原始记录和所有负载（包括未提速项目）：[MSVC 扩展比较](results/2026-09-12/extended-msvc-comparison.json)、[Clang 扩展比较](results/2026-09-12/extended-clang-comparison.json)、[MSVC 原核心比较](results/2026-09-12/core-msvc-comparison.json)、[Clang 原核心比较](results/2026-09-12/core-clang-comparison.json)。对应 `*-before.txt` / `*-after.txt` 保留全部样本和二进制标识。

探索性测量发现 double 点积的手写 SIMD 回归，已据对照实验改为编译器自动向量化；初测也遇到未改动负载的大幅计时波动，额外复测用于区分系统噪声与实现退化。探索记录保留在 `results/2026-09-12/exploratory/`。最终表格使用最终源码重新编译后的统一测量，不据单次最快值宣传性能。固定核心不能锁定频率或消除中断，小差异需谨慎解读；不宣称所有运算、编译器和硬件都提速。

42 条性能用例均检查计算结果、工作预算和计时有效性；其中 30 条与对应参考实现比较，在稳定 runner 上启用 2–6 倍的严重退化上限。阈值不是速度承诺，普通 CI 默认关闭耗时阈值。最终固定核心阈值日志：[MSVC 原有 22 条](results/2026-09-12/performance-msvc-test_performance.txt)、[MSVC 新增 20 条](results/2026-09-12/performance-msvc-test_extended_performance.txt)、[Clang 原有 22 条](results/2026-09-12/performance-clang-test_performance.txt)、[Clang 新增 20 条](results/2026-09-12/performance-clang-test_extended_performance.txt)。

## 复现

```powershell
./scripts/validate.ps1 -Compiler cl
./scripts/validate.ps1 -Compiler clang++
python scripts/test_inventory.py --build-dir build/validate-cl-Release
python scripts/generate_extended_tests.py --check
./scripts/benchmark.ps1 -Executable build/validate-cl-Release/chmath_bench_extended.exe -Processor 0 -Repeats 3
```

```sh
CXX=g++ sh scripts/validate.sh
cmake -S . -B build/sanitize-current -G Ninja -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_BUILD_TYPE=Debug -DCHMATH_SANITIZERS=ON -DCHMATH_COVERAGE=ON
cmake --build build/sanitize-current --parallel 4
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 LLVM_PROFILE_FILE=/dev/null ctest --test-dir build/sanitize-current --output-on-failure --parallel 4
mkdir -p build/sanitize-current/profiles-fresh
LLVM_PROFILE_FILE=build/sanitize-current/profiles-fresh/%p-%m.profraw build/sanitize-current/chmath_coverage
python3 scripts/coverage.py build/sanitize-current --profiles build/sanitize-current/profiles-fresh --llvm-profdata llvm-profdata-18 --llvm-cov llvm-cov-18
clang++-18 -std=c++20 -O1 -g -fsanitize=thread -fno-omit-frame-pointer -pthread -Isrc tests/test_extended_stress.cpp -o build/stress-tsan
TSAN_OPTIONS=halt_on_error=1 build/stress-tsan
```

请为不同编译器使用不同构建目录；更换编译器会使 CMake 重置缓存选项。覆盖率每次使用新的 profile 目录。相对性能阈值通过 `-DCHMATH_PERFORMANCE_CHECKS=ON` 启用，仅在稳定 Release 环境运行两个性能可执行文件。

前一版本的 300 条测试报告保留于 [2026-09-08 报告](VALIDATION-2026-09-08.md)，其数字不作为本轮结果。

## 2026-09-22 优化续作

本轮在 2026-09-12 报告之后继续推进全库性能与内存占用优化，分为两部分：**并行批处理入口**（新增能力，默认关闭）与**热点函数收敛**（更短的依赖链、更少的临时对象和超越函数调用）。公共值类型的尺寸、对齐、ABI、错误处理语义与核心零分配保证均未改变；并行入口是显式传入 `workers` 的独立重载，默认调用仍走串行路径。

### 环境与方法

- 机器：Apple Silicon（macOS 26，`hw.ncpu=16`，其中 12 个性能核心），Apple Clang 17.0.0，`-O3 -DNDEBUG`，NEON 后端。
- macOS 不支持 Windows 基准脚本的逻辑处理器绑定，因此数字取每进程 9 组样本的中位数、3 个进程再取中位数；前后源码使用相同命令、相同数据与固定随机种子。预分配缓冲区、不可内联调用边界和被消费的输出防止工作被优化掉。小于 5% 的差异不作为结论。
- 测试套件扩展到 **1019 条独立命名用例**（180 条安全、42 条性能、60 条压力）。`test_parallel` 始终以 `CHMATH_ENABLE_PARALLEL=1` 编译，在默认构建中也会真正创建线程并逐元素对比串行结果。

| 配置 | 全部测试 | 结果 |
| --- | ---: | --- |
| Release（SIMD ON，并行宏关闭） | 1023/1023 | 通过 |
| Release + `CHMATH_ENABLE_PARALLEL=ON` + `CHMATH_BUILD_BENCHMARKS=ON` | 1023/1023 | 通过 |

（1023 = 1019 条命名用例 + 堆分配计数、ODR、安装消费程序、场景示例；无覆盖率配置时不含统一覆盖率程序。）

### 热点收敛

单位为 ns/项；「前」为 2026-09-21 优化前的工作树，「后」为本轮结束。矩阵/四元数微基准为 4096 项 × 9 组样本；扩展基准为 65536 项。

| 运算 | 前 | 后 | 变化 |
| --- | ---: | ---: | ---: |
| `quaternion nlerp`（扩展基准） | 4.76 | 3.36–3.46 | 约 −28% |
| `quaternion slerp near`（扩展基准） | 5.32 | 3.72–3.85 | 约 −29% |
| `quaternion slerp wide`（扩展基准） | 16.0 | 12.4–13.1 | 约 −21% |
| `quaternion from matrix`（扩展基准） | 18.5 | 13.7–14.4 | 约 −24% |
| `chmath quaternion rotate`（核心基准） | 0.578 | 0.522 | 约 −10% |
| 4×4 float 矩阵 × 向量（`transform_point(mat4f)`） | 2.85 | 2.48 | 约 −13% |
| `mat3f * vec3f` | 0.97 | 0.88 | 约 −9% |
| `symmetric_eigen` 3×3 | 355 | 300 | 约 −15% |
| `cholesky` 3×3 | 9.8 | 8.1 | 约 −18% |

未观察到回退的负载（变化在噪声范围内）：`mat3f`/`mat4d` 矩阵乘法、Bezier 联合求值、ray AABB 查询、SoA double 变换与点积、`normalize vec3d`、`factor_qr`。

对应实现变化：

- 矩阵 × 向量、仿射复合、四元数乘法与旋转改为分量式局部累加，避免对零初始化结果元素的读改写和中间向量；4×4 float 矩阵 × 向量使用与矩阵乘法相同的 SSE2/NEON 内核，常量求值仍走通用实现。
- 四元数插值把平方范数落在 `epsilon` 内的输入视为单位，跳过近似恒等的重归一化；最终结果仍会归一化，明显缩放、零或非有限输入保持原有稳定路径。
- Cholesky 与对称特征分解逐元素校验对称性（矩阵乘法路径中 `transpose(A)*A` 的元素恰为列点积），不再生成转置矩阵；Jacobi 旋转用 `hypot` 形式的正切公式替代 `atan2` + `cos` + `sin`。
- `from_matrix` 直接在列上检查正交性与右手性，不再构造 `transpose(m)*m`；`decompose` 用三重积判断镜像，不再调用 LU 行列式；`compose` 就地缩放列。
- 这些改动同时减少了栈上临时对象：对称分解去掉两个 N×N 转置临时量，矩阵转换去掉一个 3×3 乘积，仿射复合去掉两个 3×3 临时量。`test_allocation` 仍报告核心接口 **0 次堆分配**。

### 并行批处理

`workers` 重载以 `CHMATH_ENABLE_PARALLEL=ON` 编译后启动 `std::thread` 工作线程；未启用、`workers <= 1`、批量小于 `parallel_min_chunk`（默认 32768）或平台拒绝创建线程时自动串行完成。分段按元素类型的整条缓存行对齐，整批契约在启动线程前校验，异常在所有线程汇合后重抛。1,048,576 项、16 个报告的硬件线程：

| 负载 | workers=1 | workers=16 | 加速 |
| --- | ---: | ---: | ---: |
| SoA double transform | 0.844 | 0.265 | 3.18× |
| SoA double dot | 0.346 | 0.159 | 2.17× |
| AoS double normalize | 1.063 | 0.247 | 4.30× |
| prepared ray AABB | 15.13 | 1.92 | 7.88× |
| frustum classify | 9.84 | 1.37 | 7.18× |

262,144 项时对应加速为 2.17× / 1.38× / 2.59× / 4.75× / 5.96×：变换与点积受内存带宽与线程创建成本限制，扩展性低于每项计算更重的查询。65,536 项、每项工作很少的点积可能慢于串行（默认阈值因此保持保守）；调用方对廉价内核应只在批量足够大时传 `workers`。每个并行调用自行创建线程，不使用全局线程池，多个互不相交的调用可以安全并发。

### 复现

```sh
cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCHMATH_BUILD_BENCHMARKS=ON
cmake --build build/release --parallel
ctest --test-dir build/release --output-on-failure
./build/release/chmath_bench 65536
./build/release/chmath_bench_extended 65536
./build/release/chmath_bench_parallel 1048576
cmake -S . -B build/parallel -G Ninja -DCMAKE_BUILD_TYPE=Release -DCHMATH_ENABLE_PARALLEL=ON
cmake --build build/parallel --parallel
ctest --test-dir build/parallel --output-on-failure
```

本机数字来自单台 ARM64 macOS，未固定核心或频率；它们不代表 x86、MSVC 或其他负载。CI 新增 `parallel` 作业在 Ubuntu 上用 `CHMATH_ENABLE_PARALLEL=ON` 构建并运行完整套件，并冒烟运行并行基准，但不启用耗时阈值。
