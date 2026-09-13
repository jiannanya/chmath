# 验证与性能报告

验证日期：2026-09-08。工作目录：`E:/cc/AI/tokmon/chmath`。

本轮从 33 条扩展到 **300 条独立命名、可单独运行的测试**，增加 267 条。统一运行完成 **312,312 次检查，0 失败**；不把循环次数、不同配置的重复运行或编译期断言计作独立用例。[完整输出](results/unified-tests.log)、[300 条清单（含源码位置）](results/test-inventory.json)。

## 实际构建结果

| 系统 / 编译器 | 配置 | 结果 |
| --- | --- | --- |
| Windows x64 / MSVC 19.38.33130 | Debug、Release、Release + SIMD OFF | 每种配置 **304/304 CTest 通过** |
| Windows x64 / Clang 17.0.6，MSVC 标准库 | Debug、Release、Release + SIMD OFF | 每种配置 **304/304 CTest 通过** |
| Linux x86_64 / WSL2 / GCC 13.3.0 | Debug、Release、Release + SIMD OFF | 每种配置 **304/304 CTest 通过** |
| Linux x86_64 / WSL2 / Clang 18.1.3 | Debug + ASan + UBSan + coverage | **305/305 CTest 通过** |
| Windows x64 / MSVC、Clang | Release + 相对性能阈值 ON，固定核心 | 各 **22/22 性能用例通过** |
| Windows x64 / Clang 17.0.6 | `-fno-exceptions -fno-rtti`，场景示例 | 编译、执行通过 |

304 项由 300 条独立用例、全局堆分配计数、多翻译单元 ODR、安装消费程序、场景示例组成；覆盖率配置另加统一执行程序。所有验证构建使用 C++20，项目目标启用严格告警并将告警作为错误。全部 14 个源头文件（含 1 个内部 SIMD 内核头）分别独立编译。安装测试实际执行安装，并通过 `find_package(chmath 1 CONFIG REQUIRED)` 构建、运行独立消费程序。

CMake 拒绝重复名称或少于 300 条的完整套件；运行器也检查注册表容量和重复名称。每条用例的随机种子由名称独立派生，单独与统一运行可复现。库存审计对比真实可执行文件的 `--list` 输出，已在 Windows/MSVC 与 Linux/GCC 上核对。[审计记录](results/inventory-audit.txt)。

macOS、ARM/NEON 未在本机执行。仓库提供相应 CI 路径，但本地结果不等于远程 CI 已运行。这里实际验证了 Windows/Linux x86_64、SSE2 和关闭显式 SIMD 的通用实现。

原始日志：

- [MSVC 三种配置](results/validation-windows-msvc.log)、[Clang/Windows 三种配置](results/validation-windows-clang.log)、[GCC/Linux 三种配置](results/validation-linux-gcc.log)。
- [MSVC Release 逐项输出](results/ctest-windows-msvc-release.log)、[ASan/UBSan CTest](results/ctest-linux-sanitizers.log)、[无异常/RTTI 示例](results/no-exceptions.log)。

## 测试组成

| 分组 | 独立用例 | 重点 |
| --- | ---: | --- |
| 原有核心/随机性质测试 | 33 | float/double、解析对照、矩阵残差、旋转、投影、几何、曲线、批处理 |
| 标量边界 | 30 | NaN、无穷、负零、容差、极端插值、乘除缩放 |
| 向量边界 | 30 | 零与次正规数、极端范数、维度、临界折射、非法投影 |
| 矩阵边界 | 35 | 多种尺寸、秩亏、主元、渐进/完全下溢、溢出、分解、复用求解 |
| 旋转与变换边界 | 40 | 退化轴、极端眼点、视口、左右手/深度约定、往返、非有限值 |
| 几何边界 | 35 | 退化盒/线段、切点、三角形边和顶点、极端尺度、OBB、视锥 |
| 曲线与数值边界 | 40 | 32 次 B-spline、重节点、非法权重、极端导数、根与积分预算 |
| 内存安全 | 35 | 真实保护页、只读输入、SIMD 尾部、重叠拒绝、原地与通道置换 |
| 性能 | 22 | 矩阵、归一化、旋转、AoS/SoA、曲线、分解、积分 |
| **合计** | **300** | **312,312 次运行时检查** |

安全测试在 Windows 使用 `VirtualAlloc/VirtualProtect`，在 POSIX 使用 `mmap/mprotect`。可读写区域两侧为不可访问页；逻辑缓冲区分别紧贴头部或尾部，越过最后一个合法 SIMD 分量就会触及保护页。覆盖 1/2/3/4/5/7/8/17 元素长度、完整原地、通道循环置换、只读输入、AoS/SoA；并检查所有输入/输出通道的部分重叠和长度错误，验证拒绝时输出不变。

独立分配测试替换普通及对齐的全局 `new/new[]`，循环执行矩阵、分解、旋转、求交、曲线、积分和批处理，新增 mat4 SIMD 与极端数值回退路径，结果为 **0 次堆分配**。[分配日志](results/allocation.log)。运行器、用户缓冲区和用户回调自行分配的内存不属于此统计。

## 本轮修复和优化

- `almost_equal` 拒绝非有限或负容差，避免极大异号数因 `inf <= inf` 而误判；`smoothstep` 支持端点差溢出的有限范围。
- 向量长度的 NaN 行为不再取决于分量顺序；归一化合并扫描并复用倒数；折射改善临界角计算并拒绝非法折射率；向量投影拒绝非有限输入/输出。
- 行列式在中间乘积溢出、完全下溢或渐进下溢时使用尾数/指数累积。LU 消元溢出时尝试缩放重算，避免把非零极大行列式误报为零。
- 视图方向差溢出时缩放计算；投影与反投影统一验证视口全部字段有限。
- 球体重叠处理中心差和半径和同时溢出；平面构造先归一化边，处理极大/极小三角形；线段最近点处理极端端点及极远查询点。
- B-spline 导数处理控制点差、节点宽度或倒数的中间溢出；Simpson 面板在极端范围使用指数缩放，恢复数学结果仍可表示的积分。
- mat4f 乘法新增 SSE2/NEON 内核和非对齐加载，常量求值保持通用 C++20 路径。`vec3f=12`、`quatf=16`、`mat4f=64`、`affine3f=48` 字节的布局保持不变。

这些修复不能消除一般浮点舍入误差、病态条件数或所有近退化几何问题；具体前提见 [API 文档](API.md)。

## 覆盖率

Clang 18 源码覆盖率，采用单个统一测试程序和新的 profile 目录，统计仅包含 `src/chmath`。

| 指标 | 本轮结果 |
| --- | ---: |
| 行覆盖 | **97.05%**，1876 / 1933 |
| 区域覆盖 | **96.73%**，1773 / 1833 |
| 分支覆盖 | **84.42%**，1057 / 1252 |
| 已生成函数组执行 | **244 / 244** |

上一轮行覆盖为 96.68%、分支覆盖为 82.83%；本轮新增数值保护和 SIMD 内核，统计分母也发生变化。[原始报告](results/coverage.txt)。本地逐行 HTML：`build/linux-sanitize/coverage/html/index.html`，同目录有 JSON 和未覆盖区域索引。统计仅描述实际实例化的模板和 x86 SSE2 路径，不代表所有类型/维度组合或 ARM 指令路径均已覆盖。

## 性能测量

Intel Core i9-12900K / Windows x64，单线程 CMake Release。`scripts/benchmark.ps1` 将每个进程固定到逻辑处理器 0；所有编译完成后顺序执行。未启用 fast-math、LTO 或本机专用指令开关。保留本轮修改前的编译产物作为 before，日志记录二进制 SHA256。

每进程 3 次预热、9 组计时、每组 12 次完整批处理；每版本启动 3 个独立进程，下表取三个进程中位数的中位数。点数 65536，矩阵数 4096。缓冲区预分配，不可内联工作函数及被消费的校验和防止删除计算；普通循环基线允许自动向量化。

| 运算 | 编译器 | 优化前 ns/元素 | 优化后 ns/元素 | 速度比 |
| --- | --- | ---: | ---: | ---: |
| mat4f 乘法 | MSVC | 11.951 | **3.361** | **3.56×** |
| mat4f 乘法 | Clang | 6.879 | **3.363** | **2.05×** |
| 稳定向量归一化 | MSVC | 6.480 | 6.320 | 1.03× |
| 稳定向量归一化 | Clang | 3.601 | **3.057** | **1.18×** |

矩阵乘法改善明显；Clang 归一化耗时下降约 15.1%，MSVC 归一化差异较小。其他未修改负载存在几个百分点的正负波动，不能据此声称所有运算均提速。固定核心不能锁定频率或消除中断；这是缓存友好的微基准，不包含分配、GPU 传输或完整应用流程。

- [MSVC before](results/benchmark-msvc-before.txt)、[MSVC after](results/benchmark-msvc-sse2.txt)、[MSVC 标量](results/benchmark-msvc-scalar.txt)、[MSVC 比较 JSON](results/comparison-msvc.json)。
- [Clang before](results/benchmark-clang-before.txt)、[Clang after](results/benchmark-clang-sse2.txt)、[Clang 比较 JSON](results/comparison-clang.json)。

22 条性能测试同时验证结果。10 种负载提供相对耗时阈值，覆盖四种矩阵尺寸、复用 LU、三种归一化范围、两种 SoA 大小；阈值为对应参考实现的 2–5 倍，用于检测严重退化。默认关闭时间阈值，正确性、预算和计时采集始终执行。固定核心的 Release 阈值运行在 MSVC 与 Clang 上均通过，各执行 22 条用例、206,872 次检查（含额外 10 次阈值检查）。[MSVC 阈值日志](results/performance-msvc-gated.log)、[Clang 阈值日志](results/performance-clang-gated.log)。

## 复现命令

Windows 开发者终端：

```powershell
./scripts/validate.ps1 -Compiler cl
./scripts/validate.ps1 -Compiler clang++
python scripts/test_inventory.py --build-dir build/validate-cl-Release
ctest --test-dir build/validate-cl-Release -R test_memory_safety --output-on-failure
./scripts/benchmark.ps1 -Executable build/validate-cl-Release/chmath_bench.exe -Processor 0 -Repeats 3
python scripts/benchmark_compare.py docs/results/benchmark-msvc-before.txt docs/results/benchmark-msvc-sse2.txt
```

稳定机器上单独启用性能阈值：

```powershell
cmake -S . -B build/perf -G Ninja -DCMAKE_BUILD_TYPE=Release -DCHMATH_PERFORMANCE_CHECKS=ON
cmake --build build/perf --target test_performance --parallel
./scripts/benchmark.ps1 -Executable build/perf/test_performance.exe -TestSuite -Repeats 1
```

Linux / macOS 普通验证：`CXX=g++ sh scripts/validate.sh` 或 `CXX=clang++ sh scripts/validate.sh`。Linux Clang 18 sanitizer/覆盖率示例，需要匹配的 Clang/LLVM 工具：

```sh
cmake -S . -B build/coverage -G Ninja \
  -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_BUILD_TYPE=Debug \
  -DCHMATH_COVERAGE=ON -DCHMATH_SANITIZERS=ON
cmake --build build/coverage --parallel
mkdir -p build/coverage/ctest-profiles build/coverage/unified-profiles
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
LLVM_PROFILE_FILE="$PWD/build/coverage/ctest-profiles/%p-%m.profraw" \
  ctest --test-dir build/coverage --output-on-failure --parallel 8
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
LLVM_PROFILE_FILE="$PWD/build/coverage/unified-profiles/coverage.profraw" \
  ./build/coverage/chmath_coverage
python3 scripts/coverage.py build/coverage --profiles build/coverage/unified-profiles \
  --llvm-profdata llvm-profdata-18 --llvm-cov llvm-cov-18
```

修改代码后使用新的 profile 目录。只合并统一程序的 profile 可避免独立头文件测试中未使用内联实例对统计产生干扰。
