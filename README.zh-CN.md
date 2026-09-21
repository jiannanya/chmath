# chmath

[English](README.md) | **简体中文**

面向 3D 渲染、游戏引擎和 CAD 应用的 C++20 数学库。源代码位于 `src/chmath/`，按功能模块组织，采用 header-only 形式，无第三方运行时依赖。

**2.0 使用 `chm` 命名空间**，CMake 导出目标为 `chm::chmath`。头文件路径和包名继续使用 `chmath`；迁移时将代码中的 `chmath::` 改为 `chm::`，CMake 中改用新目标和 `find_package(chmath 2 CONFIG REQUIRED)`。不提供旧命名空间别名。

核心采用固定大小的紧凑值类型、栈上临时存储和显式借用的 `std::span`。提供 `float` / `double` 常用别名；向量与矩阵也支持其他算术类型，涉及长度、旋转和分解的接口限浮点类型。float/double 批处理和 4×4 float 矩阵乘法可使用 SSE2 / NEON，并提供标量回退。

测试套件包含 **1019 个独立命名、可单独运行的用例**，其中有 **180 个内存安全、42 个性能、60 个压力用例**，另有堆分配计数、安装与链接检查。测试清单和本机实测结果见 [验证报告](docs/VALIDATION.md)。

## 模块

| 目录 | 功能 |
| --- | --- |
| `core/` | 常量、角度转换、近似比较、插值、有限值检测、避免中间溢出的乘除 |
| `vector/` | 通用 N 维向量、点积、叉积、稳定归一化、距离、投影、反射/折射 |
| `matrix/` | 通用 R×C 矩阵、乘法、转置、外积、行列式、LU、求逆、单/多右端线性方程求解 |
| `rotation/` | 四元数、轴角、固定轴 XYZ 欧拉角、矩阵转换、旋转插值、向量间旋转 |
| `transform/` | 3×4 仿射变换、TRS、分解、法线变换、look-at、透视/正交、屏幕投影 |
| `geometry/` | Ray、Segment、AABB、Sphere、Plane、Triangle、OBB；最近点、求交、SAT、视锥裁剪、缓存射线与批量查询 |
| `curves/` | 任意固定次数 Bézier、联合求值/导数、三次 Horner 缓存、曲线分割、曲面及偏导、Hermite、Catmull-Rom、B-spline、NURBS |
| `numeric/` | 稳定二次方程、Horner、多次求解可复用的 LU、Householder QR、最小二乘、Cholesky、对称特征分解、二分求根、自适应积分、补偿求和 |
| `simd/` | AoS / SoA 点与方向变换、旋转、归一化；SoA 点积/叉积，4×4 float 乘法；非对齐加载、尾部处理、原地计算、重叠检查 |
| `parallel/` | 为批处理入口的 workers 重载提供可选的 `std::thread` 拆分；串行回退、按缓存行对齐的分段、异常安全的线程汇合 |

## 快速使用

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

可按需只包含 `chmath/vector/vector.hpp` 等模块头。所有公共头文件均有独立编译检查。

通过子项目集成：

```cmake
add_subdirectory(path/to/chmath)
target_link_libraries(your_target PRIVATE chm::chmath)
```

或安装后使用：

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

不使用 CMake 时，把本项目 `src` 加入头文件搜索路径，并启用 C++20。完整场景例子见 [examples/scene.cpp](examples/scene.cpp)。

## 统一约定

- **列主序，列向量**：`m(row, column)`，内存索引为 `column * rows + row`；`A * B * p` 先执行 B。
- 所有角度使用**弧度**；四元数内存顺序为 **x、y、z、w**，乘法为 Hamilton 乘积。
- 四元数和仿射变换默认构造为单位变换；向量与普通矩阵默认构造为零。
- `vec` 通过 `v.x()` / `v.y()` / `v.z()` / `v.w()` 或 `v[i]` 访问；不使用匿名联合体或越界成员指针运算。
- 图形投影默认**右手坐标系、Z ∈ [0,1]、正向深度**；可显式切换左右手、`[-1,1]` 和反向 Z。相机右手空间看向 -Z，左手空间看向 +Z。
- `from_euler_xyz` 按固定 X、Y、Z 轴依次旋转，即 `qz * qy * qx`。
- `project` / `unproject` 的 viewport Y 向上；上方为原点的窗口坐标由调用方翻转 Y。`unproject` 接受**已求好的 MVP 逆矩阵**。
- 射线 `origin + t * direction` 的方向不要求单位长度；只有单位方向时 t 才等于距离。相切和边界接触计为命中。
- `plane` 的 normal 必须为单位向量；优先通过 `make_plane` / `plane_from_points` 构造。`rotate`、`to_matrix(quat)`、`compose` 也要求单位四元数。

## 内存与性能

常规 IEEE 754 平台上的布局由静态断言验证：

| 类型 | 字节 | 用途 |
| --- | ---: | --- |
| `vec3f` | 12 | 紧凑位置、法线、方向 |
| `vec4f` / `quatf` | 16 | 齐次坐标 / 旋转 |
| `mat3f` | 36 | 线性变换 |
| `affine3f` | 48 | 仿射变换，省略固定的最后一行 |
| `mat4f` | 64 | 投影及通用齐次变换 |

`vec3f` 对齐为 `alignof(float)`；这些类型都是标准布局、可平凡复制的值类型，不含指针、虚表或隐式额外通道。`double` 版本占用两倍空间。CPU 紧凑布局不代表满足 GPU uniform buffer 的额外对齐要求，上传时按图形 API 布局打包。

性能设计包括编译期维度、列方向连续访问、可内联的小函数、无引用计数或运行时分派，以及 SoA / mat4f SIMD 内核。mat4f 乘法在常量求值时使用通用实现，运行时使用非对齐 SIMD 加载，不增加对象对齐或尺寸。普通范围的归一化走平方范数快速路径，极大/极小输入保留缩放保护。LU 共享行缩放与主元倒数，多右端回代减少重复除法；小矩阵乘法、矩阵–向量乘法、仿射复合与四元数乘法都使用局部累加；4×4 float 矩阵–向量乘法复用矩阵乘法的 SIMD 后端。四元数插值、矩阵转换和向量间旋转复用已算结果，已经单位化的四元数插值会跳过近乎恒等的重归一化。Cholesky 和对称特征分解逐元素检查对称性而不再生成转置矩阵，Jacobi 旋转避免三角调用，`from_matrix` 直接按列检查正交性，`decompose` 用三重积而非 LU 行列式判断镜像。同一曲线或射线可通过三次多项式及射线缓存接口复用准备工作。

默认不启用 `fast-math`、全局 CPU 本机专用指令或近似倒平方根，以保留 NaN / 无穷和错误检测语义。矩阵求逆选择带选主元的 LU，优先兼顾通用性和数值稳定性；同一矩阵反复求解可复用 `factor_lu` 结果。

SIMD 路径不改变公共值类型 ABI。CMake 选项 `CHMATH_ENABLE_SIMD=OFF` 可关闭显式 SIMD；编译器仍可自动向量化普通循环。x86 SSE2 和支持 `__ARM_NEON` 的 ARM 工具链自动选择内核，其他平台执行标量实现。所有翻译单元应使用一致的宏设置；同一进程混用不同设置不受支持。

批处理接受调用方提供的内存，不分配输出。允许完整原地操作，拒绝部分重叠；SoA 的输出通道必须互不重叠。变换参数本身不能与输出缓冲区重叠。失败返回 `false`，并在写入前完成尺寸与重叠检查。

批处理入口另有可选的 workers 重载。设置 `CHMATH_ENABLE_PARALLEL=ON` 后这些重载用 `std::thread` 拆分批量；未启用、`workers <= 1` 或批量短于 `parallel_min_chunk` 时改在调用线程执行，结果完全相同。整批的尺寸与重叠契约会在任何线程启动前校验，分段按元素类型的整条缓存行对齐；平台无法创建线程时，剩余分段改由调用线程完成而不是失败。并行调用不使用全局线程池，多个调用方不会共享隐藏状态。

本机实测、基准方法和完整日志见 [docs/VALIDATION.md](docs/VALIDATION.md)。性能取决于 CPU、编译器、数据大小与缓存，不能从单个基准推导所有应用的性能。

## 错误处理与数值边界

可能失败的构造、求逆、分解和求交主要使用 `std::optional`，不分配异常对象。`try_normalize` 拒绝零向量和非有限值；便捷 `normalize(vec)` 失败时返回零向量，`normalize(quat)` 失败时返回单位四元数。

`epsilon<T>` 默认为 `32 * numeric_limits<T>::epsilon()`。`almost_equal` 同时接受有限、非负的相对与绝对容差，避免比较中间结果溢出造成误判；LU 按行缩放后检测主元，调用方可调容差。对 CAD 推荐 `double`，并按模型尺度显式设置容差。下标检查使用 `assert`，Release 下下标必须有效。普通加减乘除仍遵循 C++ 的算术规则，整数溢出、整数除零及非法下标由调用方避免。

B-spline / NURBS 是**不拥有数据的视图**，控制点、节点和权重数组必须在使用期间存活且保持不变；创建时验证输入，求值使用栈上 de Boor 工作区。次数编译期确定，限制在 0–32；NURBS 只接受严格正的有限权重。定义域外求值返回失败。内节点处采用右侧分段，定义域右端点采用左侧极限；不连续节点处导数也遵循该侧约定。

`bisect` 要求连续函数和异号括区间；`integrate` 是有深度与函数求值次数预算的自适应 Simpson 积分，误差为估计值。调用方必须检查 `converged`，回调抛出的异常按正常 C++ 规则传播，回调自行分配的内存不属于库的零分配保证。

本库提供游戏、渲染与 CAD 的通用数学基础。它不包含实体布尔运算、B-rep 拓扑、精确算术几何谓词、通用网格碰撞加速结构、稀疏/动态大矩阵、SVD 或完整曲面造型内核。浮点求交不保证所有近退化情形的精确拓扑或 watertight 三角形边界；`barycentric` 有可调相对退化容差。极端条件数、结果超出浮点范围或精度不足时仍可能失败，需要精确几何的应用应接入相应内核。

## 构建、测试和基准

最低要求：C++20 编译器、CMake 3.20；预设需要 CMake 3.21。Ninja 用于预设和验证脚本，普通 CMake 也可使用其他生成器。Windows 命令应在已加载 MSVC / Windows SDK 环境的开发者终端运行。

```sh
cmake --preset release
cmake --build --preset release --parallel
ctest --preset release
./build/release/chmath_bench
./build/release/chmath_bench_extended
./build/release/chmath_bench_parallel
```

Windows 基准可执行文件为 `build/release/chmath_bench.exe`、`chmath_bench_extended.exe` 和 `chmath_bench_parallel.exe`。可传入点数，例如 `chmath_bench 262144`。它们输出预热后九组样本的中位数和范围，包含可被编译器优化的普通 SoA 基线，不把基线强制关闭自动向量化。`chmath_bench_parallel` 始终以 `CHMATH_ENABLE_PARALLEL=1` 构建，并对比同一批量在 `workers=1` 与硬件线程数下的耗时。

使用 `parallel` 预设（`cmake --preset parallel`）或 `-DCHMATH_ENABLE_PARALLEL=ON` 编译批处理入口的 workers 重载；具体保证见「内存与性能」一节。

一键验证 Debug、Release、标量回退：

```powershell
./scripts/validate.ps1 -Compiler cl
./scripts/validate.ps1 -Compiler clang++
```

```sh
CXX=g++ sh scripts/validate.sh
CXX=clang++ sh scripts/validate.sh
```

每条命名用例在 CTest 中独立登记。CMake 会拒绝重复名称或不足 1000 条的完整测试套件，Release 中仍执行所有检查。新增参数化用例按类型、维度、角度、批量长度等实际场景独立命名；不把循环次数或同一测试在不同编译器上的重复运行计为新用例。可以只运行一类或一条测试：

```sh
ctest --test-dir build/release -L safety --output-on-failure
ctest --test-dir build/release -L stress --output-on-failure
ctest --test-dir build/release -L performance -V
./build/release/test_matrix_boundaries --case matrix_determinant_balanced_underflow_diagonal
python scripts/test_inventory.py --build-dir build/release --output build/test-inventory.json
python scripts/generate_extended_tests.py --check
```

180 条内存安全用例使用 Windows `VirtualAlloc/VirtualProtect` 或 POSIX `mmap/mprotect` 建立真实保护页并检查缓冲区契约，覆盖 SIMD 尾部、只读输入、非 SIMD 对齐地址、原地计算和失败前不写入。新增 144 条覆盖 float/double、六种操作、12 种长度。下标越界等已声明的调用方前提不会被当作合法输入执行。

60 条压力用例覆盖反复归一化、变换累积、分解复用、射线网格、样条采样、批处理流水线、病态矩阵、振荡积分和 2/4/8 线程共享只读数据。每种类型/负载有 512、4096、32768 工作单位的三个级别；积分每 128 单位执行一次有求值预算的完整积分。

42 条性能用例包含正确性对照、迭代/求值预算和预热后的计时，固定随机种子、预分配缓冲区并消费输出。默认只记录耗时；在稳定的 Release runner 上可启用 `CHMATH_PERFORMANCE_CHECKS=ON`，对其中 30 条可比较负载执行相对耗时上限。硬件负载不同的普通 CI 不启用时间阈值。Windows 可固定进程核心复现独立基准：

```powershell
./scripts/benchmark.ps1 -Executable build/release/chmath_bench.exe -Processor 0 -Repeats 3
./scripts/benchmark.ps1 -Executable build/release/chmath_bench_extended.exe -Processor 0 -Repeats 3
```

Linux / macOS 的 ASan + UBSan：

```sh
CXX=clang++ cmake --preset sanitize
cmake --build --preset sanitize --parallel
ctest --preset sanitize
```

| CMake 选项 | 默认 | 说明 |
| --- | --- | --- |
| `CHMATH_BUILD_TESTS` | 顶层工程 ON | CTest、独立头文件、多翻译单元、安装消费测试 |
| `CHMATH_BUILD_EXAMPLES` | 顶层工程 ON | 可运行的 3D 场景示例 |
| `CHMATH_BUILD_BENCHMARKS` | OFF | 无外部依赖的性能基准 |
| `CHMATH_ENABLE_SIMD` | ON | 可用平台的显式 SIMD 批处理和 mat4f 乘法 |
| `CHMATH_ENABLE_PARALLEL` | OFF | 批处理入口 workers 重载的 `std::thread` 拆分（否则串行回退） |
| `CHMATH_SANITIZERS` | OFF | Unix GCC / Clang 下启用 ASan 和 UBSan |
| `CHMATH_COVERAGE` | OFF | Clang 源码覆盖率 |
| `CHMATH_PERFORMANCE_CHECKS` | OFF | 稳定 Release runner 上启用相对性能阈值 |
| `CHMATH_WARNINGS_AS_ERRORS` | ON | 仅本项目测试、示例和基准使用严格告警 |

测试包含确定性随机性质、解析参考、边界和退化输入、极端数值、布局检查、SIMD 与标量差分、保护页/哨兵内存边界和堆分配计数。源码级覆盖率统计不会证明所有模板实例或所有输入正确；实际验证范围见报告。GitHub Actions 配置位于 `.github/workflows/ci.yml`，覆盖 Windows、Linux、macOS 和 sanitizer 构建。

详细接口索引与示例见 [docs/API.md](docs/API.md)。
