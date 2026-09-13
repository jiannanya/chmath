# chmath 接口参考

所有接口位于 `chm` 命名空间；按模块包含头文件或统一包含 `<chmath/chmath.hpp>`。以下 `T` 为浮点类型时通常选择 `float` 或 `double`。向量与矩阵不做隐式跨精度转换。

## 标量与向量

`vec<T,N>` 用 `std::array<T,N>` 紧凑存储。`vec<T,N>(s)` 将所有分量填为 s；`vec3f{1,2,3}` 按分量构造，默认构造为零。`data()` 返回连续标量数组。`vec2f/3f/4f`、`vec2d/3d/4d` 和 `vec2i/3i/4i` 为常用别名。

| 接口 | 含义 / 约束 |
| --- | --- |
| `pi<T>`, `tau<T>`, `epsilon<T>` | 编译期常量 |
| `radians`, `degrees`, `wrap_angle` | 单位转换；wrap 范围 `[-π,π)` |
| `clamp`, `saturate`, `smoothstep` | 普通 clamp 要求 lo≤hi；smoothstep 支持相等边界 |
| `lerp(a,b,t)` | 标量 / 向量插值，支持外插 |
| `almost_equal(a,b,rel,abs)` | 标量、向量或矩阵逐项比较；容差必须有限、非负；相同无穷相等，NaN 不相等 |
| `multiply_divide(a,b,c)` | 计算 a×b/c，避免可避免的中间溢出 / 下溢 |
| `dot`, `length_squared`, `distance_squared` | 快速直接运算；平方结果可能溢出 |
| `length`, `distance`, `try_normalize` | 普通范围直接计算，极大/极小输入回退到缩放范数 / 归一化 |
| `cross` | 3D 返回向量，2D 返回有向标量 |
| `hadamard`, `min`, `max`, `clamp` | 分量运算；`vector * vector` 没有含糊语义 |
| `project(vector, onto)` | 向量投影；非有限输入、onto 无法归一化或结果不可表示时失败 |
| `reflect(i,n)` | 反射，n 必须单位长度 |
| `refract(i,n,eta)` | i、n 均单位长度，eta 为有限、正的入射/出射折射率比；全反射、非法输入或非有限输出返回空 |
| `angle_between(a,b)` | `[0,π]`，退化时失败 |

基础算术与适用的组合函数使用 `constexpr`。平方根和三角函数路径在 C++20 中是运行时运算。

## 矩阵与分解

`mat<T,R,C>` 列主序存储 R×C 个元素。`mat(s)` 把主对角线填为 s，其余为零，`mat<T,N,N>::identity()` 返回单位矩阵。通过 `m(r,c)` 访问；`row` / `column` 返回副本，`set_column` 修改一列。`from_rows(std::array<vec<T,C>,R>)` 是可读的按行构造方式；标量数组构造按列主序解释。

```cpp
using namespace chm;
const auto a = mat2d::from_rows({vec2d{3, 1}, vec2d{1, 2}});
const auto x = solve(a, vec2d{9, 8}); // {2,3}
const auto lu = factor_lu(a);
if (lu) {
    const auto other = lu->solve(vec2d{4, 5}); // 已知有效分解，可反复求解
    (void)other;
}
```

| 接口 | 返回及说明 |
| --- | --- |
| `transpose`, `outer_product`, `trace` | 固定尺寸基础运算 |
| `factor_lu(A,tolerance)` | `optional<lu_factorization<T,N>>`；带行缩放的部分选主元 |
| `solve(A,b,tolerance)` | 验证 b 与有限输出，失败返回空 |
| `lu.solve(b)` | 已有分解的快速求解；要求 b 有限且输出可表示 |
| `solve(A,B,tolerance)` | B 为 N×M 矩阵，各列为一个右端；验证输入和结果，返回 `optional<mat<T,N,M>>` |
| `lu.solve(B)` | 一次处理 M 个右端，共享对角倒数；要求 B 有限、结果可表示 |
| `inverse(A,tolerance)` | 一次 LU + 多右端回代；奇异、容差内退化或非有限输出返回空 |
| `determinant(A)` | 默认不截断小主元；对中间上溢/下溢采用指数缩放重算；非有限输入返回 NaN，最终结果仍可能超出表示范围 |
| `cholesky(A,tolerance)` | 对称正定矩阵的下三角 L，满足 `A=L*transpose(L)` |
| `factor_qr(A)` | R≥C，返回完整正交 Q 与上梯形 R；Householder 算法，允许秩亏 |
| `least_squares(A,b,tolerance)` | R≥C，QR 回代求最小二乘；要求满列秩，不返回秩亏伪逆 |
| `symmetric_eigen(A,tolerance,max_sweeps)` | Jacobi 迭代；升序特征值、列特征向量、sweeps、converged |

特征向量的符号不唯一，重复特征值的基也不唯一，不应按分量精确比较。特征分解与迭代求根/积分需要同时检查返回值和 `converged`。

固定尺寸避免堆分配，但大尺寸模板会消耗栈空间和编译资源。方阵乘法、LU、Cholesky 为 O(N³)，存储为 O(N²)；QR 返回完整 R×R 的 Q，工作量 O(R²C)，不适合超大型稀疏或高瘦动态矩阵。

`mat4f * mat4f` 在支持 SSE2/NEON 的运行时使用 SIMD 内核，C++20 常量求值和其他尺寸/类型使用通用循环。对象仍只需 float 自然对齐；无需使用内部 `simd/detail/` 接口。

## 四元数与变换

`quat<T>` 默认 identity `{0,0,0,1}`，通过 `coefficients`、`imaginary()`、`x()/y()/z()/w()` 和 `data()` 访问。`q1*q2` 表示先 q2 再 q1。

| 接口 | 说明 |
| --- | --- |
| `from_axis_angle(axis,angle)` | 归一化 axis；零轴或非有限输入失败 |
| `to_axis_angle(q)` | 先归一化，输出角 `[0,π]`；identity 的轴为 +X |
| `from_euler_xyz(angles)` | 固定轴 XYZ，输入必须有限 |
| `rotation_between(a,b)` | 包含反平行处理；零向量失败 |
| `rotate(q,v)`, `to_matrix(q)` | 快速路径，要求 q 为单位四元数 |
| `from_matrix(m,tolerance)` | 验证正交性和正行列式，拒绝剪切、缩放与反射 |
| `conjugate(q)`, `inverse(q)` | 共轭；一般非单位四元数也可求逆，零四元数失败 |
| `nlerp`, `slerp` | 归一化输入、自动选择短弧；t 通常在 `[0,1]` |

`affine3<T>` 存储 3×4 矩阵，隐含第四行为 `[0,0,0,1]`。`translation(v)`、`scaling(v)`、`rotation(q)` 返回仿射变换。`compose(trs<T>{position,orientation,scale})` 为 T×R×S；`decompose` 拒绝剪切和零缩放，将反射统一放入 X 缩放。

`transform_point(affine,p)` 返回值，`transform_vector` 忽略平移；`normal_matrix` 为线性部分逆转置，`transform_normal` 随后归一化。`to_matrix(affine)` 扩展为 4×4，`to_affine(mat4)` 检查最后一行。`transform_point(mat4,p)` 执行齐次除法，w=0 或非有限结果返回空。

`look_at(eye,target,up,hand)` 返回**世界到相机**的仿射变换；眼睛与目标相同、up 平行于视线等情况失败。

```cpp
const auto projection = chm::perspective(
    chm::radians(60.0), 16.0 / 9.0, 0.1,
    std::numeric_limits<double>::infinity(),
    chm::handedness::right,
    chm::depth_range::zero_to_one,
    chm::depth_direction::reverse);
```

`perspective(fovy,aspect,near,far,hand,range,direction)` 支持正无穷 far，要求 `0<fovy<π`、aspect>0、near>0、far>near。`orthographic(left,right,bottom,top,near,far,...)` 要求所有参数有限、right>left、top>bottom、far>near。

`project(p,mvp,viewport,range)` 输出窗口坐标及 viewport 深度；它不做视锥可见性判断，可能投影位于相机后方的点。`unproject(p,inverse_mvp,viewport,range)` 执行逆操作。

`viewport<T>::valid()` 要求 x、y、width、height、min_depth、max_depth 全部有限，width/height>0 且 max_depth>min_depth；投影与反投影均执行此检查。

## 几何

`ray<T,N>`、`segment<T,N>`、`aabb<T,N>`、`sphere<T,N>` 支持通用维度（默认 3）。plane、triangle、obb、frustum 为 3D。

- `aabb` 默认 empty；`expand(point/box)` 扩展，空盒不会影响合并。`contains`、`overlaps` 包含边界。`center` 要求非空；`closest_point(box,p)` 要求有效盒。`expand(point)` 要求有限点。
- `sphere` 的 radius≥0；`obb` 的 half_extent 分量≥0，axes 列为正交单位轴。
- `plane.normal` 单位长度，方程 `dot(normal,p)+offset=0`；`signed_distance` 正值为法线侧。
- `closest_point` 支持 Segment、AABB、Plane、Triangle；Triangle 退化时退回边上的最近点。
- `area(triangle)`、`surface_area(aabb3)`、`volume(aabb3)`；平方与乘积可能溢出。
- `barycentric(triangle,p,tolerance)` 返回 a/b/c 权重；p 不在平面时对应正交投影。默认按相对面积判断退化；需要保留非常细长的三角形时可设 tolerance=0。
- `transform_bounds(affine,box)` 通过区间运算得到变换后 AABB，支持剪切和反射。

| 求交接口 | 结果 |
| --- | --- |
| `intersect(ray,aabb,t_min,t_max)` | 裁剪后的 `ray_interval{enter,exit}`，处理零/负零方向 |
| `intersect(ray,sphere,t_min,t_max)` | 体积进入 / 离开区间；射线从内部开始时 enter=t_min |
| `intersect(ray,plane,t_min,t_max,tolerance)` | 单个 t；平行或共面无唯一交点时为空 |
| `intersect(ray,triangle,t_min,t_max,cull_back_face,tolerance)` | t、重心权重、几何法线和 front_face |
| `intersect(ray,obb,t_min,t_max)` | 转到局部坐标后求区间 |
| `overlaps(obb,obb,margin)` | 15 轴 SAT；margin 是世界单位的接触余量 |

默认 t 区间为 `[0,+∞]`。返回空也可能代表非法输入或无法表示结果；接口不区分这些情况。几何对象与输入坐标应有限。

`extract_frustum(view_projection,range)` 提取向内的六个平面；零法向的无限远深度约束按 `active_mask` 忽略。平面顺序为左右、上下、较低/较高裁剪深度，反向 Z 不改变这个不等式顺序。`classify` 支持点、球、AABB，返回 outside / intersecting / inside；球和盒采用保守裁剪，可见性边界处允许额外保留对象。

`prepared_ray<T,N>::create(ray)` 验证并复制射线，缓存方向倒数。返回对象不借用原射线；`source()` 返回只读引用。`cached.intersect(box,lo,hi)` 与 `intersect(cached,box,lo,hi)` 用于同一射线反复查询 AABB。零方向单独处理；不可表示的倒数或坐标差回退到稳定除法，接近相切时使用原始射线查询复核。普通区间端点可能因倒数乘法产生末位舍入差异。

`intersect_many(cached,span<const aabb<T,N>>,span<optional<ray_interval<T>>>,lo,hi)` 写入逐盒结果；非法盒返回该元素的空结果。`classify_batch(frustum,span<const Shape>,span<containment>)` 支持点、球、AABB。两者都要求输入输出长度一致、存储不重叠；不分配内存，非法尺寸/重叠/射线区间在写入前拒绝。没有命中仍是成功完成批处理，须检查各元素结果。

## 曲线

```cpp
using namespace chm;
const std::array control{vec2d{1,0}, vec2d{1,1}, vec2d{0,1}};
const std::array<double,3> weights{1, std::sqrt(0.5), 1};
const std::array<double,6> knots{0,0,0,1,1,1};
const auto circle = nurbs_view<double,2,2>::create(control,weights,knots);
if (circle) {
    const auto midpoint = circle->evaluate(0.5); // 精确圆弧表示，结果为浮点近似
    (void)midpoint;
}
```

- `bezier<T,N,Degree>{control}`：`evaluate(t)`、`derivative()`（Degree>0）、`split(t)`；de Casteljau 算法，split 返回左右曲线，参数均重新映射为 `[0,1]`。
- `bezier.evaluate_with_derivative(t)`：返回 `curve_sample<T,N>{position,derivative}`，共享 de Casteljau 中间结果；零次曲线的导数为零。
- `cubic_polynomial<T,N>::from_bezier(bezier<T,N,3>)`：返回持有四个幂基系数的值对象，拒绝非有限控制点或不可表示的系数。提供 `evaluate`、`derivative`、`second_derivative`、`evaluate_with_derivative`，使用 Horner 求值，适合反复采样同一三次曲线。求值不检查 t 或中间范围；调用方须保证有限输入和可表示的中间量。幂基在端点附近可能出现相消，精度敏感或极端系数场景优先使用原 Bézier 求值。
- `hermite(a,tangent_a,b,tangent_b,t)`；`catmull_rom(p0,p1,p2,p3,t,tension=0)` 为均匀参数版本，插值 p1→p2。
- `bezier_patch<T,N,UDegree,VDegree>`：控制网格 `[v][u]`，`evaluate(u,v)`、`derivative_u`、`derivative_v`。
- `bspline_view<T,N,Degree>::create(control,knots)`；`evaluate(u)` 和 `derivative(u)`，`domain_min/max()`。
- `nurbs_view<T,N,Degree>::create(control,weights,knots)`；`evaluate(u)` 与 `domain_min/max()`，以齐次 de Boor 求值。
- 节点个数必须等于 `control_count+Degree+1`，非降序且有限，单个节点重数≤Degree+1，定义域长度>0；控制点至少 Degree+1 个。NURBS 权重必须有限且>0。
- 视图创建复杂度 O(control_count)，求值 O(Degree²×N)，临时存储 O(Degree×N)。Bézier 值类型直接持有控制点；借用视图不会复制整条曲线。

## 数值求根与积分

`solve_quadratic(a,b,c)` 返回零至两个升序实根；线性、重复根和恒零多项式分别处理，恒零由 `all_reals=true` 表示。无实根时返回有效但 count=0 的结果；非法输入或不可表示的根返回空。

`evaluate_polynomial(span,x)` 的系数从常数项到最高次项，空多项式返回零。

`bisect(f,lo,hi,x_tolerance,f_tolerance,max_iterations)` 返回 value、residual、iterations、converged。`integrate(f,a,b,absolute_tolerance,max_depth,max_evaluations)` 返回 value、estimated_error、evaluations、converged；支持反向区间和零长度区间，max_depth≤32。它们不建立严格数学误差界，也不保证发现采样点之间的奇点。

`compensated_sum<T>` 提供 `add`、`value`、`reset`，使用 Neumaier 补偿减少有限输入累加中的舍入误差。

## 批处理

```cpp
std::array<float, 5> x{1,2,3,4,5}, y{}, z{};
const chm::soa3<float> points{x,y,z};
const bool ok = chm::transform_points(
    chm::translation(chm::vec3f{1,2,3}), points.as_const(), points);
```

`const_soa3<T>` / `soa3<T>` 含 x/y/z 三个 span。float 使用 SSE2/NEON 四通道内核，double 使用 SSE2/AArch64 NEON 双通道内核；其他类型或不支持的平台回退到标量。`batch_backend<T>()` 返回该类型的后端名称，无模板版本返回 float 后端。ARM32 的批量归一化使用标量路径，保持普通平方根精度。

双精度点积保留可自动向量化的连续循环，不强制双通道展开；实测 Clang 可以生成更快的循环。`batch_backend<T>()` 标识该类型可用的显式批处理后端，不代表每个函数都强制使用相同指令。

| 接口 | 输入/输出与行为 |
| --- | --- |
| `transform_points(affine,input,output)` | AoS / SoA，含平移 |
| `transform_vectors(affine,input,output)` | AoS / SoA，仅线性部分 |
| `rotate_vectors(quat,input,output)` | AoS / SoA；验证单位四元数一次，共享旋转矩阵 |
| `normalize_vectors(input,output)` | AoS / SoA；零或非有限向量输出零，极端尺度稳定回退 |
| `dot_batch(a,b,span<T>)` | 两组 SoA，逐向量点积 |
| `cross_batch(a,b,soa3<T>)` | 两组 SoA，逐向量叉积 |

AoS 输入/输出类型分别为 `std::span<const vec<T,3>>` 与 `std::span<vec<T,3>>`。模板参数推导时应显式使用上述动态长度 span 类型；从 `std::array` 推导出的静态长度 span 可能需要显式转换。

所有 span 长度必须一致，返回 bool 表示尺寸与重叠契约是否成立。输入和输出可完全重合或完全分离，不能部分重叠；每个输出 SoA 通道相互分离。完整通道置换也可原地执行。旋转额外检查四元数有限且单位长度，归一化逐向量处理非法数值；其他函数不逐个检查分量是否有限，结果遵循普通浮点算术。变换参数对象不能与输出存储重叠。span 必须引用存活、可访问且正确对齐的对象；长度/重叠检查不能验证悬空指针或调用方虚构的容量。

库没有可变全局缓存。不同线程可共享不可变值对象、分解和曲线视图，但视图借用的数据也必须保持不变；并发输出需使用互不重叠的缓冲区。内存安全和压力测试覆盖这些约定，不保证违反前提的调用安全。
