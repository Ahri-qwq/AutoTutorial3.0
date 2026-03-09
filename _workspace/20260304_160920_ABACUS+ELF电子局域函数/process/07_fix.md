# 前言

本教程介绍如何使用 ABACUS 计算并可视化电子局域函数（Electron Localization Function，ELF）。

**教程目标：** 掌握 ABACUS 中 ELF 的输出设置，能对典型体系（分子、磁性金属）完成计算并用 VESTA 可视化。

**适合读者：** 已会用 ABACUS 做基本 SCF 计算，希望进一步分析化学键和电子局域化性质的用户。

**前置知识：** ABACUS 的 INPUT/STRU/KPT 基本格式，SCF 计算流程。

**教程结构：**
- 第一章：ELF 简介，介绍物理含义和 ABACUS 支持的基组
- 第二章：参数设置，说明 `out_elf` 参数和输出文件
- 第三章：案例一，水分子（自旋非极化，PW 和 LCAO 基组）
- 第四章：案例二，BCC 铁（自旋极化，PW 基组）
- 第五章：常见问题与注意事项
- 附录：参数速查、参考文献
# 一、ELF 简介

电子局域函数（Electron Localization Function，ELF）是实空间中描述电子局域化程度的标量场，最早由 Becke 和 Edgecombe 在 1990 年提出 [1]。

**数值含义：**

ELF 的取值范围是 0 到 1：

| ELF 值 | 含义 |
|--------|------|
| 接近 1 | 电子高度局域化（共价键区域、孤对电子） |
| = 0.5 | 均匀电子气（与金属性区域接近） |
| 接近 0 | 电子离域（金属键、真空区域） |

通过可视化 ELF 的三维分布，可以直观判断体系中的化学键类型——共价键、孤对电子、金属键各有特征性的 ELF 图案。

**ABACUS 支持情况：**

ABACUS 支持在以下三种计算模式下输出 ELF：
- 平面波基组（PW）的 Kohn-Sham DFT
- 原子轨道基组（LCAO）的 Kohn-Sham DFT
- 无轨道密度泛函理论（OFDFT）

本教程以 PW 和 LCAO 两种常用基组为主，通过两个典型案例介绍完整的计算和可视化流程。
# 二、参数设置

## 2.1 核心参数

在 INPUT 文件中，只需添加一行即可启用 ELF 输出：

```
out_elf 1 3
```

参数 `out_elf` 接受两个整数：

| 位置 | 含义 | 默认值 |
|------|------|--------|
| 第 1 个值 | 是否输出 ELF（1=输出，0=不输出） | 0 |
| 第 2 个值 | 输出的有效数字位数 | 3 |

ELF 计算不需要额外的输入文件，也不影响原有的 SCF 迭代流程——只需在正常的自洽计算 INPUT 中加上这一行即可。

## 2.2 输出文件

计算结束后，ELF 存储在 `OUT.${suffix}/` 目录下，具体文件取决于 `nspin` 的设置：

| nspin | 输出文件 | 内容 |
|-------|---------|------|
| 1 | `ELF.cube` | 总 ELF |
| 2 | `ELF.cube` | 总 ELF |
| 2 | `ELF_SPIN1.cube` | 自旋向上电子的 ELF |
| 2 | `ELF_SPIN2.cube` | 自旋向下电子的 ELF |

这三个文件都是标准的 Gaussian CUBE 格式，可以直接用 VESTA 打开。

## 2.3 cube 文件格式

ELF.cube 采用标准 Gaussian CUBE 格式，结构如下：

```
第 1-2 行   注释行（VESTA 不读取）
第 3 行     原子数  格点原点坐标（Bohr）
第 4-6 行   各方向 FFT 格点数  步长向量（Bohr）
第 7-N 行   原子信息（原子序数  价电子数  x  y  z，Bohr）
后续行      ELF 数据（z 变化最快，每行 6 个数据）
```

具体输出格式可参考第三章3.1节中的实际案例。
# 三、案例一：水分子 ELF（自旋非极化）

本案例对真空中的水分子（H₂O）计算 ELF，分别使用平面波（PW）和原子轨道（LCAO）两种基组。

算例文件：
- PW：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/H20-pw
- LCAO：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/H20-lcao

## 3.1 平面波计算

### 输入文件

**INPUT：**

```
INPUT_PARAMETERS
suffix              autotest
calculation         scf
ntype               2
basis_type          pw
pseudo_dir          ./

ecutwfc             100          # 截断能，单位 Ry
scf_thr             1e-6
scf_nmax            100

nspin               1            # 自旋非极化
symmetry            0

out_elf             1 3          # 输出 ELF，保留 3 位有效数字
```

**STRU：**

水分子置于约 14.82 Å 的真空立方盒内，分子坐标来自 ELF.cube 头部数据（单位 Bohr）：

```
ATOMIC_SPECIES
H   1.008    H_ONCV_PBE-1.0.upf
O  15.999    O_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.889726         # Bohr/Å

LATTICE_VECTORS
14.82   0.00    0.00
 0.00  14.82    0.00
 0.00   0.00   14.82

ATOMIC_POSITIONS
Direct

H
0.0
2
0.4274  0.6688  0.2972  0  0  0    # H1
0.5087  0.7379  0.2679  0  0  0    # H2

O
0.0
1
0.4841  0.7029  0.3264  0  0  0    # O
```

**KPT：**

水分子在真空盒内，只需 Γ 点采样：

```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

### 运行计算

```bash
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

### 输出文件

计算结束后，在 `OUT.autotest/` 目录下生成 `ELF.cube`。文件头部如下：

```
STEP: 0  Cubefile created from ABACUS. Inner loop is z, followed by y and x
1 (nspin)
3 0.0 0.0 0.0
180 0.155556 0.000000 0.000000
180 0.000000 0.155556 0.000000
180 0.000000 0.000000 0.155556
 1 1.000000 11.967714 18.726812 8.322864
 1 1.000000 14.244302 20.657387 7.503683
 8 6.000000 13.550429 19.681072 9.139011
 0.000e+00 0.000e+00 0.000e+00 0.000e+00 ...
```

- 180×180×180 FFT 格点，步长 0.155556 Bohr（约 0.082 Å）
- 对应约 28 Bohr（~14.8 Å）的真空盒子
- 三行原子信息：H（序数 1）、H（序数 1）、O（序数 8），坐标单位 Bohr
- 远离分子的真空区域 ELF≈0

### VESTA 可视化

**等高面图（3D）：**

1. 用 VESTA 打开 `ELF.cube`（File → Open）
2. 菜单 Properties → Isosurfaces
3. 添加等值面，推荐 ELF = 0.75（捕捉成键区域）
4. 调整颜色和透明度

**截面图（2D）：**

1. Properties → Sections/Isosurfaces
2. 选择 Section
3. 勾选 "Specify from three atoms"，选取 H1、O、H2 三个原子，得到 H-O-H 平面截面

**可视化结果解读：**

| 区域 | ELF 值 | 物理含义 |
|------|--------|---------|
| O-H 键区（两核之间） | ~0.7–0.8 | 共价键电子，中等局域化 |
| O 原子两侧（孤对电子） | 接近 1 | 高度局域的孤对电子 |
| 远离分子的真空区域 | ≈ 0 | 无电子分布 |

水分子是典型的极性共价分子，O 原子有两个孤对电子，ELF 图像可以直观地显示出这两个孤对电子的"兔耳"形分布。

## 3.2 原子轨道计算（LCAO）

算例地址：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/H20-lcao

### 与 PW 的差异

LCAO 计算只需在 INPUT 和 STRU 中做少量修改，STRU 中的坐标与 PW 相同：

**INPUT（LCAO 版本）：**

```
INPUT_PARAMETERS
suffix              autotest
calculation         scf
ntype               2
basis_type          lcao          # 改为 LCAO
pseudo_dir          ./
orbital_dir         ./            # 新增：指定轨道文件目录

ecutwfc             100
scf_thr             1e-6
scf_nmax            100

nspin               1
symmetry            0

out_elf             1 3
```

**STRU（LCAO 版本，新增 NUMERICAL_ORBITAL 段）：**

```
ATOMIC_SPECIES
H   1.008    H_ONCV_PBE-1.0.upf
O  15.999    O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_6au_100Ry_2s1p.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
14.82   0.00    0.00
 0.00  14.82    0.00
 0.00   0.00   14.82

ATOMIC_POSITIONS
Direct

H
0.0
2
0.4274  0.6688  0.2972  0  0  0
0.5087  0.7379  0.2679  0  0  0

O
0.0
1
0.4841  0.7029  0.3264  0  0  0
```

轨道文件命名含义：以 `O_gga_7au_100Ry_2s2p1d.orb` 为例：
- `gga`：GGA 泛函
- `7au`：截断半径 7 Bohr
- `100Ry`：与 ecutwfc = 100 Ry 匹配
- `2s2p1d`：2 个 s 轨道，2 个 p 轨道，1 个 d 轨道

### LCAO 稳定性修正

LCAO 基组存在不完备性，在远离原子核的区域，动能密度的数值行为可能出现异常，导致真空区出现非物理的非零 ELF。ABACUS 采用文献 [2] 的建议，在动能密度的分子部分加入 10⁻⁵ 的小量修正。该修正确保了远场行为的正确性，同时不影响近核区域的结果。PW 基组没有此问题。

### PW vs LCAO 结果对比

LCAO 计算得到的 ELF 与 PW 结果大体一致，O-H 键区和孤对电子均有清晰的高 ELF 峰，但由于基组不同，细节上存在差异：

- 化学键区（近核区域）：两种基组结果接近
- 真空区域：PW 结果干净（ELF≈0）；LCAO 经修正后也正常
- 孤对电子形状：LCAO 略有差异，依赖基组大小（DZP vs TZP）

对于精确的 ELF 分析，PW 基组更可靠；LCAO 计算速度更快，适合大体系的定性分析。
# 四、案例二：BCC 铁 ELF（自旋极化）

算例地址：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/bcc-Fe-pw

体心立方铁（BCC Fe）是铁磁性金属，使用 PW 基组计算自旋极化的 ELF，可以得到自旋向上和自旋向下电子各自的局域化图像。

## 4.1 输入文件

**INPUT：**

```
INPUT_PARAMETERS
suffix              autotest
calculation         scf
ntype               1
basis_type          pw
pseudo_dir          ./

ecutwfc             100          # 截断能，单位 Ry
scf_thr             1e-6
scf_nmax            200

nspin               2            # 自旋极化（共线磁矩）
smearing_method     gauss
smearing_sigma      0.01

mixing_type         broyden
mixing_beta         0.7

out_elf             1 3          # 输出 ELF（总/自旋↑/自旋↓）
```

**STRU：**

BCC Fe 传统晶胞含 2 个原子（格点参数 a = 2.866 Å），初始磁矩按元素设置为 2.3 μB：

```
ATOMIC_SPECIES
Fe  55.845  Fe_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.889726         # Bohr/Å

LATTICE_VECTORS
2.866  0.000  0.000
0.000  2.866  0.000
0.000  0.000  2.866

ATOMIC_POSITIONS
Direct

Fe
2.3              # 初始磁矩，单位 μB
2
0.000  0.000  0.000  m  1  1  1
0.500  0.500  0.500  m  1  1  1
```

初始磁矩设置为 2.3 μB（接近铁的实验磁矩 ~2.2 μB）。设置合理的初始磁矩对收敛到铁磁基态至关重要——过小的初始磁矩可能导致 SCF 收敛到非磁态。

**KPT：**

```
K_POINTS
0
Gamma
8 8 8 0 0 0
```

## 4.2 运行与输出文件

```bash
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

`nspin = 2` 时，ABACUS 输出三个 cube 文件：

| 文件名 | 内容 |
|--------|------|
| `ELF.cube` | 总 ELF（自旋向上+向下的综合） |
| `ELF_SPIN1.cube` | 自旋向上（↑）电子的 ELF |
| `ELF_SPIN2.cube` | 自旋向下（↓）电子的 ELF |

## 4.3 VESTA 可视化

对每个 cube 文件分别在 VESTA 中打开，查看 (100) 晶面截面图：

1. File → Open → 选择 `ELF.cube`（或 `ELF_SPIN1.cube` / `ELF_SPIN2.cube`）
2. Properties → Sections
3. 设置截面方向：选择 (100) 平面（设置法向量为 [1, 0, 0]）
4. 调整颜色映射（推荐蓝-白-红方案，蓝=0，红=1）

## 4.4 结果解读

对三种 ELF 的 (100) 面截面图：

**总 ELF（ELF.cube）：**
- Fe 原子核附近 ELF 接近 1（强烈局域的芯电子/内层电子）
- 原子间区域 ELF 较低（金属 d 电子相对离域）
- 整体图像具有 BCC 晶体的四重对称性

**自旋向上 ELF（ELF_SPIN1.cube）：**
- 自旋向上电子在 Fe 原子附近的局域化分布
- 反映了 d 轨道中自旋向上电子的空间分布特征

**自旋向下 ELF（ELF_SPIN2.cube）：**
- 与自旋向上 ELF 形状相似，但数值略有不同
- 自旋向上和向下的 ELF 差异体现了铁磁性——自旋向上电子更多（磁矩约 2.2 μB），其 ELF 峰值略高

BCC Fe 与水分子的对比：金属体系的 ELF 在原子间区域通常在 0.3–0.5 之间（类均匀电子气），而分子的共价键区域 ELF 可达 0.7 以上。这种对比直观地说明了金属键与共价键在电子局域化程度上的本质差异。
# 五、常见问题与注意事项

## 5.1 LCAO 计算中真空区出现非零 ELF

**现象：** 在分子或表面体系的真空区，LCAO 计算的 ELF 不为零，出现"虚假"的局域化。

**原因与处理：** 参见第三章3.2节"LCAO 稳定性修正"。ABACUS 已内置修正，通常无需手动处理；若仍有异常，切换到 PW 基组验证。

## 5.2 out_elf 第二个参数的选择

`out_elf 1 3` 中第二个参数控制输出精度（有效数字位数）：

- `3`（默认）：cube 文件较小，足够用于可视化
- `6` 或更高：适合需要精确数值的定量分析

精度越高，cube 文件越大，对于大体系（密格点、大晶胞）建议保持默认值 3。

## 5.3 如何选择 PW 还是 LCAO

| 场景 | 推荐基组 | 理由 |
|------|---------|------|
| 分子/表面体系精确分析 | PW | 无基组不完备问题，真空区干净 |
| 大体系定性分析 | LCAO | 计算成本低，ELF 近核区结果可靠 |
| 含金属元素的周期体系 | PW 或 LCAO 均可 | 验证两种基组结果一致性 |

## 5.4 ELF 计算对 SCF 参数无特殊要求

ELF 是在 SCF 收敛后从电荷密度计算得到的，不影响自洽迭代本身。因此 `out_elf 1` 可以添加到任何正常的 SCF 计算 INPUT 中，不需要修改其他参数。

## 5.5 自旋极化计算的注意事项

- 初始磁矩设置不合理，SCF 可能收敛到非磁态（磁矩接近 0）
- 建议参考实验值或文献设置初始磁矩（如 Fe 设为 2.3 μB）
- `smearing_method gauss` 配合小的 `smearing_sigma` 有助于磁性金属收敛
# 附录

## A. 参数速查表

| 参数 | 格式 | 默认值 | 说明 |
|------|------|--------|------|
| `out_elf` | int int | 0 3 | 第1值：1=输出；第2值：有效数字位数 |
| `nspin` | int | 1 | 1=非极化，2=共线磁矩 |

**输出文件速查：**

| nspin | 输出文件 | 内容 |
|-------|---------|------|
| 1 | `OUT.*/ELF.cube` | 总 ELF |
| 2 | `OUT.*/ELF.cube` | 总 ELF |
| 2 | `OUT.*/ELF_SPIN1.cube` | 自旋向上 ELF |
| 2 | `OUT.*/ELF_SPIN2.cube` | 自旋向下 ELF |

## B. 参考文献

[1] Becke A D, Edgecombe K E. A simple measure of electron localization in atomic and molecular systems[J]. The Journal of Chemical Physics, 1990, 92(9): 5397-5403.

[2] 卢天, 陈飞武. 电子定域化函数的含义与函数形式[J]. 物理化学学报, 2011, 27(12): 2786-2792.

## C. 进阶学习

- **Bader 电荷分析**：与 ELF 互补，从积分角度定量分析原子间电荷转移，ABACUS 支持输出 Bader 分析所需的电荷密度，参考 [ABACUS+Bader charge 分析教程](https://mcresearch.github.io/abacus-user-guide/)
- **差分电荷密度**：计算成键前后电荷密度之差，可配合 ELF 共同分析化学键特征
- **OFDFT 中的 ELF**：无轨道密度泛函理论（OFDFT）同样支持 `out_elf`，适用于大体系（数千原子量级）的快速扫描
