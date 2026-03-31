---
title: "ABACUS 力学性质计算"
author: "AutoTutorial 3.0"
date: "2026-03-30"
topic: "ABACUS力学性质计算"
task_type: "B"
has_case: true
draft_version: "initial"
---

# 前言

本教程是 ABACUS 系列培训的第二篇，聚焦材料的力学性质计算。读完本篇，你将掌握弹性常数和弹性模量的计算方法，能够独立完成从结构优化到弹性张量提取的完整流程。

## 适用读者

- 已完成第一篇《ABACUS 基本介绍》的学习
- 了解结构优化（relax/cell-relax）的基本概念
- 希望计算材料力学性质（弹性模量、硬度等）的研究人员

## 前置知识

- DFT 基础概念（平面波、赝势、自洽迭代）
- ABACUS 输入文件准备（INPUT、STRU、KPT）
- 结构优化的基本流程（第一篇已讲解）
- 基本的 Python 使用（用于运行后处理脚本）

## 本篇结构

| 章节 | 内容 |
|------|------|
| 第一章 | 弹性常数基础理论 |
| 第二章 | 弹性常数计算方法 |
| 第三章 | Al 弹性常数计算实战（详细案例）|
| 第四章 | 使用 pymatgen 自动化计算 |
| 第五章 | 使用 abacustest 简化流程 |

## 与第一篇的衔接

第一篇教程介绍了结构优化（relax/cell-relax）的基本概念和操作。本篇将在此基础上，利用结构优化功能计算材料的弹性性质。我们会简要回顾结构优化的关键参数，但不再重复详细讲解文件格式。

## 后续篇章

本系列共四篇：
1. ABACUS 基本介绍（已完成）
2. **力学性质计算**（本篇）
3. 能带结构计算
4. 态密度（DOS/PDOS）计算

建议按顺序阅读，每篇均以前篇为基础。
# 第一章：弹性常数基础理论

## 1.1 什么是弹性常数

弹性常数（Elastic Constants）是表征材料弹性的物理量，描述了材料在外力作用下抵抗可逆形变的能力。这是材料的重要力学性质，可以通过实验测量或理论计算获得。

在材料科学中，弹性常数有多个重要应用：
- 判断晶体的力学稳定性
- 预测材料的硬度和刚度
- 辅助新材料的计算设计
- 理解材料的各向异性特征

### 弹性变形与塑性变形

材料受力后的变形分为两类：

**弹性变形：** 外力撤除后，材料恢复原状。应力与应变呈线性关系，符合胡克定律。

**塑性变形：** 外力撤除后，材料无法恢复，发生永久变形。超出弹性极限后发生。

弹性常数描述的是材料在弹性变形范围内的行为。

## 1.2 应力与应变

### 应力张量

应力（Stress）σ 描述材料内部单位面积上受到的力，单位为 Pa（帕斯卡）或 GPa。

应力是一个二阶张量，在三维空间中可以表示为 3×3 矩阵：

```
σ = | σ_xx  σ_xy  σ_xz |
    | σ_yx  σ_yy  σ_yz |
    | σ_zx  σ_zy  σ_zz |
```

其中：
- 对角元素（σ_xx, σ_yy, σ_zz）：正应力，描述拉伸或压缩
- 非对角元素（σ_xy, σ_xz, σ_yz）：剪切应力

由于应力张量是对称的（σ_ij = σ_ji），实际上只有 6 个独立分量。

### 应变张量

应变（Strain）ε 描述材料的相对形变，是无量纲量。

应变也是一个二阶对称张量：

```
ε = | ε_xx  ε_xy  ε_xz |
    | ε_yx  ε_yy  ε_yz |
    | ε_zx  ε_zy  ε_zz |
```

其中：
- 对角元素：正应变，描述长度的相对变化
- 非对角元素：剪切应变，描述角度的变化

### 胡克定律

在线弹性范围内，应力与应变之间满足胡克定律（Hooke's Law）：

```
σ_ij = C_ijkl ε_kl
```

其中 C_ijkl 是四阶弹性刚度张量（Elastic Stiffness Tensor），即弹性常数。这个张量有 3^4 = 81 个分量，但由于对称性，独立分量大大减少。

## 1.3 Voigt 表示法

### 为什么需要 Voigt 表示法

四阶张量 C_ijkl 有 81 个分量，书写和计算都很不方便。由于应力和应变张量都是对称的，我们可以利用 Voigt 表示法将它们简化为向量形式。

### 指标映射规则

Voigt 表示法将二阶张量的两个指标映射为一个指标：

```
xx → 1
yy → 2
zz → 3
yz (或 zy) → 4
xz (或 zx) → 5
xy (或 yx) → 6
```

### 应力和应变的向量表示

应力张量简化为 6 维向量：

```
σ_1 = σ_xx
σ_2 = σ_yy
σ_3 = σ_zz
σ_4 = σ_yz
σ_5 = σ_xz
σ_6 = σ_xy
```

应变张量简化为 6 维向量（注意剪切应变的因子 2）：

```
ε_1 = ε_xx
ε_2 = ε_yy
ε_3 = ε_zz
ε_4 = 2ε_yz
ε_5 = 2ε_xz
ε_6 = 2ε_xy
```

### 弹性常数矩阵

使用 Voigt 表示法后，胡克定律变为矩阵形式：

```
| σ_1 |   | C_11  C_12  C_13  C_14  C_15  C_16 |   | ε_1 |
| σ_2 |   | C_12  C_22  C_23  C_24  C_25  C_26 |   | ε_2 |
| σ_3 | = | C_13  C_23  C_33  C_34  C_35  C_36 | × | ε_3 |
| σ_4 |   | C_14  C_24  C_34  C_44  C_45  C_46 |   | ε_4 |
| σ_5 |   | C_15  C_25  C_35  C_45  C_55  C_56 |   | ε_5 |
| σ_6 |   | C_16  C_26  C_36  C_46  C_56  C_66 |   | ε_6 |
```

弹性常数矩阵 C 是对称的（C_ij = C_ji），因此最多有 21 个独立分量。

## 1.4 晶体对称性与独立弹性常数

晶体的对称性会进一步减少独立弹性常数的数量。不同晶系的独立分量数如下：

| 晶系 | 独立弹性常数数量 | 典型材料 |
|------|-----------------|---------|
| 三斜（Triclinic） | 21 | CuSO₄·5H₂O |
| 单斜（Monoclinic） | 13 | 石膏 |
| 正交（Orthorhombic） | 9 | 硫磺 |
| 四方（Tetragonal） | 6 或 7 | TiO₂ |
| 六方（Hexagonal） | 5 | 石墨、Mg |
| 三方（Trigonal） | 6 或 7 | 石英 |
| 立方（Cubic） | 3 | Si、Al、NaCl |

### 立方晶系的弹性常数

立方晶系（如 FCC、BCC、Diamond 结构）具有最高的对称性，只有 3 个独立弹性常数：C_11、C_12、C_44。

弹性常数矩阵简化为：

```
| C_11  C_12  C_12   0     0     0   |
| C_12  C_11  C_12   0     0     0   |
| C_12  C_12  C_11   0     0     0   |
|  0     0     0    C_44   0     0   |
|  0     0     0     0    C_44   0   |
|  0     0     0     0     0    C_44 |
```

本教程的案例 Al（铝）是 FCC 结构，属于立方晶系。

## 1.5 弹性模量

弹性常数矩阵 C_ij 包含了材料弹性性质的完整信息，但不够直观。工程上常用几个标量的弹性模量来描述材料的力学性质。

### 体弹模量（Bulk Modulus, B）

体弹模量描述材料抵抗体积变化（均匀压缩或膨胀）的能力。定义为：

```
B = -V (dP/dV)
```

其中 V 是体积，P 是压强。体弹模量越大，材料越难被压缩。

对于立方晶系，体弹模量可以从弹性常数计算：

```
B = (C_11 + 2C_12) / 3
```

单位：GPa（吉帕）

### 剪切模量（Shear Modulus, G）

剪切模量描述材料抵抗剪切形变的能力，也称为刚度模量。剪切模量越大，材料越不容易发生剪切变形。

对于立方晶系：

```
G = (C_11 - C_12 + 3C_44) / 5  （Voigt 平均）
```

单位：GPa

### 杨氏模量（Young's Modulus, E）

杨氏模量描述材料在单轴拉伸或压缩时的刚性，定义为正应力与正应变的比值。杨氏模量越大，材料越不容易发生形变。

可以从体弹模量和剪切模量计算：

```
E = 9BG / (3B + G)
```

单位：GPa

### 泊松比（Poisson's Ratio, ν）

泊松比描述材料在一个方向受拉伸时，在垂直方向上收缩的程度。定义为横向应变与轴向应变之比的负值。

```
ν = (3B - 2G) / (6B + 2G)
```

泊松比是无量纲量，通常在 0 到 0.5 之间。大多数金属的泊松比在 0.25 到 0.35 之间。

### 各弹性模量的物理意义对比

| 弹性模量 | 物理意义 | 典型值（Al） |
|---------|---------|-------------|
| 体弹模量 B | 抵抗体积变化 | ~76 GPa |
| 剪切模量 G | 抵抗剪切形变 | ~26 GPa |
| 杨氏模量 E | 抵抗拉伸压缩 | ~70 GPa |
| 泊松比 ν | 横向收缩程度 | ~0.35 |

# 第二章：弹性常数计算方法

## 2.1 两种计算方法概述

从第一性原理（DFT）计算弹性常数，主要有两种方法：

### 能量-应变法（Energy-Strain Method）

**原理：** 将体系的总能量 E 按应变 ε 进行泰勒展开：

```
E(ε) = E(0) + (∂E/∂ε)ε + (1/2)(∂²E/∂ε²)ε² + ...
```

在平衡位置，一阶导数为零。弹性常数 C 与能量对应变的二阶导数相关：

```
C_ij = (1/V₀) ∂²E / ∂ε_i∂ε_j
```

**计算流程：**
1. 对结构施加不同大小的应变
2. 计算每个应变下的总能量
3. 拟合能量-应变曲线（抛物线）
4. 从曲线的曲率（二阶导数）得到弹性常数

**优点：**
- 算法简单，只需计算总能量
- 不需要计算应力张量

**缺点：**
- 需要拟合抛物线，计算点数较多
- 对能量收敛精度要求极高（通常需要 1e-8 eV 或更高）
- 对数值噪声敏感

### 应力-应变法（Stress-Strain Method）

**原理：** 直接利用应力 σ 与应变 ε 的线性关系（胡克定律）：

```
σ_i = C_ij ε_j
```

弹性常数 C 是应力对应变的一阶导数：

```
C_ij = ∂σ_i / ∂ε_j
```

**计算流程：**
1. 对结构施加不同大小的应变
2. 计算每个应变下的应力张量
3. 对应力-应变数据进行线性拟合
4. 拟合直线的斜率即为弹性常数

**优点：**
- 高信息密度：单次计算输出 6 个应力分量，而能量法只输出 1 个标量
- 计算效率高：所需应变构型数量少
- 数值稳定性好：线性拟合比二次拟合更稳定
- 对数值噪声容忍度高

**缺点：**
- 要求软件能准确计算应力张量
- 对 k 点采样和截断能的要求较高（尤其是金属）

### 方法选择

**本教程采用应力-应变法**，原因：
1. Materials Project 等主流数据库的标准方法
2. 计算效率高，适合高通量计算
3. ABACUS 对应力计算支持良好
4. pymatgen 工具默认使用此方法

## 2.2 应力-应变法详解

### 基本思路

对于立方晶系（如 Al），只有 3 个独立弹性常数：C_11、C_12、C_44。我们需要设计应变模式，使得每种应变能够独立地探测某个或某几个弹性常数。

### 应变模式设计

为了全面探测弹性常数，需要施加 6 种独立的应变状态：

| 应变类型 | 应变向量 [ε₁, ε₂, ε₃, ε₄, ε₅, ε₆] | 探测的弹性常数 |
|---------|-----------------------------------|---------------|
| 正应变 x | [δ, 0, 0, 0, 0, 0] | C_11, C_12, C_13 |
| 正应变 y | [0, δ, 0, 0, 0, 0] | C_12, C_22, C_23 |
| 正应变 z | [0, 0, δ, 0, 0, 0] | C_13, C_23, C_33 |
| 剪切应变 yz | [0, 0, 0, γ, 0, 0] | C_44 |
| 剪切应变 xz | [0, 0, 0, 0, γ, 0] | C_55 |
| 剪切应变 xy | [0, 0, 0, 0, 0, γ] | C_66 |

其中 δ 是正应变大小，γ 是剪切应变大小。

### 应变大小的选择

应变不能太大（超出线性范围），也不能太小（数值噪声影响大）。通常选择：

```
δ = ±0.5%, ±1.0%  （即 ±0.005, ±0.01）
γ = ±0.5%, ±1.0%
```

对每种应变状态，施加 4 种不同大小的应变：

```
norm_strains = [-0.010, -0.005, 0.005, 0.010]
shear_strains = [-0.010, -0.005, 0.005, 0.010]
```

6 种应变状态 × 4 种应变大小 = **24 个应变构型**

### 线性拟合原理

对于每种应变状态，我们有 4 个（应变，应力）数据点。通过线性拟合：

```
σ = C × ε
```

拟合直线的斜率即为弹性常数 C。

使用多个数据点拟合的好处：
- 减小数值噪声的影响
- 验证线性关系是否成立
- 提高结果的可靠性

### 应变矩阵的构造

对于三维晶体，应变矩阵为：

```
ε = | ε₁   ε₆/2  ε₅/2 |
    | ε₆/2  ε₂   ε₄/2 |
    | ε₅/2  ε₄/2  ε₃  |
```

变形后的晶格矢量为：

```
A' = A · (I + ε)
```

其中 A 是初始晶格矢量，I 是单位矩阵。

## 2.3 计算流程概述

完整的弹性常数计算流程包括以下步骤：

**步骤 1：结构优化（cell-relax）**
- 目的：得到平衡晶格参数和原子位置
- 计算类型：`calculation = cell-relax`
- 输出：优化后的 STRU 文件

**步骤 2：生成应变构型**
- 目的：基于平衡结构生成 24 个应变构型
- 工具：pymatgen 脚本或 abacustest
- 输出：24 个 task 文件夹，每个包含应变后的 STRU

**步骤 3：计算应力（relax）**
- 目的：对每个应变构型计算应力张量
- 计算类型：`calculation = relax`（固定晶格，允许原子弛豫）
- 输出：每个 task 的应力数据

**步骤 4：拟合弹性常数**
- 目的：从应力-应变数据拟合弹性常数
- 工具：pymatgen 或 abacustest
- 输出：弹性常数矩阵、弹性模量

### 为什么步骤 3 使用 relax 而不是 scf

在应变后的结构中，原子位置可能偏离新的平衡位置。允许原子弛豫（relax）可以：
- 消除内应力
- 得到更准确的应力张量
- 符合实验条件（准静态过程）

这种方法称为"弛豫离子模型"（relaxed-ion model），是计算弹性常数的标准做法。

# 第三章：Al 弹性常数计算实战

本章以铝（Al）的 FCC 结构为例，详细演示弹性常数计算的完整流程。Al 是典型的金属材料，属于立方晶系，只有 3 个独立弹性常数。

## 3.1 案例介绍

### 材料体系

**材料：** 铝（Aluminum, Al）
**晶体结构：** FCC（面心立方）
**空间群：** Fm-3m (225)
**晶格常数：** a ≈ 4.05 Å（实验值）
**原子数：** 4 个 Al 原子（FCC 原胞）

### 计算参数

**基组类型：** 平面波（pw）
**截断能：** 70 Ry
**k 点网格：** 6×6×6
**赝势：** Al.upf（ONCV 模守恒赝势）
**展宽方法：** Gaussian，σ = 0.015 Ry（金属体系）

### 两步计算流程

**步骤 1：晶胞优化（cell-relax）**
- 同时优化晶格参数和原子位置
- 得到平衡结构

**步骤 2：应变构型计算（relax）**
- 基于平衡结构生成 24 个应变构型
- 对每个构型计算应力
- 拟合弹性常数

## 3.2 步骤 1：晶胞优化（cell-relax）

### 3.2.1 INPUT 文件

```
INPUT_PARAMETERS

# System variables
calculation         cell-relax
symmetry            1
kspacing            0.14
precision           double

# Plane wave related variables
ecutwfc             70.0
pw_diag_nmax        20
pw_diag_ndim        2

# Electronic structure
basis_type          pw
ks_solver           dav_subspace
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-08

# Geometry relaxation
relax_method        cg
relax_nmax          60
cal_force           1
force_thr_ev        0.01
cal_stress          1
stress_thr          0.5
fixed_axes          None
```

**关键参数说明：**

- `calculation = cell-relax`：同时优化晶格和原子位置
- `kspacing = 0.14`：自动生成 k 点网格（单位：1/Bohr），对于 Al 约等于 6×6×6
- `ecutwfc = 70.0`：平面波截断能 70 Ry，对于 Al 已足够收敛
- `scf_thr = 1e-08`：高精度收敛阈值，确保应力计算准确
- `smearing_method = gauss`：Gaussian 展宽，适合金属
- `smearing_sigma = 0.015`：展宽宽度 0.015 Ry
- `mixing_type = broyden`：Broyden 混合方法，收敛性好
- `force_thr_ev = 0.01`：力收敛阈值 0.01 eV/Å
- `stress_thr = 0.5`：应力收敛阈值 0.5 kBar

### 3.2.2 STRU 文件

```
ATOMIC_SPECIES
Al 26.981500 Al.upf

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    4.05000000000     0.00000000000     0.00000000000
    0.00000000000     4.05000000000     0.00000000000
    0.00000000000     0.00000000000     4.05000000000

ATOMIC_POSITIONS
Cartesian

Al
0.000000
4
    0.00000000000     0.00000000000     0.00000000000 1 1 1
    0.00000000000     2.02500000000     2.02500000000 1 1 1
    2.02500000000     0.00000000000     2.02500000000 1 1 1
    2.02500000000     2.02500000000     0.00000000000 1 1 1
```

**文件说明：**

- `LATTICE_CONSTANT = 1.889726`：转换因子，使得 LATTICE_VECTORS 单位为 Å
- 晶格矢量：立方晶格，a = 4.05 Å（初始值，接近实验值）
- 坐标系：Cartesian（笛卡尔坐标）
- 4 个 Al 原子：FCC 结构的原胞
- `1 1 1`：允许原子在三个方向上移动

### 3.2.3 KPT 文件

```
K_POINTS
0
Gamma
6 6 6 0 0 0
```

**说明：**
- 使用 Gamma 中心的 Monkhorst-Pack 网格
- 6×6×6 k 点网格，对于 Al 金属已足够收敛

> **注意：** 也可以使用 INPUT 中的 `kspacing = 0.14` 自动生成 k 点，此时 KPT 文件内容会被忽略。

### 3.2.4 运行计算

```bash
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

计算完成后，优化后的结构保存在 `OUT.ABACUS/STRU_ION_D` 文件中。

### 3.2.5 检查收敛

查看输出日志，确认优化收敛：

```bash
grep "TOTAL-PRESSURE" OUT.ABACUS/running_cell-relax.log | tail -5
grep "LARGEST GRAD" OUT.ABACUS/running_cell-relax.log | tail -5
```

收敛判据：
- `LARGEST GRAD < 0.01 eV/Å`（力收敛）
- `TOTAL-PRESSURE < 0.5 kBar`（应力收敛）

优化后的晶格常数约为 **4.04 Å**，与实验值 4.05 Å 接近（PBE 泛函通常略微低估）。

## 3.3 步骤 2：应变构型计算（relax）

### 3.3.1 准备优化后的结构

将优化后的结构复制到新的工作目录：

```bash
mkdir elastic_calc
cd elastic_calc
cp ../OUT.ABACUS/STRU_ION_D ./STRU
```

### 3.3.2 修改 INPUT 文件

将 `calculation` 改为 `relax`（固定晶格，允许原子弛豫）：

```
INPUT_PARAMETERS

# System variables
calculation         relax
symmetry            1
precision           double

# Plane wave related variables
ecutwfc             70.0
pw_diag_nmax        20
pw_diag_ndim        2

# Electronic structure
basis_type          pw
ks_solver           dav_subspace
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-08

# Geometry relaxation
relax_method        cg
relax_nmax          60
cal_force           1
force_thr_ev        0.01
cal_stress          1
stress_thr          0.5
fixed_axes          None
```

**关键变化：**
- `calculation = relax`（不再是 cell-relax）
- 移除 `kspacing`，使用固定 k 点网格

### 3.3.3 STRU 文件（优化后）

优化后的 STRU 文件（`STRU_ION_D`）内容示例：

```
ATOMIC_SPECIES
Al 26.981500 Al.upf

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    4.04160318060     0.00000000000    -0.00000000000
    0.00000000000     4.04160318060    -0.00000000000
    0.00000000000    -0.00000000000     4.04160318060

ATOMIC_POSITIONS
Direct

Al
0.000000
4
    0.00000000000     0.00000000000     0.00000000000 1 1 1
    0.00000000000     0.50000000000     0.50000000000 1 1 1
    0.50000000000     0.00000000000     0.50000000000 1 1 1
    0.50000000000     0.50000000000     0.00000000000 1 1 1
```

**变化说明：**
- 晶格常数：从 4.05 Å 优化为 4.04160 Å
- 坐标系：从 Cartesian 变为 Direct（分数坐标）
- 原子位置：优化后的分数坐标

### 3.3.4 KPT 文件

```
K_POINTS
0
Gamma
6 6 6 0 0 0
```

保持与步骤 1 相同的 k 点设置。

### 3.3.5 生成应变构型

使用 pymatgen 脚本生成 24 个应变构型。这部分将在第四章详细介绍，这里先说明原理：

**应变构型生成原理：**
1. 读取优化后的 STRU 文件
2. 对晶格矢量施加 6 种应变模式
3. 每种应变施加 4 种大小：±0.5%, ±1.0%
4. 生成 24 个 task 文件夹（task.000 到 task.023）

**每个 task 文件夹包含：**
- INPUT：与上述相同
- KPT：与上述相同
- STRU：应变后的结构
- strain.json：记录应变大小

### 3.3.6 批量计算应力

对 24 个应变构型分别运行 ABACUS 计算：

```bash
for i in {000..023}; do
    cd task.$i
    OMP_NUM_THREADS=1 mpirun -np 8 abacus
    cd ..
done
```

每个计算完成后，应力数据保存在 `OUT.ABACUS/running_relax.log` 中。

### 3.3.7 提取应力数据

从每个 task 的输出中提取应力张量：

```bash
grep "TOTAL-STRESS" task.*/OUT.ABACUS/running_relax.log
```

应力输出格式示例：

```
TOTAL-STRESS (KBAR):
-0.12   0.00   0.00
 0.00  -0.12   0.00
 0.00   0.00  -0.12
```

单位为 kBar，需要转换为 GPa（1 kBar = 0.1 GPa）。

## 3.4 结果解读

### 3.4.1 弹性常数矩阵

通过 pymatgen 拟合得到的弹性常数矩阵（单位：GPa）：

```
C_11  C_12  C_12   0     0     0
C_12  C_11  C_12   0     0     0
C_12  C_12  C_11   0     0     0
 0     0     0    C_44   0     0
 0     0     0     0    C_44   0
 0     0     0     0     0    C_44
```

典型计算结果（Al，PBE 泛函）：
- C_11 ≈ 114 GPa
- C_12 ≈ 62 GPa
- C_44 ≈ 32 GPa

### 3.4.2 弹性模量

从弹性常数计算得到的弹性模量：

| 弹性模量 | 计算公式 | 典型值（Al） |
|---------|---------|-------------|
| 体弹模量 B | (C_11 + 2C_12) / 3 | ~79 GPa |
| 剪切模量 G | (C_11 - C_12 + 3C_44) / 5 | ~28 GPa |
| 杨氏模量 E | 9BG / (3B + G) | ~72 GPa |
| 泊松比 ν | (3B - 2G) / (6B + 2G) | ~0.34 |

### 3.4.3 与实验值对比

| 物理量 | 计算值（PBE） | 实验值 | 偏差 |
|--------|--------------|--------|------|
| C_11 | ~114 GPa | 108 GPa | +5.6% |
| C_12 | ~62 GPa | 62 GPa | 0% |
| C_44 | ~32 GPa | 28 GPa | +14% |
| 体弹模量 B | ~79 GPa | 76 GPa | +3.9% |
| 杨氏模量 E | ~72 GPa | 70 GPa | +2.9% |

PBE 泛函对 Al 的弹性常数预测较为准确，偏差在 15% 以内。

# 第四章：使用 pymatgen 自动化计算

手动生成 24 个应变构型并逐个计算比较繁琐。pymatgen 提供了自动化工具，可以简化整个流程。

## 4.1 pymatgen 简介

### 什么是 pymatgen

pymatgen（Python Materials Genomics）是一个用于材料科学计算的 Python 库，由加州大学圣地亚哥分校开发。

**主要功能：**
- 晶体结构操作和分析
- 与 Materials Project 数据库交互
- 高通量计算工作流
- 弹性常数计算

**官方网站：** https://pymatgen.org/

### 安装 pymatgen

```bash
pip install pymatgen dpdata monty numpy
```

**依赖库说明：**
- `pymatgen`：核心库
- `dpdata`：数据格式转换
- `monty`：工具函数
- `numpy`：数值计算

## 4.2 生成应变构型

### 4.2.1 准备工作

在工作目录下准备：
- `relax/`：包含优化后的 STRU 文件（STRU_ION_D）
- `INPUT`：relax 计算的 INPUT 文件
- `KPT`：k 点文件
- 赝势文件：Al.upf

### 4.2.2 使用 gene_dfm.py 脚本

执行以下命令生成应变构型：

```bash
python gene_dfm.py abacus
```

**脚本功能：**
1. 读取 `relax/STRU_ION_D` 文件
2. 生成 6 种应变模式 × 4 种应变大小 = 24 个构型
3. 创建 `task.000` 到 `task.023` 文件夹
4. 每个文件夹包含：INPUT、KPT、STRU、strain.json

**应变大小设置：**

脚本中默认的应变大小（可修改）：

```python
norm_strains = [-0.010, -0.005, 0.005, 0.010]
shear_strains = [-0.010, -0.005, 0.005, 0.010]
```

### 4.2.3 检查生成的文件

进入任意 task 文件夹查看：

```bash
cd task.000
ls
```

应该看到：INPUT、KPT、STRU、strain.json

查看 strain.json 内容：

```bash
cat strain.json
```

示例输出：

```json
{"strain": [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]}
```

表示在 x 方向施加 1% 的正应变。

## 4.3 批量计算应力

### 4.3.1 使用 shell 脚本批量提交

创建 `run_all.sh` 脚本：

```bash
#!/bin/bash
for i in {000..023}; do
    cd task.$i
    echo "Running task.$i"
    OMP_NUM_THREADS=1 mpirun -np 8 abacus > log 2>&1
    cd ..
done
```

运行脚本：

```bash
bash run_all.sh
```

### 4.3.2 检查计算完成

确认所有 task 都已完成：

```bash
ls task.*/OUT.ABACUS/running_relax.log | wc -l
```

应该输出 24。

## 4.4 计算弹性常数

### 4.4.1 使用 compute_dfm.py 脚本

所有计算完成后，执行：

```bash
python compute_dfm.py abacus
```

**脚本功能：**
1. 读取所有 task 的应力数据
2. 读取对应的应变数据（strain.json）
3. 对每种应变状态进行线性拟合
4. 计算弹性常数矩阵和弹性模量

### 4.4.2 输出结果

屏幕输出示例：

```
# Elastic Constants in GPa
114.23  61.85  61.85   0.00   0.00   0.00
 61.85 114.23  61.85   0.00   0.00   0.00
 61.85  61.85 114.23   0.00   0.00   0.00
  0.00   0.00   0.00  31.67   0.00   0.00
  0.00   0.00   0.00   0.00  31.67   0.00
  0.00   0.00   0.00   0.00   0.00  31.67
# Bulk Modulus BV = 79.31 GPa
# Shear Modulus GV = 28.12 GPa
# Youngs Modulus EV = 72.45 GPa
# Poission Ratio uV = 0.34
```

**结果说明：**
- 第 2-7 行：弹性常数矩阵（6×6）
- Bulk Modulus：体弹模量
- Shear Modulus：剪切模量
- Youngs Modulus：杨氏模量
- Poission Ratio：泊松比

### 4.4.3 保存详细结果

精度更高的结果保存在 `elastic.json` 文件中：

```bash
cat elastic.json
```

包含完整的弹性常数矩阵和所有弹性模量。

## 4.5 完整流程总结

使用 pymatgen 计算弹性常数的完整流程：

**步骤 1：结构优化**
```bash
# 准备 INPUT（cell-relax）、STRU、KPT
mkdir relax
cd relax
# 运行 ABACUS
OMP_NUM_THREADS=1 mpirun -np 8 abacus
cd ..
```

**步骤 2：生成应变构型**
```bash
python gene_dfm.py abacus
```

**步骤 3：批量计算应力**
```bash
bash run_all.sh
```

**步骤 4：计算弹性常数**
```bash
python compute_dfm.py abacus
```

**输出文件：**
- 屏幕输出：弹性常数矩阵和弹性模量
- `elastic.json`：详细结果（JSON 格式）

# 第五章：使用 abacustest 简化流程

abacustest 是 ABACUS 的前后处理工具，提供了更简洁的弹性常数计算接口。

## 5.1 abacustest 简介

abacustest 是专为 ABACUS 设计的高通量计算工具，集成了多种常用功能。

**主要功能：**
- 输入文件准备（从 CIF 生成 STRU）
- 批量任务提交
- 结果提取和可视化
- 弹性常数计算
- DOS/PDOS 绘图
- 能带绘图

**安装：**
```bash
pip install abacustest
```

**验证安装：**
```bash
abacustest -h
```

## 5.2 使用 abacustest 计算弹性常数

### 5.2.1 准备输入文件

abacustest 需要一个已优化的结构作为起点。准备以下文件：
- INPUT（relax 计算）
- STRU（优化后的结构）
- KPT
- 赝势文件

### 5.2.2 执行弹性计算

abacustest 提供了一键计算弹性常数的命令（具体命令请参考最新文档）：

```bash
abacustest model elastic -j ./
```

**功能：**
1. 自动生成 24 个应变构型
2. 批量提交 ABACUS 计算
3. 提取应力数据
4. 拟合弹性常数
5. 输出结果

> **注意：** abacustest 的命令和参数可能随版本更新而变化，使用前请执行 `abacustest model elastic -h` 查看最新帮助信息。

### 5.2.3 输出结果

abacustest 会自动输出弹性常数和弹性模量，格式与 pymatgen 类似。

## 5.3 pymatgen vs abacustest

| 对比项 | pymatgen | abacustest |
|--------|----------|-----------|
| 安装 | 需要多个依赖库 | 单个包安装 |
| 使用难度 | 需要理解脚本 | 命令行直接使用 |
| 灵活性 | 高，可自定义脚本 | 中等 |
| 适用场景 | 深度定制 | 快速计算 |
| 文档 | pymatgen 官方文档 | abacustest GitHub |

**建议：**
- 初学者：使用 abacustest，简单快速
- 高级用户：使用 pymatgen，灵活可控
- 高通量计算：两者均可，根据工作流选择

# 附录

## A.1 参考资料

**ABACUS 官方资源：**
- ABACUS 官方文档：https://abacus.deepmodeling.com/en/latest/
- ABACUS GitHub：https://github.com/deepmodeling/abacus-develop
- ABACUS 中文教程：https://mcresearch.github.io/abacus-user-guide/

**弹性常数计算：**
- ABACUS+pymatgen 计算弹性常数：https://mcresearch.github.io/abacus-user-guide/abacus-elastic.html
- Materials Project 弹性常数方法：https://docs.materialsproject.org/methodology/materials-methodology/elasticity

**工具文档：**
- pymatgen 官方文档：https://pymatgen.org/
- abacustest GitHub：https://github.com/pxlxingliang/abacus-test

## A.2 常见问题

**Q1：为什么我的弹性常数与实验值差异较大？**

可能原因：
- 截断能或 k 点不够收敛（金属尤其敏感）
- 泛函选择（PBE 对某些材料误差较大）
- 结构未充分优化
- 应变大小不合适（超出线性范围）

**Q2：金属体系计算弹性常数需要注意什么？**

- k 点密度要足够高（建议 8×8×8 或更密）
- 使用 Methfessel-Paxton 展宽（`smearing_method mp`）
- 展宽宽度适中（0.01-0.02 Ry）
- SCF 收敛阈值要严格（1e-8 或更小）

**Q3：如何判断计算结果是否可靠？**

- 检查应力-应变曲线的线性度（R² > 0.99）
- 验证弹性常数的力学稳定性判据
- 与实验值或文献值对比
- 做收敛性测试（ecutwfc、k 点）

**Q4：计算很慢怎么办？**

- 使用更多并行核心
- 减小 ecutwfc（在收敛范围内）
- 使用 LCAO 基组（大体系）
- 在云平台（如 Bohrium）上计算

## A.3 进阶学习

完成本篇后，建议学习：

1. **收敛性测试**
   - ecutwfc 对弹性常数的影响
   - k 点密度对金属体系的影响

2. **其他晶系的弹性常数**
   - 六方晶系（5 个独立常数）
   - 正交晶系（9 个独立常数）

3. **高级功能**
   - 温度相关的弹性性质
   - 压力下的弹性常数
   - 各向异性分析

4. **后续教程**
   - 第三篇：能带结构计算
   - 第四篇：态密度（DOS/PDOS）计算
