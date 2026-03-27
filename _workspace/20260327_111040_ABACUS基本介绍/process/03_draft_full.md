# 前言

本教程是 ABACUS 系列培训的第一篇，面向刚接触第一性原理计算的初学者。读完本篇，你将能够独立准备 ABACUS 的输入文件、运行一次基本计算、读懂输出结果。

## 适用读者

- 希望使用 ABACUS 开展材料计算研究的科研人员和学生
- 有一定 DFT 概念基础，但没有实际操作经验的用户
- 从其他软件（如 VASP、Quantum ESPRESSO）迁移到 ABACUS 的用户

## 前置知识

- 基本的 Linux 命令行操作（文件创建、目录切换、命令执行）
- 密度泛函理论（DFT）的基本概念（平面波基组、赝势、自洽迭代）
- 不要求编程能力

## 本篇结构

| 章节 | 内容 |
|------|------|
| 第一章 | ABACUS 软件背景与功能概览 |
| 第二章 | INPUT 文件：计算参数详解 |
| 第三章 | KPT 文件：布里渊区采样 |
| 第四章 | STRU 文件：结构信息 |
| 第五章 | 结构优化实战（relax/cell-relax）|
| 第六章 | 输出文件解读 |
| 第七章 | abacustest：输入准备与结果提取 |

## 后续篇章

本系列共四篇：
1. **ABACUS 基本介绍**（本篇）
2. 力学性质计算：弹性模量
3. 能带结构计算
4. 态密度（DOS/PDOS）计算

建议按顺序阅读，每篇均以前篇为基础。
# 第一章：ABACUS 软件简介

## 1.1 软件背景

ABACUS（Atomic-orbital Based Ab-initio Computation at UStc），中文名**原子算筹**，是国内自主研发的开源第一性原理材料计算软件。软件基于密度泛函理论（Density Functional Theory，DFT），采用模守恒赝势和周期性边界条件，支持平面波（Plane Wave，PW）和数值原子轨道（Numerical Atomic Orbital，NAO/LCAO）两种基组。

ABACUS 托管于 DeepModeling 开源社区，代码完全公开，任何人均可参与贡献：

- **GitHub 仓库：** https://github.com/deepmodeling/abacus-develop
- **官方文档：** https://abacus.deepmodeling.com/en/latest/
- **中文教程：** https://mcresearch.github.io/abacus-user-guide/

软件免费开源，可在 Linux 集群、云计算平台（如 Bohrium）上运行。

## 1.2 核心功能

ABACUS 目前支持以下主要计算类型：

| 计算类型 | 关键字 | 说明 |
|----------|--------|------|
| 电子自洽迭代 | `scf` | 求解电子基态，获得总能量、电荷密度 |
| 非自洽计算 | `nscf` | 基于已有电荷密度计算能带/DOS |
| 结构优化 | `relax` | 固定晶胞，优化原子位置 |
| 晶胞优化 | `cell-relax` | 同时优化晶胞参数和原子位置 |
| 分子动力学 | `md` | 模拟原子随时间的运动轨迹 |

除基本功能外，ABACUS 还支持：
- DFT+U（处理强关联体系）
- 范德华修正（vdW correction）
- 隐式溶剂模型
- 杂化泛函（HSE06 等）
- 自旋极化与非共线磁矩计算
- 与 DeePMD-kit 联用的机器学习分子动力学

## 1.3 两种基组：PW 与 LCAO

ABACUS 同时支持两种基组，各有适用场景：

| 对比项 | 平面波（PW） | 数值原子轨道（LCAO） |
|--------|-------------|----------------------|
| 控制参数 | `basis_type pw` | `basis_type lcao` |
| 额外文件 | 仅赝势（.upf） | 赝势（.upf）+ 轨道（.orb）|
| 精度控制 | ecutwfc 截断能 | 轨道文件精度 |
| 计算速度 | 中等 | 通常更快（尤其大体系）|
| 适用场景 | 通用，精度基准 | 大体系、高通量计算 |

初学者建议先用 PW 基组入手，熟悉后再尝试 LCAO。

## 1.4 ABACUS 计算文件体系

运行一个 ABACUS 计算，至少需要准备以下三个文件：

```
计算目录/
├── INPUT        # 计算参数控制（必需）
├── KPT          # 布里渊区 k 点采样（必需）
├── STRU         # 晶体结构信息（必需）
├── *.upf        # 赝势文件（必需）
└── *.orb        # 数值轨道文件（LCAO 基组才需要）
```

计算完成后，输出文件保存在：

```
计算目录/
├── OUT.ABACUS/           # 输出文件夹（ABACUS 为默认后缀）
│   ├── INPUT             # 包含所有参数（含默认值）的完整 INPUT
│   ├── running_scf.log   # 详细运行日志
│   ├── KPT.info          # k 点信息
│   └── ...               # 其他输出文件
└── time.json             # 各模块计时信息
```

下面三章分别介绍这三个核心输入文件。
# 第二章：INPUT 文件详解

INPUT 文件是 ABACUS 计算的核心控制文件，决定了计算类型、使用的基组、精度要求以及输出内容。

## 2.1 文件格式与基本规则

INPUT 文件必须以 `INPUT_PARAMETERS` 作为第一行（有效内容），之前的内容会被忽略。参数格式为：

```
关键字   值
```

注意：**不使用等号**（区别于 VASP 的 INCAR）。

其他规则：
- 以 `#` 或 `/` 开头的行为注释，整行忽略
- 同一参数出现多次时，**取最后一次**的值（区别于 VASP 取第一次）
- 布尔值可用 `1`/`0`、`True`/`False`、`T`/`F`，大小写均可
- 文件名固定为 `INPUT`，不可更改

一个最简 INPUT 示例（Si 的 SCF 计算）：

```
INPUT_PARAMETERS
ntype           1
pseudo_dir      ./
ecutwfc         50
basis_type      pw
calculation     scf
```

## 2.2 通用参数

这些参数在几乎所有计算中都需要设置：

### suffix

```
suffix          ABACUS
```

输出文件夹的后缀名。计算完成后会生成 `OUT.ABACUS/` 文件夹。建议设为体系名称，如 `suffix Si`，则输出到 `OUT.Si/`。

### ntype

```
ntype           2
```

体系中元素的种类数。有多少种元素就填几。**必须与 STRU 文件中的元素种数一致。**

### pseudo_dir / orbital_dir

```
pseudo_dir      ./
orbital_dir     ./
```

赝势文件和轨道文件所在的目录。可以使用相对路径或绝对路径。`./` 表示当前目录。使用 LCAO 基组时，`orbital_dir` 必须设置。

### calculation

```
calculation     scf
```

最重要的参数之一，决定计算类型：

| 值 | 含义 |
|----|------|
| `scf` | 自洽电子结构计算，获取总能量和基态电荷密度 |
| `nscf` | 非自洽计算，用于能带/DOS（需已有电荷密度）|
| `relax` | 固定晶胞，优化原子位置 |
| `cell-relax` | 同时优化晶胞参数和原子位置 |
| `md` | 分子动力学模拟 |

### basis_type

```
basis_type      pw
```

基组类型：
- `pw`：平面波，通用性强，精度基准
- `lcao`：数值原子轨道，速度更快，需额外提供轨道文件

## 2.3 精度控制参数

### ecutwfc

```
ecutwfc         50
```

平面波截断能，单位 **Rydberg（Ry）**。控制平面波基组的大小——截断能越大，基组越完备，计算越精确，但计算量也越大。

- 默认值：50 Ry
- 典型范围：50–100 Ry（轻元素偏低，重元素或含 d/f 轨道偏高）
- **即使使用 LCAO 基组，也需要设置** `ecutwfc`（局部赝势和力的计算仍用平面波）

### scf_thr

```
scf_thr         1e-6
```

SCF 自洽迭代的收敛阈值，单位 **Ry**。当相邻两步迭代的电荷密度差小于该值时，认为收敛。

- 常规计算：`1e-6`（默认）
- 高精度计算（如力、应力）：`1e-7` 或更小
- 快速预判：`1e-4`（不推荐用于最终结果）

### scf_nmax

```
scf_nmax        100
```

SCF 迭代的最大步数，默认 100。若体系较难收敛（如磁性金属），可适当增大至 200–300。

## 2.4 常用功能参数

### smearing_method / smearing_sigma

展宽（smearing）方法，用于处理费米面附近的电子占据，改善 SCF 收敛性：

```
smearing_method     gauss
smearing_sigma      0.01
```

`smearing_method` 可选值：

| 值 | 适用场景 |
|----|----------|
| `gauss` | 绝缘体、半导体（默认）|
| `mp` | 金属（推荐）|
| `mp2` | 金属（高阶 Methfessel-Paxton）|
| `fd` | 有限温度模拟（Fermi-Dirac）|
| `fixed` | 绝缘体，固定占据数 |

`smearing_sigma` 为展宽宽度，单位 Ry，默认 0.015 Ry：
- 金属：0.005–0.02 Ry
- 绝缘体：0.01 Ry（实际影响不大）

### symmetry

```
symmetry        1
```

是否使用晶体对称性：
- `1`：打开对称性（加速计算，默认）
- `0`：关闭晶体对称性，保留时间反演对称性
- `-1`：关闭所有对称性

做能带计算（nscf）时通常设为 `0`。

### cal_force / cal_stress

```
cal_force       1
cal_stress      1
```

是否计算原子受力和晶胞应力。做结构优化时必须开启（relax/cell-relax 会自动开启）。

### nspin

```
nspin           1
```

自旋极化设置：
- `1`：无自旋极化（默认）
- `2`：共线自旋极化（铁磁/反铁磁体系）
- `4`：非共线磁矩（需设置 `noncolin 1`）

## 2.5 结构优化专用参数

使用 `calculation relax` 或 `calculation cell-relax` 时需要额外设置：

### force_thr_ev

```
force_thr_ev    0.01
```

原子受力的收敛阈值，单位 **eV/Å**。所有原子中的最大受力小于该值时认为力收敛。

- 常规优化：0.01–0.05 eV/Å
- 高精度：0.005 eV/Å 或更小

### stress_thr

```
stress_thr      0.5
```

晶胞应力的收敛阈值，单位 **kBar**（`cell-relax` 专用）。

- 常规优化：0.5–5 kBar
- 高精度：0.1 kBar

### relax_nmax

```
relax_nmax      100
```

结构优化的最大离子迭代步数，默认 50。复杂体系可增大至 100–200。

### out_stru

```
out_stru        1
```

是否在每个离子步后输出 STRU 文件，默认 0。设为 1 时，每步优化结果保存为 `OUT.suffix/STRU_ION_D`，便于追踪优化过程。

## 2.6 输出控制参数

### out_chg

```
out_chg         1
```

是否输出电荷密度文件（`SPIN1_CHG.cube` 等）。做能带计算时需先用 SCF 输出电荷密度，再用 NSCF 读入。

### init_chg

```
init_chg        atomic
```

初始电荷密度的来源：
- `atomic`：从原子电荷密度叠加（默认）
- `file`：从已有的 `SPIN1_CHG` 等文件读入（NSCF 时使用）

## 2.7 完整示例

### 示例1：Si 晶体 SCF 计算（PW 基组）

```
INPUT_PARAMETERS

# 通用参数
suffix          Si
ntype           1
pseudo_dir      ./
calculation     scf
basis_type      pw

# 精度控制
ecutwfc         50          # 单位：Ry
scf_thr         1e-6        # 单位：Ry
scf_nmax        100

# 展宽（Si 是半导体）
smearing_method gauss
smearing_sigma  0.01

# 对称性
symmetry        1
```

### 示例2：MgO 晶体结构优化（LCAO 基组）

```
INPUT_PARAMETERS

# 通用参数
suffix          MgO
ntype           2
pseudo_dir      ./
orbital_dir     ./
calculation     cell-relax
basis_type      lcao

# 精度控制
ecutwfc         100         # 单位：Ry
scf_thr         1e-6
scf_nmax        100

# 展宽（MgO 是绝缘体）
smearing_method gauss
smearing_sigma  0.01

# 对称性
symmetry        1

# 结构优化收敛判据
force_thr_ev    0.01        # 单位：eV/Å
stress_thr      0.5         # 单位：kBar
relax_nmax      100
out_stru        1
```

> **参数设置原则：** 参数不是越多越好。对于大多数场景，默认值已经是合理选择。只设置你明确知道含义的参数。完整的参数列表见官方文档：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html
# 第三章：KPT 文件详解

KPT 文件用于设置布里渊区的 k 点采样，类似于 VASP 的 KPOINTS 文件。

## 3.1 布里渊区采样基础

晶体具有平移周期性，在倒空间（k 空间）中对应一个布里渊区。求解 Kohn-Sham 方程时，需要在布里渊区内进行积分，实际计算中用有限个 k 点离散近似。

k 点采样的关键原则：
- k 点越密，计算越精确，但耗时越多
- 金属体系需要较密的 k 点（费米面附近状态变化快）
- 绝缘体/半导体对 k 点密度要求相对较低
- 扩胞后体积增大，倒空间收缩，k 点可相应减少

## 3.2 KPT 文件格式

标准的自动生成格式如下：

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

逐行说明：

**第1行：** `K_POINTS`，固定关键字，不可更改。

**第2行：** k 点总数。设为 `0` 时，表示自动生成 Monkhorst-Pack 网格，由第3、4行控制。

**第3行：** 采样方法：
- `Gamma`：以 Γ 点为中心的 Monkhorst-Pack 方法（最常用）
- `MP`：标准 Monkhorst-Pack 方法（网格不以 Γ 为中心）

**第4行：** `n1 n2 n3 s1 s2 s3`
- `n1 n2 n3`：沿三个倒格矢方向各细分的份数
- `s1 s2 s3`：网格的偏移量（通常取 `0 0 0`，不偏移）

## 3.3 k 点设置建议

### 根据晶系调整

| 晶系 | 建议 |
|------|------|
| 立方 | 各向同性，如 `4 4 4` |
| 六方/三方 | 面内与 c 轴方向不同，如 `6 6 4` |
| 正方/斜方 | 根据倒格矢长度比例调整 |
| 低对称 | 参考文献或做收敛性测试 |

### 根据体系类型调整

- **金属：** k 点适当加密，如 `8 8 8` 或更密
- **绝缘体/半导体：** `4 4 4` 到 `6 6 6` 通常足够
- **大体系（> 50 原子）：** 可以用较稀的 k 点，甚至只用 Γ 点

### gamma_only 模式

对于原子数较多的大体系，可以在 INPUT 中设置：

```
gamma_only      1
```

此时只使用 Γ 点计算，KPT 文件内容被忽略。好处是计算内存和速度显著提升（LCAO 基组尤为明显）。

## 3.4 完整示例

### 示例1：立方体系（如 Si、MgO）

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

### 示例2：六方体系（如 Mg、石墨）

```
K_POINTS
0
Gamma
6 6 4 0 0 0
```

### 示例3：大体系，只用 Γ 点

KPT 文件内容可以忽略，在 INPUT 中设置 `gamma_only 1` 即可。或者写：

```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

> **注意：** KPT 文件中不支持写注释（`#` 不会被识别为注释）。格式必须严格按照上述结构。
# 第四章：STRU 文件详解

STRU 文件包含晶体结构的全部几何信息，以及赝势和轨道文件的路径。类似于 VASP 的 POSCAR + POTCAR 合并体。

## 4.1 文件结构概览

STRU 文件由若干"块"组成，各块按固定顺序排列：

```
ATOMIC_SPECIES          # 元素种类与赝势
NUMERICAL_ORBITAL       # 轨道文件（LCAO 才需要）
LATTICE_CONSTANT        # 晶格常数（缩放因子）
LATTICE_VECTORS         # 晶格矢量
ATOMIC_POSITIONS        # 原子坐标
```

PW 基组不需要 `NUMERICAL_ORBITAL` 块，其余块相同。

## 4.2 ATOMIC_SPECIES 块

格式：

```
ATOMIC_SPECIES
元素名   原子质量   赝势文件名
```

示例（MgO，两种元素）：

```
ATOMIC_SPECIES
Mg   24.305   Mg_ONCV_PBE-1.0.upf
O    15.999   O_ONCV_PBE-1.0.upf
```

- **元素名**：与 `ATOMIC_POSITIONS` 块中一致
- **原子质量**：用于分子动力学，SCF/relax 计算中不影响结果
- **赝势文件名**：文件需放在 INPUT 中 `pseudo_dir` 指定的目录下

**赝势文件获取：** http://abacus.ustc.edu.cn/pseudo/list.htm
推荐使用 ONCV 模守恒赝势（PBE 泛函），命名规范为 `元素_ONCV_PBE-1.0.upf`。

## 4.3 NUMERICAL_ORBITAL 块（LCAO 专用）

使用 `basis_type lcao` 时，需要额外指定每个元素的轨道文件：

```
NUMERICAL_ORBITAL
Mg_gga_8au_100Ry_4s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb
```

- 每个元素对应一行，顺序与 `ATOMIC_SPECIES` 中一致
- 轨道文件名的含义：`元素_gga_截断半径au_截断能Ry_轨道配置.orb`
  - `8au`：轨道截断半径为 8 a.u.（原子单位）
  - `100Ry`：对应截断能 100 Ry
  - `4s2p1d`：轨道配置（s、p、d 轨道的数量）
- 轨道文件需与所用赝势匹配（同一套 PBE 泛函）
- 文件放在 INPUT 中 `orbital_dir` 指定的目录下

**PW 基组时删去此块即可，其余不变。**

## 4.4 LATTICE_CONSTANT 块

```
LATTICE_CONSTANT
1.8897261258369282
```

这是一个缩放因子。后续的 `LATTICE_VECTORS` 中的数值乘以此常数，得到实际的晶格矢量长度（单位：Bohr）。

常用的对应关系：

| LATTICE_CONSTANT 值 | LATTICE_VECTORS 单位 |
|---------------------|----------------------|
| `1.8897261258369282` | Å（埃）|
| `1.0` | Bohr（原子单位）|

最常见的做法是将 LATTICE_CONSTANT 设为 `1.8897261258369282`，然后 LATTICE_VECTORS 中直接填写以 **Å** 为单位的晶格矢量。

## 4.5 LATTICE_VECTORS 块

```
LATTICE_VECTORS
4.27957   0.00000   0.00000
0.00000   4.27957   0.00000
0.00000   0.00000   4.27957
```

三行，每行三个数，构成 3×3 矩阵，描述三个晶格基矢的方向和（相对）长度。

- 实际晶格矢量 = 矩阵数值 × LATTICE_CONSTANT（单位 Bohr）
- 支持非正交晶格（如六方、三斜）
- 若 LATTICE_CONSTANT 设为 `1.8897261258369282`，则直接填 Å 数值

六方晶格示例（MgB₂，a = 3.083 Å，c = 3.524 Å）：

```
LATTICE_VECTORS
3.083000   0.000000   0.000000
-1.541500  2.669500   0.000000
0.000000   0.000000   3.524000
```

## 4.6 ATOMIC_POSITIONS 块

```
ATOMIC_POSITIONS
Direct
元素名
初始磁矩
原子数
x  y  z  mx  my  mz
...
```

逐项说明：

**第1行：** 坐标系类型
- `Direct`：晶格坐标（分数坐标，0 到 1 之间）
- `Cartesian`：笛卡尔坐标，单位为 LATTICE_CONSTANT（Bohr 或 Å）

**元素名：** 与 ATOMIC_SPECIES 中一致。

**初始磁矩：** 该元素所有原子的初始磁矩值（Bohr-mag）。非磁性计算填 `0.0`。

**原子数：** 该元素的原子个数。

**坐标行：** 每个原子一行，格式为：
```
x  y  z  mx  my  mz
```
- `x y z`：坐标值
- `mx my mz`：各方向是否允许移动（`1` = 可移动，`0` = 固定），用于 relax 计算

## 4.7 完整示例

### 示例1：Si 晶体（PW 基组，Diamond 结构，8原子超胞）

```
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
5.43070   0.00000   0.00000
0.00000   5.43070   0.00000
0.00000   0.00000   5.43070

ATOMIC_POSITIONS
Direct
Si
0.0
8
0.000  0.000  0.000  1  1  1
0.250  0.250  0.250  1  1  1
0.500  0.500  0.000  1  1  1
0.750  0.750  0.250  1  1  1
0.500  0.000  0.500  1  1  1
0.750  0.250  0.750  1  1  1
0.000  0.500  0.500  1  1  1
0.250  0.750  0.750  1  1  1
```

### 示例2：MgO 晶体（LCAO 基组，岩盐结构，8原子超胞）

```
ATOMIC_SPECIES
Mg   24.305   Mg_ONCV_PBE-1.0.upf
O    15.999   O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Mg_gga_8au_100Ry_4s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
4.27957   0.00000   0.00000
0.00000   4.27957   0.00000
0.00000   0.00000   4.27957

ATOMIC_POSITIONS
Direct
Mg
0.0
4
0.0   0.0   0.0   0  0  0
0.0   0.5   0.5   0  0  0
0.5   0.0   0.5   0  0  0
0.5   0.5   0.0   0  0  0
O
0.0
4
0.5   0.0   0.0   0  0  0
0.5   0.5   0.5   0  0  0
0.0   0.0   0.5   0  0  0
0.0   0.5   0.0   0  0  0
```

注意：MgO 示例中原子坐标后的 `0 0 0` 表示结构优化时固定所有原子（也可改为 `1 1 1` 允许移动）。

> **获取晶体结构：** 可从 Materials Project（https://materialsproject.org）、ICSD 等数据库下载 CIF 文件，再通过工具（如 ASE、abacustest）转换为 STRU 格式。下一章的案例中会展示 abacustest 的转换用法。
# 第五章：结构优化实战

前几章分别介绍了三个输入文件。本章以一个具体例子，展示三个文件如何配合完成一次结构优化计算。

## 5.1 结构优化的物理意义

实验制备或数据库中的晶体结构并非总处于 DFT 意义下的能量极小值。结构优化（几何弛豫）通过调整原子位置和晶格参数，使体系达到理论上的稳定构型，是大多数材料计算工作的第一步。

**relax 与 cell-relax 的区别：**

| 类型 | 优化对象 | 需要额外设置 |
|------|----------|--------------|
| `relax` | 原子位置（晶格固定）| `force_thr_ev` |
| `cell-relax` | 原子位置 + 晶格矢量 | `force_thr_ev` + `stress_thr` |

通常建议先做 `cell-relax` 同时弛豫晶格和原子，再视情况做 `relax` 精细优化原子位置。

## 5.2 完整案例：Si 晶体 cell-relax

以金刚石结构硅（8 原子超胞）为例，演示一次完整的结构优化。

### INPUT 文件

```
INPUT_PARAMETERS

suffix          Si
ntype           1
pseudo_dir      ./
calculation     cell-relax
basis_type      pw

ecutwfc         50
scf_thr         1e-6
scf_nmax        100

smearing_method gauss
smearing_sigma  0.01

symmetry        1

force_thr_ev    0.01
stress_thr      0.5
relax_nmax      100
out_stru        1
```

### KPT 文件

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

### STRU 文件

```
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
5.50000   0.00000   0.00000
0.00000   5.50000   0.00000
0.00000   0.00000   5.50000

ATOMIC_POSITIONS
Direct
Si
0.0
8
0.000  0.000  0.000  1  1  1
0.250  0.250  0.250  1  1  1
0.500  0.500  0.000  1  1  1
0.750  0.750  0.250  1  1  1
0.500  0.000  0.500  1  1  1
0.750  0.250  0.750  1  1  1
0.000  0.500  0.500  1  1  1
0.250  0.750  0.750  1  1  1
```

注意：初始晶格常数故意设为 5.50 Å（实验值约 5.43 Å），目的是观察 cell-relax 对晶格的优化效果。

### 运行计算

```bash
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

## 5.3 判断优化是否收敛

每完成一步离子迭代，输出中会显示：

```
TOTAL-PRESSURE: -2.070e+00 KBAR
ETOT DIFF (eV)       : -1.234e-03
LARGEST GRAD (eV/A)  : 8.521e-03
```

**收敛判据（cell-relax 需同时满足）：**

| 量 | 含义 | 收敛条件 |
|----|------|----------|
| `LARGEST GRAD (eV/A)` | 所有原子中最大受力 | < `force_thr_ev`（本例 0.01 eV/Å）|
| `TOTAL-PRESSURE (KBAR)` | 晶胞总压强（应力迹的三分之一）| < `stress_thr`（本例 0.5 kBar）|

两个条件同时满足时，计算收敛并自动停止。`relax` 只需满足力的判据。

## 5.4 查看优化结果

**优化轨迹：** 设置 `out_stru 1` 后，每步优化的结构保存在 `OUT.Si/STRU_ION_D`，可用于追踪优化过程。

**最终结构：** 优化收敛后，最终的晶格参数和原子位置会出现在 `running_cell-relax.log`（或屏幕输出）中：

```
 Cell Parameters:  (Angstrom)
  a = 5.43020   b = 5.43020   c = 5.43020
  alpha = 90.000   beta = 90.000   gamma = 90.000
```

优化后的晶格常数（约 5.43 Å）与实验值吻合。最终结构也以 STRU 格式保存在 `OUT.Si/` 下，可直接用于后续计算。
# 第六章：输出文件解读

ABACUS 计算完成后，会在当前目录生成输出文件夹和日志文件。理解这些输出是判断计算是否成功、提取物理量的关键。

## 6.1 输出文件结构

计算完成后，目录结构如下：

```
计算目录/
├── INPUT
├── KPT
├── STRU
├── time.json              # 各模块计时信息（JSON 格式）
└── OUT.suffix/            # 主要输出文件夹（suffix 对应 INPUT 中设置）
    ├── INPUT              # 包含所有参数（含默认值）
    ├── running_scf.log    # 详细运行日志（与屏幕输出基本相同）
    ├── kpoints            # k 点信息
    ├── STRU_READIN_ADJUST.cif  # 读入的结构（CIF 格式）
    └── ...                # 其他输出文件
```

**time.json：** 记录各模块耗时，用于性能分析。

**OUT.suffix/INPUT：** 包含所有实际使用的参数，包括用户未设置的默认值，便于确认计算设置。

**OUT.suffix/running_scf.log：** 最重要的输出文件，包含完整的计算过程和结果。

## 6.2 屏幕输出解读

运行 ABACUS 时，屏幕会实时显示计算进度。以下逐部分解读（以 SCF 计算为例）：

### 6.2.1 软件信息

```
                              ABACUS v3.10.1

               Atomic-orbital Based Ab-initio Computation at UStc

               Website: http://abacus.ustc.edu.cn/
               Documentation: https://abacus.deepmodeling.com/
               Repository: https://github.com/deepmodeling/abacus-develop
               Commit: 1234abcd (2024-01-15)

 Fri Jan 15 10:30:00 2024
 MAKE THE DIR         : OUT.Si/
```

显示软件版本、官网、文档地址、代码提交版本和计算开始时间。

### 6.2.2 参数回显

```
 READING GENERAL INFORMATION
                           global_out_dir : OUT.Si/
                           global_in_card : INPUT
                               pseudo_dir : ./
                              pseudo_type : auto
                                    ntype : 1
...
```

回显读入的参数，用于确认设置是否正确。若有警告（如赝势价电子数异常），会在此处显示。

### 6.2.3 SCF 迭代过程

```
 DONE(0.123 SEC) : INIT SCF
 ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)
 GE1    -1.234567e+02  0.000000e+00   5.432e-02  1.234e+00
 GE2    -1.235012e+02  -4.450000e-02  2.123e-02  1.198e+00
 GE3    -1.235089e+02  -7.700000e-03  8.765e-03  1.201e+00
 GE4    -1.235102e+02  -1.300000e-03  3.210e-03  1.195e+00
 GE5    -1.235105e+02  -3.000000e-04  1.234e-03  1.203e+00
 GE6    -1.235106e+02  -1.000000e-04  4.567e-04  1.199e+00
 GE7    -1.235106e+02  -1.000000e-05  1.234e-04  1.197e+00
 GE8    -1.235106e+02  -1.000000e-06  3.456e-05  1.201e+00
 GE9    -1.235106e+02  -1.000000e-07  8.765e-06  1.195e+00
 GE10   -1.235106e+02  -1.000000e-08  2.345e-06  1.198e+00
 GE11   -1.235106e+02  -1.000000e-09  6.789e-07  1.200e+00
```

各列含义：

| 列 | 含义 |
|----|------|
| ITER | 迭代步数（GE = Geometry Electronic）|
| ETOT(eV) | 当前步的总能量（单位：eV）|
| EDIFF(eV) | 与上一步的能量差（单位：eV）|
| DRHO | 电荷密度变化（无量纲）|
| TIME(s) | 本步耗时（秒）|

**收敛判断：** 当 `DRHO < scf_thr` 时，SCF 收敛。本例中 `scf_thr = 1e-6 Ry ≈ 1.36e-5 eV`，第11步 DRHO 已小于阈值，收敛成功。

### 6.2.4 计算结果汇总

```
 --------------------------------------------
 !FINAL_ETOT_IS -1.235106e+02 eV
 --------------------------------------------

 Fermi Energy (eV) = 5.4321
```

**总能量：** `!FINAL_ETOT_IS` 后的数值是体系的总能量（单位：eV），这是 DFT 计算最核心的结果。

**费米能级：** 金属体系的费米能级，或半导体/绝缘体的价带顶能量（需结合能带判断）。

### 6.2.5 各模块耗时统计

```
 --------------------------------------------------------------------------------
 CLASS_NAME         NAME                TIME(Sec)  CALLS   AVG(Sec) PER(%)
 --------------------------------------------------------------------------------
 Driver             driver_line         27.239     1       27.24    100.00
 PW_Basis           setup_struc_factor  0.123      1       0.12     0.45
 ORB_control        read_orb_first      1.234      1       1.23     4.53
 ...
 --------------------------------------------------------------------------------
```

显示各模块的耗时、调用次数、平均耗时和占比，用于性能分析和优化。

## 6.3 running_scf.log 文件

`OUT.suffix/running_scf.log` 与屏幕输出内容基本相同，但保存在文件中便于后续查阅。包含：
- 完整的参数设置
- SCF 迭代详细过程
- 最终能量和费米能级
- 各模块耗时统计

## 6.4 提取关键结果

### 提取总能量

```bash
grep "!FINAL_ETOT_IS" OUT.Si/running_scf.log
```

输出：
```
 !FINAL_ETOT_IS -1.235106e+02 eV
```

### 提取费米能级

```bash
grep "Fermi Energy" OUT.Si/running_scf.log
```

### 提取力（relax 计算）

```bash
grep "TOTAL-FORCE" OUT.Si/running_relax.log
```

## 6.5 常见警告与错误

### SCF 不收敛

```
 SCF NOT CONVERGED after 100 iterations
```

**可能原因：**
- 初始结构不合理（原子距离过近）
- smearing 参数不当（金属体系 sigma 过小）
- mixing 参数需要调整（降低 mixing_beta）
- ecutwfc 或 k 点不足

**解决方法：**
- 检查 STRU 文件，确保结构合理
- 金属体系使用 `smearing_method mp`，增大 `smearing_sigma`
- 降低 `mixing_beta`（如从 0.4 降至 0.1）
- 增大 `scf_nmax` 至 200

### 内存不足

```
 Error: Memory allocation failed
```

**解决方法：**
- 减小 ecutwfc
- 减少 k 点数量
- 使用更多节点/核心分配内存
- LCAO 基组比 PW 节省内存

### k 点设置不当

```
 Warning: k-point mesh too coarse
```

**解决方法：**
- 增加 k 点密度（如从 2×2×2 增至 4×4×4）
- 做 k 点收敛性测试

> **提示：** 遇到问题时，先查看 `OUT.suffix/warning.log`（如果存在），其中记录了所有警告信息。
# 第七章：使用 abacustest 准备输入与提取结果

abacustest 是 ABACUS 的前后处理工具，可以快速准备输入文件、批量提交任务、提取计算结果。

## 7.1 abacustest 简介

**功能定位：**
- 从 CIF/POSCAR 等格式生成 ABACUS 输入文件
- 批量设置计算参数
- 提取和可视化计算结果
- 高通量计算任务管理

**安装：**

```bash
# 方法1：通过 pip 安装
pip install abacustest

# 方法2：从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

**验证安装：**

```bash
abacustest -h
```

## 7.2 准备输入文件

### 7.2.1 基本用法

```bash
abacustest model inputs prepare -h
```

### 7.2.2 从 CIF 文件生成 STRU

假设有一个 Si.cif 文件：

```bash
abacustest model inputs prepare -f Si.cif --ftype cif
```

会在当前目录生成：
- `STRU`：结构文件
- `INPUT`：基本参数文件
- `KPT`：k 点文件

### 7.2.3 自定义参数

```bash
abacustest model inputs prepare \
  -f Si.cif \
  --ftype cif \
  --input "calculation=scf,basis_type=pw,ecutwfc=50" \
  --kpt "4 4 4 0 0 0"
```

参数说明：
- `--input`：设置 INPUT 参数（逗号分隔）
- `--kpt`：设置 k 点网格
- `--pp`：指定赝势文件路径
- `--orb`：指定轨道文件路径（LCAO）

## 7.3 提取计算结果

### 7.3.1 提取总能量

计算完成后，可以用 abacustest 提取关键结果：

```bash
abacustest model post -j ./
```

会自动识别 `OUT.suffix/` 文件夹，提取总能量、费米能级等信息。

### 7.3.2 提取结构信息

对于 relax/cell-relax 计算，提取优化后的结构：

```bash
abacustest model post -j ./ --extract-stru
```

会生成优化后的 STRU 文件。

## 7.4 实战示例：从 CIF 到计算结果

完整流程演示：

**步骤1：下载 CIF 文件**

从 Materials Project 下载 Si 的 CIF 文件（mp-149.cif）。

**步骤2：生成输入文件**

```bash
abacustest model inputs prepare \
  -f mp-149.cif \
  --ftype cif \
  --input "calculation=scf,basis_type=pw,ecutwfc=50,scf_thr=1e-6" \
  --kpt "4 4 4 0 0 0"
```

**步骤3：运行 ABACUS**

```bash
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

**步骤4：提取结果**

```bash
abacustest model post -j ./
```

输出示例：
```
Total Energy: -123.456 eV
Fermi Energy: 5.432 eV
```

## 7.5 进阶功能预告

abacustest 还支持更多高级功能，将在后续教程中介绍：
- DOS/PDOS 绘图（`abacustest model dos-pdos`）
- 能带绘图（`abacustest model band`）
- 弹性常数计算（`abacustest model elastic`）

> **提示：** abacustest 功能丰富，完整文档见 GitHub：https://github.com/pxlxingliang/abacus-test
# 附录

## A.1 参考资料

**官方资源：**
- ABACUS 官方文档：https://abacus.deepmodeling.com/en/latest/
- ABACUS GitHub 仓库：https://github.com/deepmodeling/abacus-develop
- ABACUS 中文教程：https://mcresearch.github.io/abacus-user-guide/
- 赝势和轨道文件下载：http://abacus.ustc.edu.cn/pseudo/list.htm

**社区资源：**
- DeepModeling 社区：https://github.com/deepmodeling
- Bohrium 云平台：https://bohrium.dp.tech/
- abacustest 工具：https://github.com/pxlxingliang/abacus-test

## A.2 常见问题

**Q1：如何选择赝势和轨道文件？**

推荐使用 ONCV 模守恒赝势（PBE 泛函）：
- 赝势：`元素_ONCV_PBE-1.0.upf`
- 轨道（LCAO）：`元素_gga_8au_100Ry_轨道配置.orb`

轨道文件需与赝势匹配（同一套泛函）。

**Q2：ecutwfc 如何设置？**

- 轻元素（H、C、N、O）：50–70 Ry
- 常见金属（Al、Fe、Cu）：60–80 Ry
- 重元素或含 d/f 轨道：80–100 Ry
- 建议做收敛性测试确定

**Q3：k 点如何设置？**

- 小体系（< 10 原子）：6×6×6 或更密
- 中等体系（10–50 原子）：4×4×4
- 大体系（> 50 原子）：2×2×2 或 gamma_only
- 金属体系需要更密的 k 点

**Q4：计算很慢怎么办？**

- 使用 LCAO 基组代替 PW
- 减小 ecutwfc（在精度允许范围内）
- 减少 k 点（大体系可用 gamma_only）
- 增加并行核心数

**Q5：SCF 不收敛怎么办？**

- 检查结构是否合理（原子距离不能过近）
- 金属体系使用 `smearing_method mp`
- 降低 `mixing_beta`（如 0.1–0.2）
- 增大 `scf_nmax`

## A.3 进阶学习方向

完成本篇后，建议按以下顺序学习：

1. **收敛性测试**（本系列未深入）
   - ecutwfc 收敛性测试
   - k 点收敛性测试

2. **力学性质计算**（系列第2篇）
   - 弹性常数计算
   - 体弹模量、剪切模量

3. **电子结构计算**（系列第3、4篇）
   - 能带结构
   - 态密度（DOS/PDOS）

4. **高级功能**
   - DFT+U（强关联体系）
   - 杂化泛函（HSE06）
   - 磁性材料计算
   - 分子动力学

5. **高通量计算**
   - 使用 abacustest 批量计算
   - 与 Bohrium 云平台联用
