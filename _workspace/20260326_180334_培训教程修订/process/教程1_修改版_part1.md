---
title: "ABACUS 输入输出体系：三文件入门与收敛性测试"
author: "AutoTutorial 3.0"
date: "2026-03-25"
topic: "ABACUS 输入输出体系"
task_type: "A"
has_case: false
word_count: ~4500
chapters: 3
---

# ABACUS 输入输出体系：三文件入门与收敛性测试

> 本文由 AutoTutorial 3.0 提供。

---

## 前言

本教程介绍 ABACUS 的三个核心输入文件：`STRU`、`KPT`、`INPUT`。掌握这三个文件的写法，是用 ABACUS 做任何第一性原理计算的前提。

**教程目标**

- 理解 STRU、KPT、INPUT 各自承担什么职责
- 学会为平面波（PW）和数值原子轨道（LCAO）两种基组写 STRU
- 掌握 k 点采样的两种格式：MP 网格（SCF 用）和路径 k 点（能带用）
- 理解 INPUT 中的关键参数分组，重点掌握 mixing 参数的物理含义
- 通过铝（Al）FCC 晶体的收敛性测试，学会确定 ecutwfc 和 k 点密度

**适用读者**

初次使用 ABACUS 的用户，或熟悉其他 DFT 软件（VASP、QE）需要迁移的用户。建议具备基本的 DFT 概念（自洽场计算、布里渊区采样）。

**前置知识**

- 晶体结构基本概念（晶格常数、分数坐标）
- 自洽场（SCF）计算的基本流程

**教程结构**

| 章节 | 主题 | 案例 |
|------|------|------|
| 第一章 | STRU 结构文件 | Al FCC（PW 和 LCAO 两种写法） |
| 第二章 | KPT k 点文件 | MP 网格 + 路径 k 点 + Al k 点收敛 |
| 第三章 | INPUT 计算参数 | Si 完整示例 + Al ecutwfc 收敛测试 |
| 附录 | 输出文件概览 | running_scf.log 关键信息 |

---

## 第一章 STRU：晶体结构文件
# ABACUS 软件背景

ABACUS（Atomic-orbital Based Ab-initio Computation at UStc，原子算筹）是中国科学技术大学开发的开源第一性原理计算软件。该软件基于密度泛函理论（DFT），支持平面波（PW）和数值原子轨道（LCAO）两种基组，可用于材料的电子结构、力学性质、光学性质等多种物理性质的计算。

## 软件特点

**开源与社区驱动**
- ABACUS 是完全开源的软件，代码托管在 GitHub（https://github.com/deepmodeling/abacus-develop）和 Gitee 镜像仓库
- 由 DeepModeling 开源社区维护，接受全球开发者贡献
- 文档完善，提供中英文双语支持（https://abacus.deepmodeling.com）

**双基组支持**
- 平面波（PW）基组：精度高，适合周期性体系，计算成本随体系大小快速增长
- 数值原子轨道（LCAO）基组：计算效率高，适合大体系（数百原子），精度可通过增加基函数控制

**丰富的功能**
- 电子结构计算：自洽场（SCF）、能带结构、态密度（DOS/PDOS）
- 结构优化：固定晶胞优化（relax）、变胞优化（cell-relax）
- 力学性质：弹性常数、声子谱
- 分子动力学：NVE、NVT、NPT 系综
- 高级功能：DFT+U、隐式溶剂模型、实时含时密度泛函理论（RT-TDDFT）等

**模守恒赝势**
- ABACUS 使用 ONCV（Optimized Norm-Conserving Vanderbilt）赝势，UPF 格式
- 推荐使用 SG15 ONCV 赝势库（http://abacus.ustc.edu.cn/pseudo/list.htm）

## 与其他软件的对比

| 特性 | ABACUS | VASP | Quantum ESPRESSO |
|------|--------|------|------------------|
| 开源 | ✅ 完全开源 | ❌ 商业软件 | ✅ 开源 |
| 基组 | PW + LCAO | PW + PAW | PW |
| 大体系效率 | ✅ LCAO 高效 | 中等 | 中等 |
| 学习曲线 | 中等 | 较陡 | 较陡 |
| 中文文档 | ✅ 完善 | 有限 | 有限 |

对于熟悉 VASP 的用户，ABACUS 的输入文件体系与 VASP 类似：
- STRU ≈ POSCAR（结构文件）
- INPUT ≈ INCAR（参数文件）
- KPT ≈ KPOINTS（k点文件）

## 适用场景

ABACUS 特别适合以下场景：
- 需要开源软件进行学术研究和教学
- 大体系计算（数百原子，使用 LCAO 基组）
- 需要中文文档和社区支持
- 探索新算法和方法开发（开源代码便于二次开发）

---
## 第一章 STRU：晶体结构文件

`STRU` 文件存储晶体的结构信息，包括元素种类、赝势/轨道文件名、晶格参数和原子坐标。ABACUS 中必须将这个文件命名为 `STRU`，不可更改。

### 1.1 文件整体结构

STRU 文件由若干以关键字开头的块组成，按顺序排列。对于平面波（PW）计算，必选块有四个；LCAO 计算需额外添加 `NUMERICAL_ORBITAL` 块：

```
# 注释行（可以用 # 开头）

ATOMIC_SPECIES          # 元素种类和赝势文件（必选）

NUMERICAL_ORBITAL       # 轨道文件（LCAO 专用，PW 不需要）

LATTICE_CONSTANT        # 晶格常数（必选）

LATTICE_VECTORS         # 晶格向量（必选）

ATOMIC_POSITIONS        # 原子坐标（必选）
```

各块的顺序不能颠倒。

### 1.2 各块详解

#### ATOMIC_SPECIES

```
ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf
```

每行格式：`元素名  原子量  赝势文件名`

- 元素名：需与 ATOMIC_POSITIONS 中的标签一致
- 原子量：用于分子动力学中的质量，SCF 计算结果与此值无关
- 赝势文件名：ABACUS 使用模守恒（ONCV）赝势，UPF 格式。文件放在 INPUT 中 `pseudo_dir` 指定的目录下。推荐使用 [SG15 ONCV 赝势库](http://abacus.ustc.edu.cn/pseudo/list.htm)

多元素体系每行写一个：

```
ATOMIC_SPECIES
Mg  24.305  Mg_ONCV_PBE-1.0.upf
O   15.999  O_ONCV_PBE-1.0.upf
```

#### NUMERICAL_ORBITAL（LCAO 专用）

```
NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb
```

每行一个轨道文件，顺序与 ATOMIC_SPECIES 中的元素对应。文件放在 INPUT 中 `orbital_dir` 指定目录下。PW 计算删去这整个块。

轨道文件命名格式：`元素_泛函_截断半径_能量截断_基组规模.orb`
- 例：`Al_gga_7au_100Ry_4s4p1d.orb` 表示 GGA 泛函、截断半径 7 Bohr、100 Ry 能量截断、4s4p1d 基函数

#### LATTICE_CONSTANT

```
LATTICE_CONSTANT
1.8897259886
```

晶格常数的单位因子（单位：Bohr）。LATTICE_VECTORS 中的数值乘以这个因子才是实际长度。

常用换算：
- `1.8897259886` Bohr = 1.0 Å（最常用写法：把长度写成 Å，LATTICE_CONSTANT 写 1.8897259886）
- 直接写 `1.0`：则 LATTICE_VECTORS 中的单位为 Bohr

#### LATTICE_VECTORS

```
LATTICE_VECTORS
4.03460  0.00000  0.00000
0.00000  4.03460  0.00000
0.00000  0.00000  4.03460
```

三行分别是晶格向量 **a₁、a₂、a₃**，单位为 LATTICE_CONSTANT。对于立方晶系，三个向量互相垂直且等长。

实际晶格常数 = LATTICE_VECTORS 中的值 × LATTICE_CONSTANT。上例中 Al FCC 晶格常数 = 4.03460 × 1.8897259886 Bohr = 4.03460 Å。

#### ATOMIC_POSITIONS

```
ATOMIC_POSITIONS
Direct            # 坐标类型

Al                # 元素名（与 ATOMIC_SPECIES 对应）
0                 # 初始磁矩（Bohr mag，非磁体系写 0）
4                 # 原子数目
0.0  0.0  0.0  0 0 0   # x y z  move_x move_y move_z
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

**坐标类型**（第一行）：
- `Direct`：分数坐标（推荐）。坐标值在 0~1 之间，以晶格向量为基
- `Cartesian`：直角坐标，单位为 LATTICE_CONSTANT
- `Cartesian_angstrom`：直角坐标，单位为 Å

**move_x move_y move_z**：控制原子在结构优化时是否可以移动（1=可移动，0=固定）。SCF 计算中这三个数没有实际作用，习惯写 `0 0 0` 或 `1 1 1` 均可。

多元素体系，各元素块依次排列（元素名→磁矩→原子数→坐标行）。

### 1.3 Al FCC 完整示例：PW 写法

铝（Al）面心立方（FCC）结构，原胞含 4 个原子，实验晶格常数 4.050 Å，计算优化值约 4.035 Å。

```
# Al FCC，PW 计算，4 原子原胞

ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897259886

LATTICE_VECTORS
4.03460  0.00000  0.00000
0.00000  4.03460  0.00000
0.00000  0.00000  4.03460

ATOMIC_POSITIONS
Direct
Al
0
4
0.0  0.0  0.0  0 0 0
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

注意：PW 写法不包含 `NUMERICAL_ORBITAL` 块。

### 1.4 Al FCC 完整示例：LCAO 写法

```
# Al FCC，LCAO 计算，4 原子原胞

ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.8897259886

LATTICE_VECTORS
4.03460  0.00000  0.00000
0.00000  4.03460  0.00000
0.00000  0.00000  4.03460

ATOMIC_POSITIONS
Direct
Al
0
4
0.0  0.0  0.0  0 0 0
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

与 PW 写法唯一的区别是添加了 `NUMERICAL_ORBITAL` 块。同时，INPUT 中需要将 `basis_type` 改为 `lcao`，并填写 `orbital_dir` 路径。

### 1.5 PW 与 LCAO 写法对比

| 项目 | PW | LCAO |
|------|----|----|
| NUMERICAL_ORBITAL 块 | 不需要 | 必须 |
| INPUT 中 basis_type | pw | lcao |
| INPUT 中 orbital_dir | 不需要 | 必须 |
| ecutwfc 意义 | 控制基组大小 | 控制辅助平面波格点密度 |
| 适用场景 | 收敛性测试、高精度 | 大体系、O(N) 线性标度 |

---

## 第二章 KPT：布里渊区 k 点

`KPT` 文件定义布里渊区的采样方案。k 点选得太少会引起能量误差；选得太多则浪费计算资源。文件同样必须命名为 `KPT`。

ABACUS 支持两种格式：**MP 网格**（用于 SCF/relax/MD）和**路径 k 点**（用于能带计算）。

### 2.1 MP 网格格式

Monkhorst-Pack 方法在布里渊区生成均匀网格，是自洽计算的标准选择。

```
K_POINTS          # 固定关键字，第一行
0                 # 0 表示自动生成网格
Gamma             # 采样方法：Gamma 或 MP
4  4  4  0 0 0    # 沿三个倒格矢方向的细分数 + 偏移
```

**第三行：Gamma 与 MP 的区别**
- `Gamma`：将 Γ 点（0,0,0）包含在网格内（推荐默认选择）
- `MP`：标准 Monkhorst-Pack，网格可能不包含 Γ 点

**第四行：6 个数字**
- 前三个：沿 **a₁\*、a₂\*、a₃\*** 方向的细分数。`4 4 4` 表示 4×4×4 的均匀网格
- 后三个：网格偏移量（通常写 `0 0 0`）

**另一种写法：kspacing 参数**

不使用 KPT 文件，直接在 INPUT 中设置 `kspacing 0.1`（单位：Bohr⁻¹），ABACUS 会自动生成对应的 k 网格。等效细分数约为 |**b**| / kspacing，其中 **b** 为倒格矢长度。kspacing 越小，k 点越密。

### 2.2 路径 k 点格式（能带计算）

非自洽计算（`calculation nscf`）时，用线段模式在高对称路径上采样 k 点：

```
K_POINTS
8             # 高对称点段数（行数 = 此数字）
Line
0.00000000  0.00000000  0.00000000  25  # Γ
0.50000000  0.00000000  0.50000000  9   # X
0.62500000  0.25000000  0.62500000  1   # U|K
0.37500000  0.37500000  0.75000000  27  # K
0.00000000  0.00000000  0.00000000  22  # Γ
0.50000000  0.50000000  0.50000000  18  # L
0.50000000  0.25000000  0.75000000  12  # W
0.50000000  0.00000000  0.50000000  1   # X
```

每行格式：`kx  ky  kz  N`，其中 `kx ky kz` 是高对称点的分数坐标（约化坐标），`N` 是从该点到下一个高对称点之间插入的 k 点数（最后一行的 N 通常写 1）。

上例是 FCC 晶体的 Γ-X-U-K-Γ-L-W-X 路径，来自 `atomkit` 工具自动生成。坐标依赖于晶体对称性，建议用 atomkit 或 ASE 自动生成，避免手写出错。

**生成工具**：`atomkit` 的"Generate K-Mesh & K-Path"功能可根据 STRU 文件自动确定高对称路径并输出 KPT 文件：

```bash
echo -e "3\n301\n3\n101 STRU\n0.06" | atomkit
# 3: Generate K-Mesh & K-Path
# 301: 输出 ABACUS 格式
# 3: Kpath for Bulk Structure
# 101 STRU: 读取 STRU 文件
# 0.06: kspacing
```

### 2.3 Al k 点收敛测试

k 点数量需要收敛测试。以 Al FCC 为例，保持 ecutwfc=50 Ry 不变，逐步增大 k 点细分数（n × n × n，n 从 1 到 16），从 `istate.info` 读取最高占据态本征值随 k 点密度的变化：

| k 细分数 n | 最高占据态能量 (eV) |
|-----------|-------------------|
| 1 | -93.8589 |
| 2 | -93.3732 |
| 3 | -93.3681 |
| 4 | -93.3671 |
| 5 | -55.0161 |
| 6 | -55.0160 |
| 7 | -55.0160 |
| 8 | -54.8113 |
| 12 | -54.8056 |
| 14 | -54.8055 |
| 16 | -54.6115 |

> 数据来源：Al 元素晶体的自洽迭代计算与平面波收敛测试及 k 点收敛性测试（ABACUS 3.2.0，Al_ONCV_PBE-1.2.upf）
>
> **注**：表中数值为 `istate.info` 中的最高占据态能量（即费米能级附近的能带本征值），不是 `FINAL_ETOT_IS` 所报告的总能量。两者物理含义不同：最高占据态能量用于判断 k 点是否充分采样费米面；总能量用于判断计算是否收敛到正确的基态。收敛测试时两种指标均可使用，但单位和数量级差异很大，不能混淆。

n=1 到 4 时最高占据态能量约为 -93 eV，说明 k 点严重不足，布里渊区采样过于粗糙，费米面未被正确描述。n=5 后跳至 -55 eV 附近，这是 k 网格密度达到能正确描述 Γ 点附近能带色散的阈值。n=6 之后变化已很小（<0.01 eV），趋于收敛。

**收敛判断**：相邻 k 点密度的能量差 < 1 meV/atom 即认为收敛。从数据看，n≥6 时可认为基本收敛，**实际 Al 计算推荐使用 8×8×8 或更密**。

**金属体系 k 点密度要求更高**。与绝缘体/半导体相比，金属的费米面需要更密的 k 网格才能准确描述，通常需要至少 8×8×8。

KPT 写法：

```
K_POINTS
0
Gamma
8  8  8  0 0 0
```

---

## 第三章 INPUT：计算控制参数

`INPUT` 文件控制计算的全部行为：计算类型、精度、算法选择。文件必须命名为 `INPUT`（不可更改）。

### 3.1 文件语法

```
INPUT_PARAMETERS        # 第一行必须是这个关键字，之前的内容被忽略

# 以 # 或 / 开头的行是注释，会被忽略

suffix   Si             # 参数名与值之间用空格分隔
ecutwfc  60             # 参数名不区分大小写
```

- 第一行必须是 `INPUT_PARAMETERS`，否则文件不被识别
- 参数不区分大小写（`ecutwfc` 和 `ECUTWFC` 等价）
- 布尔值支持多种写法：`1`/`0`、`true`/`false`、`T`/`F`
- 运行结束后，`OUT.suffix/INPUT` 文件会列出所有实际使用的参数值（含默认值），可用于核查

**设置原则**：参数不是越多越好，未设置的参数使用默认值。默认值通常已能满足常见计算，只在有明确需要时才修改。

### 3.2 参数分组与常用参数

INPUT 参数可分为五组：

#### 组一：基本参数（General）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| suffix | ABACUS | 输出目录后缀，结果保存在 OUT.suffix/ |
| calculation | scf | 计算类型，见下表 |
| basis_type | pw | 基组类型：pw 或 lcao |
| pseudo_dir | ./ | 赝势文件目录 |
| orbital_dir | ./ | 轨道文件目录（lcao 专用） |
| symmetry | 1 | 是否利用晶体对称性（推荐开启） |
| ntype | — | 元素种类数（通常可省略，ABACUS 自动从 STRU 读取） |

**calculation 可选值**：

| 值 | 含义 |
|----|------|
| scf | 自洽电子结构计算（Self-Consistent Field） |
| relax | 离子弛豫（固定晶胞，优化原子位置） |
| cell-relax | 变胞弛豫（同时优化晶胞参数和原子位置） |
| md | 分子动力学 |
| nscf | 非自洽计算（读入已有电荷密度，计算能带/DOS） |

#### 组二：SCF 迭代

| 参数 | 默认值 | 说明 |
|------|--------|------|
| ecutwfc | — | 平面波截断能（Ry），PW 计算必填，见 3.5 节 |
| scf_nmax | 100 | SCF 最大迭代步数 |
| scf_thr | 1e-6 (pw) / 1e-7 (lcao) | 电荷密度收敛阈值（Ry），两步之间密度差 Δρ < scf_thr 时判断收敛 |

**scf_thr 建议值**：
- 一般计算：`1e-6`（PW）或 `1e-7`（LCAO）
- 需要高精度（如声子、弹性常数）：`1e-8`

#### 组三：KS 方程求解

| 参数 | 默认值 | 说明 |
|------|--------|------|
| nbands | 自动 | 计算的能带数，默认根据价电子数自动设置，需要空带（能带、DOS）时手动调大 |
| ks_solver | cg (pw) / genelpa (lcao) | KS 方程的对角化算法 |
