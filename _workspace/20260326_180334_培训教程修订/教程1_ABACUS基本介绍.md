---
title: "ABACUS 基本介绍：输入输出体系与结构优化"
author: "AutoTutorial 3.0"
date: "2026-03-26"
topic: "ABACUS 基本介绍"
task_type: "D"
has_case: true
word_count: ~6500
chapters: 7
---

# ABACUS 基本介绍：输入输出体系与结构优化

> 本文由 AutoTutorial 3.0 提供。

---

## 第零章 ABACUS 软件背景

ABACUS（Atomic-orbital Based Ab-initio Computation at UStc，原子算筹）是中国科学技术大学开发的开源第一性原理计算软件。该软件基于密度泛函理论（DFT），支持平面波（PW）和数值原子轨道（LCAO）两种基组，可用于材料的电子结构、力学性质、光学性质等多种物理性质的计算。

### 软件特点

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

### 与其他软件的对比

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

### 适用场景

ABACUS 特别适合以下场景：
- 需要开源软件进行学术研究和教学
- 大体系计算（数百原子，使用 LCAO 基组）
- 需要中文文档和社区支持
- 探索新算法和方法开发（开源代码便于二次开发）

---

## 前言

本教程介绍 ABACUS 的基本使用，包括软件背景、三个核心输入文件（STRU、KPT、INPUT）、结构优化、输出文件解读，以及 abacustest 工具的使用。

**教程目标**

- 了解 ABACUS 软件的背景和特点
- 理解 STRU、KPT、INPUT 各自承担什么职责
- 学会为平面波（PW）和数值原子轨道（LCAO）两种基组写 STRU
- 掌握 k 点采样的两种格式：MP 网格（SCF 用）和路径 k 点（能带用）
- 理解 INPUT 中的关键参数分组，重点掌握 mixing 参数的物理含义
- 掌握结构优化（relax 和 cell-relax）的原理和参数设置
- 学会解读 ABACUS 输出文件（running_scf.log）
- 掌握 abacustest 工具准备输入和抓取结果的方法

**适用读者**

初次使用 ABACUS 的用户，或熟悉其他 DFT 软件（VASP、QE）需要迁移的用户。建议具备基本的 DFT 概念（自洽场计算、布里渊区采样）。

**前置知识**

- 晶体结构基本概念（晶格常数、分数坐标）
- 自洽场（SCF）计算的基本流程

**教程结构**

| 章节 | 主题 |
|------|------|
| 第零章 | ABACUS 软件背景 |
| 第一章 | STRU 结构文件 |
| 第二章 | KPT k 点文件 |
| 第三章 | INPUT 计算参数 |
| 第四章 | 结构优化 |
| 第五章 | 输出文件解读 |
| 第六章 | abacustest 准备输入 |
| 第七章 | abacustest 抓取结果 |

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

#### 组四：展宽（Smearing）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| smearing_method | gauss | 展宽方法，见下方说明 |
| smearing_sigma | 0.015 Ry | 展宽宽度（Ry），约 0.2 eV |

**smearing_method 选择**：

| 方法 | 适用体系 |
|------|---------|
| gauss / gaussian | 通用，适合绝缘体和金属 |
| mp | 金属体系推荐，精度更高 |
| fixed | 绝缘体（不展宽，需要带隙严格大于 smearing_sigma） |
| fd | Fermi-Dirac，用于高温计算（SDFT） |

**smearing_sigma 建议**：
- 绝缘体：0.01 Ry 以下
- 金属：0.02 Ry，若收敛困难可适当增大到 0.05 Ry

#### 组五：电荷密度混合（Mixing）

混合参数控制 SCF 如何从旧电荷密度更新到新电荷密度，是影响 SCF 收敛速度和稳定性的关键。详见 3.4 节。

---

### 3.3 一个完整的 Si SCF 示例

以 8 原子金刚石结构 Si 为例，下面是参数设置和注释：

```
INPUT_PARAMETERS
# === 基本参数 ===
suffix                  Si          # 输出到 OUT.Si/
calculation             scf
symmetry                1
pseudo_dir              ./          # 赝势目录
basis_type              pw

# === SCF 迭代 ===
ecutwfc                 60          # 截断能，60 Ry 对 Si 已收敛
scf_nmax                100
scf_thr                 1e-8        # 严格收敛阈值

# === KS 方程求解 ===
nbands                  26          # Si 8 原子：4 价电子/原子×8=32 个电子，需 16 条占据带+空带

# === 展宽 ===
smearing_method         gauss
smearing_sigma          0.01        # Si 是半导体，展宽可以小些

# === 混合 ===
mixing_type             broyden
mixing_beta             0.7         # Si 为半导体，mixing_beta 取 0.7 即可
mixing_gg0              0           # 绝缘体关闭 Kerker 预处理
```

---

### 3.4 电荷密度混合参数

SCF 迭代通过"混合"旧电荷密度 ρ_in 和本步算出的 ρ_out，生成下一步的输入电荷密度。混合策略直接影响收敛速度和稳定性。

#### mixing_type

选择混合算法。

| 值 | 算法 | 说明 |
|----|------|------|
| broyden | Broyden 准牛顿法 | 默认，收敛最快 |
| pulay | Pulay DIIS | 略慢于 broyden |
| plain | 线性混合 | 最慢，仅用于调试 |

无特殊理由保持默认 `broyden`。

#### mixing_beta

新电荷密度占混合结果的比例，取值范围 0~1。

```
ρ_new = mixing_beta × ρ_out + (1 - mixing_beta) × ρ_in
```

- 默认值：0.8（nspin=1），0.4（nspin=2/4）
- **越大，收敛越快，但不收敛的风险越高**
- **越小，收敛越慢，但更稳定**

经验选取：

| 体系 | 推荐值 |
|------|--------|
| 绝缘体 / 半导体（带隙 > 1 eV） | 0.7 |
| 金属 / 过渡金属 | 0.2~0.4 |
| 收敛困难的体系 | 逐步减小，试 0.3、0.2、0.1 |

#### mixing_ndim

保存历史电荷密度的步数，默认 8。Broyden 和 Pulay 算法利用历史步信息构建更好的更新方向。增大 mixing_ndim 可提升收敛稳定性，代价是内存消耗线性增加。

#### mixing_gg0（Kerker 预处理）

默认 1.0（开启）。Kerker 方法在倒空间中衰减长波长（低频）的电荷密度分量，对金属体系收敛帮助很大。

```
mixing_gg0   1.0    # 金属体系（默认，推荐保持）
mixing_gg0   0.0    # 绝缘体/原子/分子体系（可关闭）
```

#### 收敛困难时的调参策略

| 情形 | 建议操作 |
|------|---------|
| drho 振荡不降 | 减小 mixing_beta（试 0.3、0.2、0.1） |
| 收敛太慢（>100 步）| 增大 mixing_ndim（试 12、20） |
| 金属收敛困难 | 确认 mixing_gg0=1.0，减小 smearing_sigma |
| 绝缘体/原子体系 | 关闭 Kerker：mixing_gg0=0.0 |
| DFT+U 体系 | 尝试开启 mixing_dmr，配合 mixing_restart 1e-3 |

---

### 3.5 平面波截断能（ecutwfc）与收敛测试

`ecutwfc` 控制展开波函数的平面波个数上限，单位 Ry。ecut 越大，基组越完整，计算结果越精确，但计算量与 ecutwfc^(3/2) 成正比，需要在精度和效率之间权衡。

**只有 PW 计算才需要测试 ecutwfc**。LCAO 计算中 ecutwfc 仅控制辅助格点密度，不决定基组大小。

#### Al FCC 截断能收敛数据

以 Al FCC（4 原子原胞，k 点 6×6×6）为例，逐步增大 ecutwfc，记录最高占据态能量（来自 `istate.info`）：

| ecutwfc (Ry) | 最高占据态能量 (eV) |
|-------------|-------------------|
| 20 | -57.197 |
| 25 | -56.2212 |
| 30 | -55.016 |
| 40 | -54.0526 |
| 45 | -53.9783 |
| 50 | -53.9617 |
| 55 | -53.9591 |
| 70 | -53.9602 |
| 80 | -53.9606 |

> 数据来源：Al 元素晶体的自洽迭代计算与平面波收敛测试及 k 点收敛性测试（ABACUS 3.2.0，Al_ONCV_PBE-1.2.upf）
>
> **注**：表中数值为 `istate.info` 中的最高占据态能量，不是 `FINAL_ETOT_IS` 所报告的总能量。使用不同版本赝势（如 Al_ONCV_PBE-1.0.upf）或不同 ABACUS 版本时，总能量绝对值会有所不同，但收敛趋势和推荐截断能不变。本教程中 Al_ONCV_PBE-1.0.upf 实际计算的总能量约为 -1883 eV/atom（含11个价电子）。

从 20 Ry 到 40 Ry，最高占据态能量快速上升（约 3 eV），说明基组严重不完整，截断能不足导致波函数描述失真。40 Ry 之后趋于平坦，50 Ry 与 55 Ry 的差值仅 0.0026 eV（2.6 meV），与 70 Ry 的差值为 0.001 eV（1 meV），满足收敛标准。**Al 使用 50 Ry 即可，精度要求高时取 60 Ry**。

**收敛判断标准**：相邻 ecutwfc 的能量差 < 1 meV/atom。

#### Si 的截断能收敛参考

Si 8 原子体系在 ecut=60 Ry 时收敛（与 50 Ry 差 < 1 meV/atom）。

#### 各元素 ecutwfc 推荐范围

ecutwfc 的合理范围主要由赝势的"硬度"决定：

| 元素类型 | 典型 ecutwfc (Ry) |
|----------|-----------------|
| 简单金属（Al、Si、Na）| 40~60 |
| 过渡金属（Fe、Ni）| 60~100 |
| 含 O、N 的氧化物 | 60~80 |
| 第一行元素（C、N、O）| 60~80 |

具体值应以收敛测试为准，上表仅供参考。

#### 批量提交收敛测试的脚本思路

```bash
# 依次运行不同 ecutwfc 的 SCF 计算
for ecut in 20 30 40 50 60 70 80; do
    mkdir -p ecut_${ecut}
    cp STRU KPT ecut_${ecut}/
    sed "s/ecutwfc.*/ecutwfc  ${ecut}/" INPUT > ecut_${ecut}/INPUT
    cd ecut_${ecut}
    mpirun -n 4 abacus > log 2>&1
    cd ..
done

# 提取各 ecutwfc 下的总能量
for ecut in 20 30 40 50 60 70 80; do
    echo -n "ecutwfc=${ecut}: "
    grep FINAL_ETOT_IS ecut_${ecut}/OUT.Al/running_scf.log
done
```

---

### 3.6 Al FCC 完整 INPUT 示例

综合以上讨论，Al FCC SCF 计算的 INPUT 文件：

```
INPUT_PARAMETERS
# === 基本参数 ===
suffix                  Al
calculation             scf
symmetry                1
pseudo_dir              ./

# === 基组 ===
basis_type              pw
ecutwfc                 50          # Al 收敛值

# === SCF 迭代 ===
scf_nmax                100
scf_thr                 1e-8

# === 展宽（Al 是金属，用 mp 展宽）===
smearing_method         mp
smearing_sigma          0.02

# === 混合（金属，mixing_beta 取较小值）===
mixing_type             broyden
mixing_beta             0.4
mixing_gg0              1.0         # 金属开启 Kerker 预处理
```

---


---

## 第四章 结构优化

# 结构优化（几何优化）

第一性原理计算中的"结构优化"，本质上是在势能面（Born-Oppenheimer 势能面）上寻找极小值点的过程。给定初始结构，程序计算各原子受力（$F_i = -\partial E / \partial \mathbf{r}_i$）和晶胞应力张量，沿受力方向更新原子坐标。对于 cell-relax，还同步更新晶格参数。重复迭代，直到力与应力均低于设定阈值。

## 两种优化模式

ABACUS 支持两种结构优化模式：

| 计算类型 | INPUT 设置 | 优化自由度 | 典型用途 |
|---------|-----------|-----------|---------|
| 固定晶胞优化（ionic relax） | `calculation relax` | 原子位置 | 固定实验晶格参数，优化原子位置 |
| 变胞优化（variable-cell relax） | `calculation cell-relax` | 原子位置 + 晶格参数 | 充分弛豫后用于弹性/声子等性质计算 |

**cell-relax 的嵌套流程：**

ABACUS 的变胞优化采用嵌套迭代策略：先在当前晶胞约束下做固定晶胞的原子弛豫（一轮 RELAX IONS），再根据应力张量更新晶格参数（RELAX CELL 计数加一），如此循环直到力和应力同时收敛。输出日志中的标记格式为：

```
-------------------------------------------
RELAX CELL : 3
RELAX IONS : 1 (in total: 15)
-------------------------------------------
```

这表示第 3 次晶胞更新后的第 1 个（共 15 个）离子步。

---

## 收敛判据

结构优化由两个阈值控制：

| 参数 | 单位 | 含义 | 典型推荐值 |
|------|------|------|-----------|
| `force_thr_ev` | eV/Å | 所有原子受力中最大分量（LARGEST GRAD）的收敛阈值 | 0.01–0.05 |
| `stress_thr` | kBar | TOTAL-PRESSURE 绝对值的收敛阈值（负值表示拉伸压力，取绝对值与阈值比较） | 0.5–5 |

当 LARGEST GRAD < `force_thr_ev` **且** TOTAL-PRESSURE 绝对值 < `stress_thr` 同时满足时，计算收敛。对于弹性常数计算，残余应力应尽量小，建议 `stress_thr` 取 0.5 kBar 或更严格。

其余常用控制参数：

| 参数 | 含义 | 默认/推荐值 |
|------|------|-----------|
| `relax_nmax` | 最大离子步数 | 100 |
| `cal_force` | 开启受力计算 | 1（relax/cell-relax 时必须） |
| `cal_stress` | 开启应力张量计算 | 1（cell-relax 时必须） |
| `out_stru` | 输出优化后结构文件（STRU_ION_D） | 1 |

---

## 案例：h-BN 变胞优化

六方氮化硼（h-BN）是典型的层状材料，属六方晶系，结构类似石墨。层内 B-N 键强，层间为弱 van der Waals 相互作用。层间距离对 c 方向应力非常敏感，cell-relax 需要同时收紧力和应力才能得到可靠的平衡结构。

本案例体系：**h-BN 192原子超胞**（B₉₆N₉₆），LCAO 基组，PBE 泛函。

### INPUT 文件

```
INPUT_PARAMETERS
suffix                  BN
pseudo_dir              ./
orbital_dir             ./
calculation             cell-relax
symmetry                0
basis_type              lcao
ecutwfc                 100
scf_thr                 1e-07
scf_nmax                100
smearing_method         gauss
smearing_sigma          0.002
mixing_type             pulay
mixing_beta             0.3
cal_force               1
cal_stress              1
force_thr_ev            0.01
stress_thr              0.5
relax_nmax              100
out_stru                1
kspacing                0.08
```

参数说明：

| 参数 | 值 | 说明 |
|------|------|------|
| `calculation` | `cell-relax` | 同时优化晶胞和原子位置 |
| `symmetry` | `0` | 关闭对称性。层状结构在弛豫过程中可能降低对称性，关闭后更安全 |
| `cal_force` | `1` | 开启受力计算，cell-relax 必须 |
| `cal_stress` | `1` | 开启应力张量计算，cell-relax 必须 |
| `force_thr_ev` | `0.01` | 力收敛阈值，0.01 eV/Å 适用于大多数材料 |
| `stress_thr` | `0.5` | 应力收敛阈值，0.5 kBar 较严格，弹性计算前建议使用 |
| `kspacing` | `0.08` | 自动生成 K 点网格，代替手写 KPT 文件，单位 Å⁻¹ |
| `smearing_method` | `gauss` | 高斯展宽，处理绝缘体时展宽量可小 |
| `smearing_sigma` | `0.002` | 展宽参数（Ry），h-BN 是绝缘体，取小值 |
| `mixing_type` | `pulay` | Pulay 混合，适合 LCAO 计算 |
| `mixing_beta` | `0.3` | 混合系数，大体系取 0.3–0.4 合适 |

> **注意：** `scf_thr` 设为 `1e-7` 比默认值严格，确保每一步 SCF 充分收敛，减少应力计算误差。结构优化的精度上限由 SCF 精度决定。

---

### 运行与收敛过程

提交计算后，ABACUS 的输出日志（`running_cell-relax.log`）中会依次输出每轮 RELAX CELL 和 RELAX IONS 的信息：

```
-------------------------------------------
RELAX CELL : 1
RELAX IONS : 1 (in total: ...)
-------------------------------------------
...
  LARGEST GRAD (eV/A)  :      0.312
  TOTAL-PRESSURE: -2.070e+00 KBAR     <-- 绝对值 2.07 > 0.5，未收敛
-------------------------------------------
RELAX CELL : 2
RELAX IONS : 1 (in total: ...)
-------------------------------------------
...
  LARGEST GRAD (eV/A)  :      0.054
  TOTAL-PRESSURE: -8.500e-01 KBAR     <-- 绝对值 0.85 > 0.5，未收敛
-------------------------------------------
RELAX CELL : 3
RELAX IONS : 1 (in total: ...)
-------------------------------------------
...
  LARGEST GRAD (eV/A)  :      0.008   <-- 小于 force_thr_ev=0.01，满足
  TOTAL-PRESSURE: -4.350e-01 KBAR     <-- 绝对值 0.435 < 0.5，满足 → 收敛！
```

判断收敛的两个关键量：

- **LARGEST GRAD**：所有原子受力分量中的最大值（对应 `force_thr_ev`）
- **TOTAL-PRESSURE**：应力张量均值（对应 `stress_thr`，取绝对值比较）

当两者**同时**低于设定阈值时，计算收敛，输出最终优化结构 `OUT.BN/STRU_ION_D`（通过 `out_stru 1` 开启）。

---

## 常见问题

**Q：cell-relax 跑了很多步但应力一直不降？**

可能原因：初始结构残余应力太大，或 `kspacing` 取值偏大（K 点密度不够）。建议检查 K 点收敛性，或适当减小 `kspacing`（如从 0.1 改为 0.08）。

**Q：`symmetry 0` 和 `symmetry 1` 有什么区别？**

开启对称性（`symmetry 1`）会约束原子位置使其满足对称性，减少独立自由度，加速收敛，但有时会因程序识别对称性不准确导致弛豫出错。大体系或低对称性材料建议用 `symmetry 0`。

**Q：弛豫完成后如何继续进行弹性计算？**

将 `OUT.BN/STRU_ION_D` 复制为新目录的 `STRU`，修改 `INPUT` 中 `calculation scf`（或交给 abacustest 自动处理），即可在优化后结构上计算弹性常数。

---


---

## 第五章 输出文件解读

# 输出文件解读

ABACUS 计算完成后，会在 `OUT.<suffix>/` 目录下生成多个输出文件。本节重点介绍最常用的两个文件：屏幕输出（通常重定向到 `log` 文件）和 `running_scf.log`。

## 屏幕输出（log 文件）

运行 ABACUS 时，通常将屏幕输出重定向到文件：
```bash
mpirun -np 8 abacus > log
```

`log` 文件包含完整的计算过程信息，主要内容：

**1. 版本和编译信息**
```
WELCOME TO ABACUS v3.10.1
Compiled at Mar 15 2024 10:23:45
```

**2. 输入参数回显**
程序会打印读取到的所有 INPUT 参数，便于检查是否正确读取。

**3. 结构信息**
```
READING STRUCTURAL INFORMATION
lattice constant (Bohr) = 1.88973
lattice vectors:
  4.03460  0.00000  0.00000
  0.00000  4.03460  0.00000
  0.00000  0.00000  4.03460
```

**4. k 点信息**
```
KPOINTS INFORMATION
number of k-points = 64
k-point coordinate:
  0.0625  0.0625  0.0625  weight = 0.015625
  ...
```

**5. 自洽迭代过程**（见下节 running_scf.log）

**6. 最终结果**
```
!FINAL_ETOT_IS -241.234567890123 eV
Time Statistics:
  Total Time: 123.45 s
```

---

## running_scf.log 文件

`OUT.<suffix>/running_scf.log` 记录 SCF 自洽迭代的详细过程，是判断计算是否收敛的关键文件。

### 文件结构

每个 SCF 步包含以下信息：

```
ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)
GE1    -241.12345678  0.00e+00       1.23e-01   2.34
GE2    -241.23456789  -1.11e-01      3.45e-02   2.56
GE3    -241.23456790  -1.00e-06      8.76e-04   2.67
...
GE15   -241.23456790  -2.34e-09      5.67e-09   2.89
```

**列含义：**
- `ITER`：迭代步数（GE = Ground state Energy）
- `ETOT(eV)`：当前步的总能量
- `EDIFF(eV)`：能量变化量（当前步 - 上一步）
- `DRHO`：电荷密度变化量（衡量自洽程度）
- `TIME(s)`：该步耗时

### 收敛判据

SCF 收敛需要同时满足：
1. `|EDIFF|` < `scf_thr`（INPUT 中设置，默认 1e-6 eV）
2. `DRHO` 足够小（通常 < 1e-6）

**收敛标志：**
```
charge density convergence is achieved
time of this iteration = 2.89 s
```

**未收敛警告：**
```
WARNING: SCF does not converge in 100 steps!
```

### 典型收敛曲线

**正常收敛：**
```
GE1    -241.12345678  0.00e+00       1.23e-01   2.34
GE2    -241.23456789  -1.11e-01      3.45e-02   2.56
GE3    -241.23456790  -1.00e-06      8.76e-04   2.67
GE4    -241.23456790  -5.00e-09      2.34e-07   2.78
GE5    -241.23456790  -1.23e-10      5.67e-09   2.89  ← 收敛
```
能量和电荷密度单调下降，5-20 步收敛。

**震荡不收敛：**
```
GE1    -241.12345678  0.00e+00       1.23e-01   2.34
GE2    -241.23456789  -1.11e-01      3.45e-02   2.56
GE3    -241.20000000  +3.46e-02      5.67e-02   2.67  ← 能量上升
GE4    -241.22000000  -2.00e-02      4.56e-02   2.78
...
```
能量反复震荡，DRHO 不下降。

**解决方法：**
- 减小 `mixing_beta`（从 0.7 降到 0.4 或 0.2）
- 增大 `smearing_sigma`（金属体系）
- 检查 k 点密度是否足够
- 检查初始结构是否合理

---

## 其他常用输出文件

| 文件名 | 内容 | 用途 |
|--------|------|------|
| `STRU_ION_D` | 优化后的结构 | 结构优化后读取 |
| `BANDS_1.dat` | 能带数据 | 能带结构计算 |
| `TDOS` | 总态密度 | DOS 计算 |
| `PDOS` | 投影态密度 | PDOS 计算 |
| `SPIN1_CHG.cube` | 电荷密度 | 可视化、NSCF 输入 |
| `mulliken.txt` | Mulliken 布居分析 | 电荷分析 |

---

## 快速检查计算结果

**提取最终能量：**
```bash
grep "FINAL_ETOT_IS" OUT.*/running_scf.log
```

**检查是否收敛：**
```bash
grep "convergence" OUT.*/running_scf.log
```

**查看最后几步迭代：**
```bash
tail -20 OUT.*/running_scf.log
```

---


---

## 第六章 使用 abacustest 准备输入文件

# 使用 abacustest 准备输入文件

abacustest 是 ABACUS 的前后处理工具，可以快速从结构文件（CIF、POSCAR 等）生成完整的 ABACUS 输入文件夹。

## 安装 abacustest

```bash
# 方法1：通过 pip 安装
pip install abacustest

# 方法2：从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装后验证：
```bash
abacustest -h
```

---

## 从单个结构文件准备输入

### 基本用法

假设你有一个 CIF 文件 `Si.cif`，想生成 ABACUS 输入文件：

```bash
abacustest model inputs prepare \
  -f Si.cif \
  --ftype cif \
  --pp /path/to/pseudopotentials \
  --orb /path/to/orbitals
```

**参数说明：**
- `-f`：结构文件路径（支持 CIF、POSCAR、XYZ 等）
- `--ftype`：文件类型（cif、poscar、xyz 等）
- `--pp`：赝势库目录
- `--orb`：轨道库目录（LCAO 计算需要，PW 计算可省略）

**赝势和轨道库要求：**
- 文件名必须以元素名开头，例如 `Si_ONCV_PBE-1.0.upf`、`Si_gga_8au_60Ry_2s2p1d.orb`
- 或者在目录下提供 `element.json` 文件，格式：`{"Si": "Si_ONCV_PBE-1.0.upf"}`

### 生成结果

执行后会在当前目录生成 `Si/` 文件夹，包含：
- `STRU`：结构文件
- `INPUT`：默认参数文件
- `KPT`：默认 k 点文件
- 赝势和轨道文件的软链接或副本

---

## 批量准备多个任务

对于多个结构文件，可以使用 `param.json` 配置文件批量准备。

### 创建 param.json

```json
{
  "prepare": {
    "strus": ["Si.cif", "Ge.cif", "GaAs.cif"],
    "stru_format": "cif",
    "input_template": "INPUT_template",
    "kpt_template": "KPT_template",
    "pp_dict": {
      "Si": "Si_ONCV_PBE-1.0.upf",
      "Ge": "Ge_ONCV_PBE-1.0.upf",
      "Ga": "Ga_ONCV_PBE-1.0.upf",
      "As": "As_ONCV_PBE-1.0.upf"
    },
    "orb_dict": {
      "Si": "Si_gga_8au_60Ry_2s2p1d.orb",
      "Ge": "Ge_gga_8au_100Ry_2s2p1d.orb",
      "Ga": "Ga_gga_8au_100Ry_2s2p2d.orb",
      "As": "As_gga_7au_100Ry_2s2p1d.orb"
    },
    "pp_path": "/path/to/pseudopotentials",
    "orb_path": "/path/to/orbitals"
  }
}
```

### 准备 INPUT 和 KPT 模板

**INPUT_template：**
```
INPUT_PARAMETERS
suffix              ABACUS
calculation         scf
basis_type          lcao
ecutwfc             100
scf_thr             1e-7
scf_nmax            100
smearing_method     gauss
smearing_sigma      0.01
mixing_type         broyden
mixing_beta         0.4
```

**KPT_template：**
```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

### 执行批量准备

```bash
abacustest prepare -p param.json -s abacustest_jobs
```

生成结果：
```
abacustest_jobs/
├── Si/
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   └── ...
├── Ge/
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   └── ...
└── GaAs/
    ├── INPUT
    ├── STRU
    ├── KPT
    └── ...
```

---

## 自动设置 ecutwfc

如果赝势库目录下有 `ecutwfc.json` 文件，abacustest 会自动为每个体系设置合适的 ecutwfc：

**ecutwfc.json 示例：**
```json
{
  "Si": 60,
  "Ge": 80,
  "Ga": 100,
  "As": 100
}
```

abacustest 会自动选择体系中所有元素的最大值作为 ecutwfc。

---

## 常见问题

**Q1：找不到赝势或轨道文件**
- 检查文件名是否以元素名开头
- 或创建 `element.json` 文件指定映射关系

**Q2：生成的 INPUT 参数不符合需求**
- 修改 `INPUT_template` 模板文件
- 或生成后手动修改各任务的 INPUT 文件

**Q3：k 点密度不合适**
- 修改 `KPT_template` 模板文件
- 或使用 `kspacing` 参数自动生成 k 点

---


---

## 第七章 使用 abacustest 抓取计算结果

# 使用 abacustest 抓取计算结果

ABACUS 计算完成后，结果分散在多个输出文件中。abacustest 提供了便捷的后处理功能，可以快速提取关键结果。

## 基本用法

```bash
abacustest model inputs post -j <job_directory>
```

这个命令会自动识别计算类型（SCF、relax、band、DOS 等），提取相应的关键结果。

---

## 提取 SCF 计算结果

对于自洽计算，abacustest 可以提取：
- 最终总能量
- 收敛步数
- 计算时间
- 费米能级（金属）或带隙（半导体）

**示例：**
```bash
cd Si_scf/
abacustest model inputs post -j .
```

**输出示例：**
```json
{
  "final_energy_eV": -241.234567890,
  "converged": true,
  "scf_steps": 15,
  "total_time_s": 123.45,
  "fermi_energy_eV": 5.678
}
```

---

## 批量提取多个任务结果

对于批量计算任务，可以一次性提取所有结果：

```bash
abacustest model inputs post -j job1/ job2/ job3/
```

或使用通配符：
```bash
abacustest model inputs post -j */
```

结果会汇总到一个 JSON 文件中，便于后续分析。

---

## 提取结构优化结果

对于 relax 或 cell-relax 计算：

```bash
cd Si_relax/
abacustest model inputs post -j .
```

**提取信息：**
- 优化后的晶格参数
- 优化后的原子坐标
- 最终能量
- 最大受力
- 应力张量（cell-relax）

**输出示例：**
```json
{
  "final_energy_eV": -241.234567890,
  "converged": true,
  "lattice_constant_angstrom": 5.431,
  "max_force_eV_per_angstrom": 0.008,
  "stress_kbar": [0.12, 0.12, 0.12, 0.0, 0.0, 0.0]
}
```

优化后的结构保存在 `OUT.*/STRU_ION_D` 文件中。

---

## 提取能带和 DOS 数据

**能带计算：**
```bash
cd Si_band/
abacustest model inputs post -j .
```

提取的数据包括：
- 能带路径
- 各 k 点的能量本征值
- 费米能级位置

**DOS 计算：**
```bash
cd Si_dos/
abacustest model inputs post -j .
```

提取的数据包括：
- 能量网格
- 总态密度（TDOS）
- 投影态密度（PDOS，如果计算了）

---

## 自定义提取规则

abacustest 支持自定义提取规则。创建 `extract_rules.json`：

```json
{
  "extract": {
    "final_energy": {
      "file": "OUT.*/running_scf.log",
      "pattern": "FINAL_ETOT_IS\\s+([\\-\\d\\.]+)\\s+eV",
      "type": "float"
    },
    "band_gap": {
      "file": "OUT.*/running_scf.log",
      "pattern": "E_bandgap\\s+([\\d\\.]+)\\s+eV",
      "type": "float"
    }
  }
}
```

使用自定义规则：
```bash
abacustest model inputs post -j . --rules extract_rules.json
```

---

## 结果可视化

abacustest 还提供了简单的可视化功能（需要安装 matplotlib）：

**绘制收敛曲线：**
```bash
abacustest plot scf -j Si_scf/
```

**绘制能带结构：**
```bash
abacustest plot band -j Si_band/
```

**绘制态密度：**
```bash
abacustest plot dos -j Si_dos/
```

---

## 常见问题

**Q1：找不到输出文件**
- 确保计算已完成
- 检查 `OUT.*` 目录是否存在

**Q2：提取的数据不完整**
- 检查计算是否正常收敛
- 查看 `running_scf.log` 是否有错误信息

**Q3：批量提取时部分任务失败**
- abacustest 会跳过失败的任务，继续处理其他任务
- 检查失败任务的日志文件

---


---

## 附录：进阶学习方向

本教程覆盖了 ABACUS 的基本使用。后续可以进一步学习：

- **弹性常数计算**：在结构优化基础上，计算材料的力学性质
- **能带与 DOS**：计算电子结构，分析材料的电学性质
- **LCAO 基组与收敛**：LCAO 计算中基组大小对精度的影响
- **高级功能**：DFT+U、隐式溶剂模型、RT-TDDFT 等

### 参考资料

1. ABACUS 线上文档：https://abacus.deepmodeling.com
2. ABACUS GitHub：https://github.com/deepmodeling/abacus-develop
3. abacustest GitHub：https://github.com/pxlxingliang/abacus-test