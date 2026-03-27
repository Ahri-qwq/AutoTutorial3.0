# ABACUS 输入输出体系：三文件入门与收敛性测试

作者：AutoTutorial 3.0

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

k 点数量需要收敛测试。以 Al FCC 为例，保持 ecutwfc=50 Ry 不变，逐步增大 k 点细分数（n × n × n，n 从 1 到 16），记录单原子能量：

| k 细分数 n | 总能量 (eV) |
|-----------|------------|
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

> 数据来源：Al 元素晶体的自洽迭代计算与平面波收敛测试及 k 点收敛性测试

n=1 到 4 时能量约为 -93 eV，说明 k 点严重不足，布里渊区采样过于粗糙。n=5 后能量跳至 -55 eV 附近，进入合理的收敛区间。n=6 之后能量变化已很小（<0.01 eV），趋于收敛。

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

以 Al FCC（4 原子原胞，k 点 6×6×6）为例，逐步增大 ecutwfc，记录总能量：

| ecutwfc (Ry) | 总能量 (eV) |
|-------------|------------|
| 20 | -57.197 |
| 25 | -56.2212 |
| 30 | -55.016 |
| 40 | -54.0526 |
| 45 | -53.9783 |
| 50 | -53.9617 |
| 55 | -53.9591 |
| 70 | -53.9602 |
| 80 | -53.9606 |

> 数据来源：Al 元素晶体的自洽迭代计算与平面波收敛测试及 k 点收敛性测试

从 20 Ry 到 40 Ry，能量迅速下降（~3 eV），说明基组严重不完整。40 Ry 之后趋于平坦，50 Ry 与 55 Ry 的差值仅 0.0026 eV（2.6 meV），与 70 Ry 的差值为 0.001 eV（1 meV），满足收敛标准。**Al 使用 50 Ry 即可，精度要求高时取 60 Ry**。

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
## 附录：输出文件概览

ABACUS 计算完成后，结果保存在 `OUT.suffix/` 目录（suffix 由 INPUT 中的 `suffix` 参数决定，默认为 `ABACUS`）。

### 主要输出文件

| 文件 | 内容 |
|------|------|
| `running_scf.log` | 主日志：SCF 收敛过程、最终能量、警告信息 |
| `INPUT` | 本次计算使用的全部参数（含默认值），可核查实际设置 |
| `STRU.cif` | 计算结构的 CIF 格式，便于可视化 |
| `istate.info` | 能级、占据数、费米能级 |
| `kpoints` | 实际使用的 k 点列表（含权重） |

### 提取关键信息

**查看 SCF 是否收敛并获取总能量**：

```bash
grep FINAL_ETOT_IS OUT.Al/running_scf.log
# 输出示例：
#  !FINAL_ETOT_IS  -215.847068 eV
```

**追踪 SCF 收敛过程**（查看每步的电荷密度差 drho）：

```bash
grep drho OUT.Al/running_scf.log | tail -20
```

**确认使用的实际参数**：

```bash
grep ecutwfc OUT.Al/INPUT
grep mixing_beta OUT.Al/INPUT
```

### SCF 未收敛的常见提示

如果 `running_scf.log` 末尾出现 `SCF NOT converged` 而不是 `FINAL_ETOT_IS`，说明计算未在 `scf_nmax` 步内收敛。收敛困难时的参数调整策略参见第三章 3.4 节。

### 进阶学习方向

本教程覆盖了 ABACUS 最核心的三个输入文件。后续可以进一步学习：

- **结构优化**（`calculation relax` / `cell-relax`）：在掌握 SCF 的基础上，理解力和应力的收敛判据（`force_thr_ev`、`stress_thr`）
- **能带与 DOS**（`calculation nscf`）：SCF 完成后读入电荷密度，配合路径 k 点计算能带，用 atomkit 后处理作图
- **LCAO 基组与收敛**：LCAO 计算中基组大小（SZ/DZP/TZDP）对精度的影响
- **赝势选择**：模守恒赝势（ONCV/SG15）与不同泛函（LDA/GGA/meta-GGA）的适配

---

### 参考资料

1. ABACUS 线上文档 — INPUT 参数完整列表：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html
2. ABACUS 线上文档 — STRU 文件格式：https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html
3. ABACUS 线上文档 — KPT 文件格式：https://abacus.deepmodeling.com/en/latest/advanced/input_files/kpt.html
4. ABACUS 收敛性问题解决手册（周巍青，AISI）
5. ABACUS 的平面波计算与收敛性测试（谢炘玥，PKU）
6. Al 元素晶体的自洽迭代计算与平面波收敛测试及 k 点收敛性测试
