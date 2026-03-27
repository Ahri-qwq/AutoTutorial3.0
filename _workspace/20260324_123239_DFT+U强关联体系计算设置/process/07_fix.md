---
title: "ABACUS 使用教程｜DFT+U 强关联体系计算：以反铁磁 NiO 为例"
author: "AutoTutorial 3.0"
date: "2026-03-24"
topic: "DFT+U 强关联体系计算"
task_type: "C"
has_case: true
word_count: ~2800
chapters: 6
---

# ABACUS 使用教程｜DFT+U 强关联体系计算：以反铁磁 NiO 为例

## 一、介绍

NiO 是过渡金属氧化物中最经典的强关联体系之一。实验上它是一个绝缘体，能隙约 4 eV，并具有 II 型反铁磁序。若用标准的 LDA 或 GGA 泛函计算，往往得到远小于实验值的能隙，甚至错误地预测为金属态。

原因在于 d 轨道电子的局域化。LDA/GGA 对电子间相互作用的描述是均匀化的，无法准确反映局域 d 电子之间强烈的库仑排斥。DFT+U 方法通过引入 Hubbard U 参数，在 d（或 f）轨道上附加一项平均场型的 Hartree-Fock 修正，抑制自相互作用带来的非物理离域化，从而更准确地描述强关联体系。

ABACUS 在 LCAO 基组（`basis_type lcao`）下实现了 DFT+U 功能，理论细节见参考文献 [1]。本教程以反铁磁 NiO 的 SCF 计算为主线，完整演示：

- 如何在 STRU 文件中设置反铁磁初始磁矩
- 如何在 INPUT 中开启 DFT+U 并配置关键参数
- 如何读取和理解 occupation matrix 输出
- 如何使用 occupation matrix control（omc）功能

**适用读者：** 熟悉 ABACUS 基本 SCF 计算流程，希望对过渡金属体系启用 DFT+U 的用户。

---

## 二、案例背景：为什么 NiO 需要 DFT+U

### 2.1 NiO 的电子结构特点

NiO 晶体中，Ni 以 +2 价存在，电子构型为 [Ar] 3d⁸。8 个 d 电子填充在 5 条 d 轨道上，d 轨道未满，轨道间的局域库仑排斥（on-site Coulomb interaction）很强。这种强局域化的 d 电子是 LDA/GGA 难以准确描述的典型情形。

### 2.2 纯 GGA 的问题

标准 GGA-PBE 计算中，NiO 的 d 轨道能带过度展宽（过离域），导致：

- 预测的能隙显著低于实验值（~4 eV），甚至给出金属态
- 磁矩被低估

这是泛函本身的系统性误差，不是收敛问题。

### 2.3 DFT+U 如何修正

DFT+U 的能量泛函形式为：

$$E_{\rm DFT+U} = E_{\rm DFA} + E_U - E_{\rm dc}$$

其中 $E_{\rm DFA}$ 是标准 LDA/GGA 的能量，$E_U$ 是 Hubbard 修正项，$E_{\rm dc}$ 是双计数项（避免对已在 $E_{\rm DFA}$ 中部分包含的相互作用重复计数）。

修正的核心量是 occupation matrix $n_{I,mm'}^\sigma$，描述原子 $I$ 上 d（或 f）轨道的占据情况：

$$n_{I, m m^{\prime}}^\sigma=\frac{1}{N_{\mathbf{k}}} \sum_{n \mathbf{k}} f_{n \mathbf{k}}^\sigma\left\langle\psi_{n \mathbf{k}}^\sigma\left|\hat{P}_{I, m m^{\prime}}^\sigma\right| \psi_{n \mathbf{k}}^\sigma\right\rangle$$

U 修正对高占据轨道和低占据轨道产生不同方向的能量移动，将它们分开，从而打开能隙、改善磁矩描述。

---

## 三、输入文件准备

本案例计算 II 型反铁磁 NiO，原胞包含 2 个 Ni 原子和 2 个 O 原子，两个 Ni 原子的磁矩大小相等、方向相反。

### 3.1 STRU 文件：反铁磁结构的关键设置

处理反铁磁的关键技巧：**将初始磁矩不同的同种元素定义为不同的 atomic species**。虽然 Ni1 和 Ni2 共用相同的赝势和轨道文件，但分开定义后可以分别设置磁矩方向，也便于在 DFT+U 参数中独立控制。

完整 STRU 文件如下：

```
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf
Ni2 58.693 Ni_ONCV_PBE-1.0.upf
O   15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Ni_gga_9au_100Ry_4s2p2d1f.orb
Ni_gga_9au_100Ry_4s2p2d1f.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
9.6226  0.0000  0.0000    #lattice vector A1 (Bohr)
7.9999  5.3468  0.0000    #lattice vector A2 (Bohr)
7.9999  2.4270  4.7629    #lattice vector A3 (Bohr)

ATOMIC_POSITIONS
Direct

Ni1
2.0
1
0.000000000  0.000000000  0.000000000  0  0  0

Ni2
-2.0
1
0.500000000  0.500000000  0.500000000  0  0  0

O
0.0
2
0.250000000  0.250000000  0.250000000  0  0  0
0.750000000  0.750000000  0.750000000  0  0  0
```

几点说明：

- `ATOMIC_POSITIONS` 采用 Direct（分数坐标）格式
- 每个 species 块中，磁矩写在原子数之前（如 `Ni1` 块的 `2.0`），这是**按 species 统一设置磁矩**的方式
- Ni1 设为 +2.0 μB，Ni2 设为 -2.0 μB，O 设为 0.0
- 原子坐标后的 `0 0 0` 表示固定该原子位置（SCF 计算不需要移动原子）
- Ni1 和 Ni2 共用同一套赝势（`Ni_ONCV_PBE-1.0.upf`）和轨道文件（`Ni_gga_9au_100Ry_4s2p2d1f.orb`）

### 3.2 INPUT 文件：DFT+U 参数逐一解析

完整 INPUT 文件如下：

```
INPUT_PARAMETERS
suffix              NiO
calculation         scf
basis_type          lcao
nspin               2

ecutwfc             100
scf_thr             1e-7
scf_nmax            100
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.4

#Parameter DFT+U
dft_plus_u          1
orbital_corr        2 2 -1
hubbard_u           5.0 5.0 0.0

#输出控制
out_bandgap         1
out_mul             1
out_chg             1
```

**DFT+U 核心三参数：**

| 参数 | 含义 | 本例设置 |
|------|------|----------|
| `dft_plus_u` | DFT+U 总开关；设为 0 时即使配置了其他+U参数也不生效 | 1 |
| `orbital_corr` | 数组，长度为 ntype；指定每种 species 施加+U的 l 量子数，-1 表示不施加 | 2 2 -1 |
| `hubbard_u` | 数组，长度同 orbital_corr；每种 species 的 U 值，单位 eV | 5.0 5.0 0.0 |

本例 `ntype = 3`（Ni1、Ni2、O），`orbital_corr = 2 2 -1` 的含义：

- Ni1：l=2（d 轨道）施加 +U，U = 5.0 eV
- Ni2：l=2（d 轨道）施加 +U，U = 5.0 eV
- O：-1，不施加 +U

**输出控制参数：**

| 参数 | 作用 |
|------|------|
| `out_bandgap 1` | 在日志中输出每个自旋通道的能隙 |
| `out_mul 1` | 输出 Mulliken 布居分析（mulliken.txt），含原子磁矩 |
| `out_chg 1` | 输出电荷密度，并将末步 occupation matrix 写入 onsite.dm |

### 3.3 KPT 文件

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

4×4×4 的 Gamma 中心 k 点网格，适用于 NiO 的 SCF 计算。

### 3.4 文件树总览

```
ABACUS_DFT+U/
├── INPUT
├── STRU
├── KPT
├── Ni_ONCV_PBE-1.0.upf
├── Ni_gga_9au_100Ry_4s2p2d1f.orb
├── O_ONCV_PBE-1.0.upf
└── O_gga_7au_100Ry_2s2p1d.orb
```

---

## 四、运行与结果分析

### 4.1 运行计算

进入工作目录后执行：

```bash
cd ABACUS_DFT+U
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

计算完成后，输出文件位于 `OUT.NiO/` 目录。

### 4.2 关键结果提取

#### 总能量

在 `OUT.NiO/running_scf.log` 中搜索 `FINAL_ETOT`：

```
--------------------------------------------
!FINAL_ETOT_IS -9255.7279034240546025 eV
--------------------------------------------
```

#### 能隙

由于设置了 `out_bandgap 1`，搜索 `E_bandgap`，取最后一次出现的值：

```
E_bandgap    +0.205369322748    +2.794192983776
```

两列分别对应 spin up 和 spin down 通道的能隙（单位 eV）。spin down 通道的能隙约 2.79 eV，相比纯 GGA 有明显改善。

#### 磁矩

搜索 `absolute magnetism`，取最后一次出现：

```
      total magnetism (Bohr mag/cell) = 0.00000000
   absolute magnetism (Bohr mag/cell) = 3.35321634
```

- `total magnetism = 0`：符合反铁磁体系整体无净磁矩的预期
- `absolute magnetism ≈ 3.35`：体系内存在大小相等、方向相反的局域磁矩

**原子分辨磁矩（Mulliken 分析）：**

由于设置了 `out_mul 1`，查看 `OUT.NiO/mulliken.txt`，搜索 `Magnetism`：

```
Total Magnetism on atom  Ni1     1.8268646
Total Magnetism on atom  Ni2    -1.8268646
Total Magnetism on atom  O      -3.6718263e-13
Total Magnetism on atom  O       1.7330755e-13
```

两个 Ni 原子的磁矩大小相等（约 ±1.83 μB）、方向相反，O 原子磁矩接近零，确认得到了反铁磁基态。

### 4.3 Occupation Matrix 的读取与解析

DFT+U 计算的特色输出是 **occupation matrix**，记录每个 +U 原子上 d 轨道的占据情况，是理解 +U 修正效果的关键信息。

在 `running_scf.log` 中搜索以 `L(S)DA+U` 开头的块。每个 SCF 步骤都会输出一次，查看最后一次即可。格式如下：

**头部信息**（说明哪些 species 施加了 +U）：

```
atom_type=0  L=2  chi=0    U=5ev
atom_type=1  L=2  chi=0    U=5ev
```

与 INPUT 中的设定一致：atom_type 0（Ni1）和 1（Ni2）均在 l=2（d 轨道）上施加 U=5 eV。

**各原子的 occupation matrix**：

```
atoms  0          // 原子编号：第一个 Ni 原子（Ni1）
L  2              // +U 的 l 量子数（d 轨道）
zeta  0           // 基组 zeta 指标；Ni 只有一个 d 基组，故 zeta=0
spin  0           // spin up 的 5×5 矩阵
// [5×5 matrix]
spin  1           // spin down 的 5×5 矩阵
// [5×5 matrix]

atoms  1          // 第二个 Ni 原子（Ni2），格式同上
L  2
zeta  0
spin  0
// [5×5 matrix]
spin  1
// [5×5 matrix]
```

d 轨道共有 5 个（m = -2, -1, 0, +1, +2），因此每个自旋通道的 occupation matrix 是 5×5 的实对称矩阵。对角元代表各 d 轨道的占据数，取值在 0 到 1 之间。

**onsite.dm 输出文件：**

由于设置了 `out_chg 1`，ABACUS 将最后一步 SCF 收敛后的 occupation matrix 写入 `OUT.NiO/onsite.dm`，格式与日志中的输出完全相同。该文件在 occupation matrix control 功能中会用到。

---

## 五、进阶：Occupation Matrix Control

### 5.1 背景

DFT+U 计算存在**多解性**问题。对于同一体系，不同的初始条件可能收敛到不同的局域极小值，对应不同的 occupation matrix，从而给出不同的能量和物理性质。

Occupation Matrix Control（`omc`）功能允许用户在 SCF 迭代中固定或引导 occupation matrix，主要用途：

1. 将体系引导到特定电子态（避免陷入非目标极小值）
2. 验证 SCF 结果的稳定性

### 5.2 三种模式

| `omc` 值 | 行为 |
|----------|------|
| 0 | 不使用 omc，标准 DFT+U 计算（默认） |
| 1 | SCF 第一步读入 `initial_onsite.dm`，后续正常更新 occupation matrix |
| 2 | 读入 `initial_onsite.dm`，整个 SCF 过程始终使用此 occupation matrix，不更新 |

### 5.3 使用 omc = 2 验证 NiO 计算

**第一步：** 将收敛后的 occupation matrix 复制到工作目录：

```bash
cp OUT.NiO/onsite.dm ./initial_onsite.dm
```

**第二步：** 在 INPUT 中添加 `omc 2`：

```
#Parameter DFT+U
dft_plus_u    1
orbital_corr  2 2 -1
hubbard_u     5.0 5.0 0.0
omc           2
```

**第三步：** 重新运行 ABACUS。

查看新的 `running_scf.log`，确认每一步的 occupation matrix 保持不变。最终的总能量和磁矩与标准 SCF 完全一致，说明之前的计算已收敛到稳定的反铁磁基态。

### 5.4 initial_onsite.dm 文件格式

`initial_onsite.dm` 的格式与 `onsite.dm` 完全相同，也与日志中 `L(S)DA+U` 块的格式一致。用户可以手动编辑该文件，指定自定义的 occupation matrix，以引导体系进入特定电子态。

### 5.5 进阶：通过 Yukawa Potential 自动确定 U 值

如果不确定 U 的取值，ABACUS 支持在 SCF 中自动计算：

| 参数 | 说明 |
|------|------|
| `yukawa_potential 1` | 开启 Yukawa 势方法，程序内自动计算 U 值 |
| `yukawa_lambda` | 手动指定 screening length；不设则由电子密度 on-the-fly 计算 |

理论细节见参考文献 [1]。

---

## 六、参数速查

**DFT+U 相关参数总览：**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `dft_plus_u` | int | 0 | DFT+U 总开关；1 为开启 |
| `orbital_corr` | int[] | -1 | 各 species 施加+U的 l 量子数；-1 不施加 |
| `hubbard_u` | double[] | 0.0 | 各 species 的 U 值（eV） |
| `omc` | int | 0 | Occupation matrix control 模式（0/1/2） |
| `yukawa_potential` | int | 0 | 是否通过 Yukawa 势自动计算 U |
| `yukawa_lambda` | double | — | Yukawa 势 screening length（手动设置时使用） |
| `onsite_radius` | double | — | 投影函数的截断半径（Bohr） |
| `out_chg` | int | 0 | 设为 1 时输出 onsite.dm |

> `orbital_corr` 和 `hubbard_u` 的数组长度必须等于 ntype（即 STRU 中 ATOMIC_SPECIES 的种数）。

---

## 七、结语

过渡金属氧化物、稀土化合物等强关联体系是第一性原理计算的常见挑战。ABACUS 在 LCAO 基组下实现了 DFT+U 功能，三个核心参数（`dft_plus_u`、`orbital_corr`、`hubbard_u`）即可启用，上手门槛低。

本教程以反铁磁 NiO 为例，演示了完整的计算流程，以及 occupation matrix control 功能的使用方式。遇到问题欢迎在 [ABACUS GitHub](https://github.com/deepmodeling/abacus-develop) 提交 Issue。

---

**参考文献**

[1] Qu X, Xu P, Jiang H, He L, Ren X. DFT+U within the framework of linear combination of numerical atomic orbitals. *The Journal of Chemical Physics*. 2022;156(23):234104. [doi:10.1063/5.0090122](https://doi.org/10.1063/5.0090122)
