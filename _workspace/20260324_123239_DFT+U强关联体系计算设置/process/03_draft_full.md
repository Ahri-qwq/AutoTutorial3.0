## 一、介绍

NiO 是过渡金属氧化物中最经典的强关联体系之一。实验上它是一个绝缘体，能隙约 4 eV，并具有 II 型反铁磁序。但若用标准的 LDA 或 GGA 泛函进行计算，往往得到的能隙远小于实验值，甚至错误地预测为金属态。

问题的根源在于 d 轨道电子的局域化。LDA/GGA 对电子间相互作用的描述是均匀化的，无法准确反映局域 d 电子之间强烈的库仑排斥。DFT+U 方法通过引入 Hubbard U 参数，在 d（或 f）轨道上附加一项平均场型的 Hartree-Fock 修正，抑制自相互作用带来的非物理离域化，从而更准确地描述强关联体系。

ABACUS 在 LCAO 基组（`basis_type lcao`）下实现了 DFT+U 功能，理论细节见参考文献 [1]。本教程以反铁磁 NiO 的 SCF 计算为主线，完整演示：

- 如何在 STRU 文件中设置反铁磁初始磁矩
- 如何在 INPUT 中开启 DFT+U 并配置关键参数
- 如何读取和理解 occupation matrix 输出
- 如何使用 occupation matrix control（omc）功能固定 occupation matrix

**适用读者：** 熟悉 ABACUS 基本 SCF 计算流程，希望了解如何对过渡金属体系启用 DFT+U 的用户。
## 二、案例背景：为什么 NiO 需要 DFT+U

### 2.1 NiO 的电子结构特点

NiO 晶体中，Ni 以 +2 价存在，电子构型为 [Ar] 3d⁸。8 个 d 电子填充在 5 条 d 轨道上，d 轨道未满，轨道间的局域库仑排斥（on-site Coulomb interaction）很强。这种强局域化的 d 电子是 LDA/GGA 难以准确描述的典型情形。

### 2.2 纯 GGA 的问题

在标准 GGA-PBE 计算中，NiO 的 d 轨道能带过度展宽（过离域），导致：

- 预测的能隙显著低于实验值（~4 eV），甚至给出金属态
- 磁矩被低估

这不是收敛问题，而是泛函本身的系统性误差。

### 2.3 DFT+U 如何修正

DFT+U 的能量泛函可以写成：

$$E_{\rm DFT+U} = E_{\rm DFA} + E_U - E_{\rm dc}$$

其中 $E_{\rm DFA}$ 是标准 LDA/GGA 的能量，$E_U$ 是 Hubbard 修正项，$E_{\rm dc}$ 是双计数项（避免对已在 $E_{\rm DFA}$ 中部分包含的相互作用重复计数）。

修正的核心是 occupation matrix $n_{I,mm'}^\sigma$，它描述原子 $I$ 上 $d$（或 $f$）轨道的占据情况：

$$n_{I, m m^{\prime}}^\sigma=\frac{1}{N_{\mathbf{k}}} \sum_{n \mathbf{k}} f_{n \mathbf{k}}^\sigma\left\langle\psi_{n \mathbf{k}}^\sigma\left|\hat{P}_{I, m m^{\prime}}^\sigma\right| \psi_{n \mathbf{k}}^\sigma\right\rangle$$

U 修正会对高占据和低占据的轨道产生不同方向的能量移动，将它们分开，从而打开能隙。
## 三、输入文件准备

### 3.1 STRU 文件：反铁磁结构的关键设置

本案例计算 II 型反铁磁 NiO，原胞中两个 Ni 原子的磁矩大小相等、方向相反。

在 ABACUS 的 STRU 文件中，处理反铁磁的常用做法是**将不同初始磁矩的同种元素定义为不同的 atomic species**。虽然它们使用相同的赝势和轨道文件，但通过分别定义 Ni1 和 Ni2，可以对每种 species 统一设置磁矩方向：

```
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf
Ni2 58.693 Ni_ONCV_PBE-1.0.upf
O   15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Ni_gga_9au_100Ry_4s2p2d1f.orb
Ni_gga_9au_100Ry_4s2p2d1f.orb
O_gga_7au_100Ry_2s2p1d.orb
```

在 `ATOMIC_POSITIONS` 部分，分别为 Ni1 和 Ni2 指定初始磁矩：

```
ATOMIC_POSITIONS
Cartesian_angstrom

Ni1
2.0     // 初始磁矩 +2.0 μB
...     // 原子坐标

Ni2
-2.0    // 初始磁矩 -2.0 μB
...     // 原子坐标

O
0.0
...
```

这里 `Ni1` 和 `Ni2` 共用同一套赝势（`Ni_ONCV_PBE-1.0.upf`）和轨道文件（`Ni_gga_9au_100Ry_4s2p2d1f.orb`），差别仅在于初始磁矩符号相反。

> **说明：** 使用双 species 方案的好处是，DFT+U 的参数（`orbital_corr`、`hubbard_u`）按 species 顺序指定，可以独立控制每种 Ni 的 U 值。两个 Ni species 使用相同的 U 值是物理上正确的——它们本质上是同一种元素。

---

### 3.2 INPUT 文件：DFT+U 参数逐一解析

完整的 INPUT 文件如下（仅列出与 DFT+U 相关的关键部分）：

```
INPUT_PARAMETERS
suffix              NiO
calculation         scf
basis_type          lcao
nspin               2

#Parameter DFT+U
dft_plus_u          1
orbital_corr        2 2 -1
hubbard_u           5.0 5.0 0.0

#输出控制
out_bandgap         1
out_mul             1
out_chg             1
```

**核心三参数：**

| 参数 | 含义 | 本例设置 |
|------|------|----------|
| `dft_plus_u` | DFT+U 总开关，1 为开启，0 为关闭（即使设置了其他+U参数也不会生效） | 1 |
| `orbital_corr` | 数组，长度为 ntype（atomic species 数），指定每种 species 施加 +U 的 l 量子数；-1 表示不施加 | 2 2 -1 |
| `hubbard_u` | 数组，长度与 orbital_corr 相同，单位 eV，指定每种 species 的 U 值 | 5.0 5.0 0.0 |

本例中 `ntype = 3`（Ni1、Ni2、O），`orbital_corr = 2 2 -1` 表示：
- Ni1：在 l=2（d 轨道）施加 +U，U = 5.0 eV
- Ni2：在 l=2（d 轨道）施加 +U，U = 5.0 eV
- O：不施加 +U（-1）

**输出控制参数：**

| 参数 | 作用 |
|------|------|
| `out_bandgap 1` | 在日志中输出每个自旋通道的能隙 |
| `out_mul 1` | 输出 Mulliken 布居分析（mulliken.txt），含原子磁矩 |
| `out_chg 1` | 输出电荷密度及末步的 onsite.dm 文件 |

---

### 3.3 KPT 文件

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

4×4×4 的 Gamma 中心 k 点网格，适用于 NiO 的 SCF 计算。

---

**文件树总览：**

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
## 四、运行与结果分析

### 4.1 运行计算

进入工作目录后，执行：

```bash
cd ABACUS_DFT+U
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

计算完成后，输出文件位于 `OUT.NiO/` 目录。

---

### 4.2 关键结果提取

#### 总能量

在 `OUT.NiO/running_scf.log` 中搜索 `FINAL_ETOT`：

```
--------------------------------------------
!FINAL_ETOT_IS -9255.7279034240546025 eV
--------------------------------------------
```

#### 能隙

由于设置了 `out_bandgap = 1`，搜索 `E_bandgap`，取最后一次出现的值：

```
E_bandgap    +0.205369322748    +2.794192983776
```

两列分别对应 spin up 和 spin down 通道的能隙（单位 eV）。NiO 是反铁磁绝缘体，spin down 通道的能隙约 2.79 eV，与实验相比已有明显改善（纯 GGA 通常给出远小的值甚至金属态）。

#### 磁矩

搜索 `absolute magnetism`，取最后一次出现：

```
      total magnetism (Bohr mag/cell) = 0.00000000
   absolute magnetism (Bohr mag/cell) = 3.35321634
```

- **total magnetism = 0**：符合反铁磁体系整体无净磁矩的预期
- **absolute magnetism ≈ 3.35**：体系内存在大小相等、方向相反的局域磁矩

**原子分辨磁矩（Mulliken 分析）：**

由于设置了 `out_mul = 1`，查看 `OUT.NiO/mulliken.txt`，搜索 `Magnetism`：

```
Total Magnetism on atom  Ni1     1.8268646
Total Magnetism on atom  Ni2    -1.8268646
Total Magnetism on atom  O      -3.6718263e-13
Total Magnetism on atom  O       1.7330755e-13
```

两个 Ni 原子的磁矩大小相等（约 ±1.83 μB）、方向相反，O 原子磁矩接近于零，确认得到了反铁磁基态。

---

### 4.3 Occupation Matrix 的读取与解析

DFT+U 计算的一个特色输出是 **occupation matrix**，它记录了每个 +U 原子上 d 轨道的占据情况，是理解 +U 修正效果的重要信息。

在 `running_scf.log` 中搜索以 `L(S)DA+U` 开头的块，格式如下：

**头部信息**（说明哪些 species 施加了 +U 及其参数）：

```
atom_type=0  L=2  chi=0    U=5ev
atom_type=1  L=2  chi=0    U=5ev
```

与 INPUT 文件中的设定一致：atom_type 0（Ni1）和 1（Ni2）均在 l=2（d 轨道）上施加 U=5 eV。

**各原子的 occupation matrix**：

```
atoms  0          // 原子编号：第一个 Ni 原子（Ni1）
L  2              // +U 的 l 量子数（d 轨道）
zeta  0           // 基组 zeta 指标；Ni 只有一个 d 基组，故 zeta=0
spin  0           // spin up 的 5×5 矩阵
// [5×5 matrix]
spin  1           // spin down 的 5×5 矩阵
// [5×5 matrix]

atoms  1          // 第二个 Ni 原子（Ni2）
L  2
zeta  0
spin  0
// [5×5 matrix]
spin  1
// [5×5 matrix]
```

d 轨道共有 5 个（m = -2, -1, 0, +1, +2），因此每个自旋通道的 occupation matrix 是 5×5 的实对称矩阵。对角元代表各 d 轨道的占据数，0 到 1 之间。

**onsite.dm 输出文件：**

由于设置了 `out_chg 1`，ABACUS 会将最后一步 SCF 收敛后的 occupation matrix 写入 `OUT.NiO/onsite.dm`，文件格式与日志中的输出完全相同。这个文件在使用 occupation matrix control 功能时会用到（见下一章）。
## 五、进阶：Occupation Matrix Control

### 5.1 背景与动机

DFT+U 计算存在一个已知问题：**多解性**。对于同一个体系，不同的初始磁矩或初始电荷密度可能收敛到不同的局域极小值，对应不同的 occupation matrix，从而得到不同的能量和物理性质。

Occupation Matrix Control（`omc`）功能允许用户在 SCF 迭代过程中固定 occupation matrix，可以帮助：

1. 将体系引导到特定的电子态（避免陷入非目标极小值）
2. 验证收敛后的结果是否稳定

### 5.2 三种模式

ABACUS 提供三种 `omc` 模式：

| `omc` 值 | 行为 |
|----------|------|
| 0 | 不使用 occupation matrix control，标准 DFT+U 计算（默认） |
| 1 | SCF 第一步读入 `initial_onsite.dm`，后续正常更新 occupation matrix |
| 2 | 读入 `initial_onsite.dm`，整个 SCF 过程中始终使用此 occupation matrix，不更新 |

### 5.3 使用流程

以验证 NiO 计算为例，使用 omc = 2：

**第一步：** 从上一步计算的 `OUT.NiO/onsite.dm` 复制到工作目录，改名为 `initial_onsite.dm`：

```bash
cp OUT.NiO/onsite.dm ./initial_onsite.dm
```

**第二步：** 在 INPUT 文件中添加 `omc 2`：

```
#Parameter DFT+U
dft_plus_u    1
orbital_corr  2 2 -1
hubbard_u     5.0 5.0 0.0
omc           2
```

**第三步：** 重新运行 ABACUS，在 `running_scf.log` 中可以确认整个 SCF 过程中 occupation matrix 保持不变。

**验证结果：** 使用 omc = 2 固定收敛后的 occupation matrix 重新计算，最终得到的总能量和磁矩与标准 SCF 完全一致，说明之前的 SCF 已收敛到稳定的基态。

### 5.4 initial_onsite.dm 文件格式

`initial_onsite.dm` 的格式与 `onsite.dm` 完全相同，也与日志中 `L(S)DA+U` 块的输出格式一致。用户可以手动编辑此文件，指定自定义的 occupation matrix，以引导体系进入特定的电子态。

### 5.5 进阶参数：Yukawa Potential 自动确定 U 值

除上述参数外，ABACUS 还支持通过 Yukawa 势在 SCF 中自动计算 U 值，无需手动设定：

| 参数 | 说明 |
|------|------|
| `yukawa_potential 1` | 开启 Yukawa 势方法，程序内自行计算 U 值 |
| `yukawa_lambda` | 手动指定 Yukawa 势的 screening length；不设则由电子密度 on-the-fly 计算 |

这种方法适合不确定 U 值的场合，理论细节见参考文献 [1]。
## 六、结语

过渡金属氧化物、稀土化合物等强关联体系是第一性原理计算中的常见挑战。ABACUS 在 LCAO 基组下实现了 DFT+U 功能，参数设置直观，通过 `dft_plus_u`、`orbital_corr`、`hubbard_u` 三个参数即可快速上手。

本教程以反铁磁 NiO 为例，展示了从输入文件准备到结果分析的完整流程，以及 occupation matrix control 功能的使用方式。使用过程中如有问题，欢迎在 [ABACUS GitHub](https://github.com/deepmodeling/abacus-develop) 提交 Issue。

---

**参考文献**

[1] Qu X, Xu P, Jiang H, He L, Ren X. DFT+U within the framework of linear combination of numerical atomic orbitals. *The Journal of Chemical Physics*. 2022;156(23):234104. [doi:10.1063/5.0090122](https://doi.org/10.1063/5.0090122)
