---
title: "使用 ABACUS + PYATB 计算材料光学性质：以晶态 SiO₂ 为例"
author: "AutoTutorial 3.0"
date: "2026-03-02"
topic: "光学性质（介电函数/吸收谱）计算"
task_type: "C"
has_case: true
case_material: "晶态 SiO₂（β-方英石），8 Si + 16 O"
word_count: 2122
lines: 660
chapters: 5
---

# 使用 ABACUS + PYATB 计算材料光学性质：以晶态 SiO₂ 为例

# 前言

本教程以晶态二氧化硅（β-方英石 SiO₂）为例，演示 **ABACUS + PYATB** 联合计算材料介电函数的完整流程。同时提供 Python 后处理脚本，展示如何从介电函数进一步得到折射率、消光系数、吸收系数等线性光学性质。

## 适用读者

- 熟悉 DFT 第一性原理计算基本流程
- 了解 ABACUS LCAO 基组计算的基本操作
- 希望计算材料光学性质的研究人员

## 前置知识

- ABACUS 基本输入文件（INPUT、STRU、KPT）的编写
- SCF（自洽场）计算的基本概念
- Python 基础（用于后处理）

## 教程结构

| 步骤 | 内容 | 工具 |
|------|------|------|
| 第一章 | 物理原理与计算方法 | — |
| 第二章 | 案例材料介绍 | — |
| 第三章 | 自洽计算，输出紧束缚矩阵 | ABACUS |
| 第四章 | 计算介电函数 | PYATB |
| 第五章 | 后处理，得到各光学性质 | Python |

## 计算环境

本教程使用以下环境：

- **镜像**：`registry.dp.tech/dptech/prod-19853/abacus-pyatb-open:v0.0.1`
- **推荐配置**：`c16_m32_cpu`（16 核 32 GB）

> **说明**：该镜像已预装 ABACUS 和 PYATB，可直接在 Bohrium 平台使用。
# 一、计算原理

## 1.1 介电函数与线性光学性质

介电函数 $\epsilon(\omega)$ 是复数函数，描述材料对电磁场的响应：

$$
\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)
$$

**实部 $\epsilon_1$** 与材料的折射、色散特性相关；**虚部 $\epsilon_2$** 与光的吸收相关，其起始上升位置对应材料的光学带隙。

由介电函数可以导出一系列线性光学性质：

| 物理量 | 符号 | 计算公式 |
|--------|------|----------|
| 折射率 | $n$ | $\displaystyle n = \sqrt{\frac{\sqrt{\epsilon_1^2+\epsilon_2^2}+\epsilon_1}{2}}$ |
| 消光系数 | $\kappa$ | $\displaystyle \kappa = \sqrt{\frac{\sqrt{\epsilon_1^2+\epsilon_2^2}-\epsilon_1}{2}}$ |
| 吸收系数 | $\alpha$ | $\displaystyle \alpha = \frac{\sqrt{2}\,\omega}{c}\sqrt{\sqrt{\epsilon_1^2+\epsilon_2^2}-\epsilon_1}$ |
| 能量损失函数 | $L$ | $\displaystyle L = \operatorname{Im}\!\left(\frac{-1}{\epsilon}\right) = \frac{\epsilon_2}{\epsilon_1^2+\epsilon_2^2}$ |
| 反射率 | $R$ | $\displaystyle R = \frac{(n-1)^2+\kappa^2}{(n+1)^2+\kappa^2}$ |

> **单位说明**：吸收系数 $\alpha$ 的单位为 cm⁻¹。将能量 $E$（eV）换算为角频率时，
> $\omega = E/\hbar$，其中 $\hbar \approx 4.136\times10^{-15}\ \text{eV·s}$；
> 光速取 $c = 3\times10^{10}\ \text{cm/s}$。

## 1.2 ABACUS + PYATB 计算方法

### 紧束缚模型的构建

ABACUS 采用数值原子轨道（NAO）作为基函数展开 Kohn-Sham 波函数。对于给定 $\mathbf{k}$ 点，Kohn-Sham 方程在 NAO 基下变为广义本征值问题：

$$
H(\mathbf{k})\,C_n(\mathbf{k}) = E_{n\mathbf{k}}\,S(\mathbf{k})\,C_n(\mathbf{k})
$$

其中哈密顿矩阵 $H(\mathbf{k})$ 和重叠矩阵 $S(\mathbf{k})$ 通过傅里叶变换由实空间矩阵得到：

$$
H_{\nu\mu}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\,H_{\nu\mu}(\mathbf{R}), \quad
S_{\nu\mu}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\,S_{\nu\mu}(\mathbf{R})
$$

ABACUS 完成自洽计算后，通过设置 `out_mat_hs2 = 1` 和 `out_mat_r = 1`，可以将以下实空间矩阵以稀疏格式输出：
$H_{\nu\mu}(\mathbf{R})$（哈密顿量）、$S_{\nu\mu}(\mathbf{R})$（重叠矩阵）以及偶极矩阵 $r_{\nu\mu,a}(\mathbf{R})$（$a=x,y,z$）。
PYATB 读取这三组矩阵，在任意 $\mathbf{k}$ 点重建哈密顿量，用于后续物性计算。

### Kubo-Greenwood 公式与介电函数

光电导率张量由 Kubo-Greenwood 公式给出：

$$
\sigma_{\alpha\beta}(\hbar\omega) = -\frac{ie^2\hbar}{NV_{\text{cell}}} \sum_{\mathbf{k}}\sum_{n,m}
\frac{f_{n\mathbf{k}}-f_{m\mathbf{k}}}{E_{n\mathbf{k}}-E_{m\mathbf{k}}}
\frac{\langle n\mathbf{k}|v_\alpha|m\mathbf{k}\rangle\langle m\mathbf{k}|v_\beta|n\mathbf{k}\rangle}
{\hbar\omega + E_{n\mathbf{k}} - E_{m\mathbf{k}} + i\eta}
$$

其中 $f_{n\mathbf{k}}$ 为 Fermi-Dirac 占据函数，$\eta$ 为展宽参数，速度矩阵元 $\langle n\mathbf{k}|v_\alpha|m\mathbf{k}\rangle$ 由偶极矩阵 $r_{\nu\mu,a}(\mathbf{k})$ 和本征矢量计算得到。

在 PYATB 中，介电函数虚部由下式直接计算（避免 $\omega=0$ 处的奇点处理问题）：

$$
\epsilon_2^{\alpha\beta}(\omega) = \frac{e^2\pi}{\epsilon_0\hbar}
\int\frac{d\mathbf{k}}{(2\pi)^3}\sum_{n,m}f_{nm}\,r^\alpha_{nm}r^\beta_{mn}\,\delta(\omega_{nm}-\omega)
$$

介电函数实部 $\epsilon_1$ 则通过 Kramers-Kronig 变换由 $\epsilon_2$ 求得：

$$
\epsilon_1^{\alpha\beta}(\omega) = \delta_{\alpha\beta} + \frac{2}{\pi}\mathcal{P}\int_0^\infty
\frac{\omega'\,\epsilon_2^{\alpha\beta}(\omega')}{\omega'^2-\omega^2}\,d\omega'
$$

整体计算流程如下：

```
ABACUS SCF 自洽计算
  │  (basis_type = lcao, out_mat_hs2 = 1, out_mat_r = 1)
  ↓
输出三组矩阵文件
  │  data-HR-sparse_SPIN0.csr  (哈密顿量)
  │  data-SR-sparse_SPIN0.csr  (重叠矩阵)
  │  data-rR-sparse.csr        (偶极矩阵)
  ↓
PYATB 计算介电函数
  │  (Kubo-Greenwood + Kramers-Kronig)
  ↓
输出 dielectric_function_*.dat
  ↓
Python 后处理
  └→ 折射率、消光系数、吸收系数、能量损失函数
```
# 二、案例介绍：晶态 SiO₂

本教程选用**β-方英石（β-cristobalite）SiO₂** 单胞作为示例材料。β-方英石是 SiO₂ 的一种高温相，具有立方对称结构（空间群 $Fd\bar{3}m$）。

SiO₂ 是典型的宽带隙绝缘体，光学透明窗口覆盖紫外到近红外波段。其光学性质在光学涂层、光纤、半导体器件等领域有重要应用价值。

## 晶体结构

| 参数 | 数值 |
|------|------|
| 晶系 | 立方 |
| 晶格常数 $a$ | 7.12 Å |
| 原胞中原子数 | 24（8 个 Si + 16 个 O） |
| 空间群 | $Fd\bar{3}m$ |

单胞中每个 Si 原子以四面体方式与 4 个 O 原子成键，形成 SiO₄ 四面体网络。

## 计算条件概述

| 项目 | 设置 |
|------|------|
| 计算软件 | ABACUS v3.4.3（SCF） + PYATB（介电函数） |
| 交换关联泛函 | GGA-PBE |
| 基组 | LCAO，数值原子轨道 2s2p1d（7 au 截断） |
| 赝势 | ONCV PBE |
| 截断能 | 100 Ry |
| k 点网格 | 6×6×6 Gamma 中心 |
| 介电函数 k 网格 | 20×20×20 |
| 能量范围 | 0 ~ 30 eV |

> **注意**：GGA-PBE 会低估带隙，得到的光学带隙偏低。若需要准确的带隙，应采用杂化泛函（如 HSE06）或 GW 方法，但计算成本更高。本案例重点演示计算流程。
# 三、ABACUS 自洽计算

本章完成 ABACUS 的 SCF 计算，目标是得到收敛的电子结构，并输出构建紧束缚模型所需的三组矩阵文件。

进入计算目录：

```bash
cd ~/abacus-pyatb_tutorial/silica_PrimaryCell
```

## 3.1 INPUT 文件

INPUT 文件控制所有计算参数。对于 PYATB 接口，有两个参数是必须的：`out_mat_hs2`（输出哈密顿量和重叠矩阵）和 `out_mat_r`（输出偶极矩阵）；`out_chg` 为可选项，输出电荷密度文件。

```
# INPUT
INPUT_PARAMETERS

suffix                  silica
calculation             scf
esolver_type            ksdft
symmetry                0           # 关闭对称性，确保矩阵文件完整输出
init_chg                atomic

pseudo_dir              ./
orbital_dir             ./

basis_type              lcao        # 必须使用 LCAO 基组
ks_solver               genelpa     # 适合 LCAO 的并行本征值求解器
smearing_method         gaussian
smearing_sigma          0.01        # Ry，绝缘体可适当减小
mixing_type             broyden
mixing_beta             0.1
ecutwfc                 100         # Ry，截断能
scf_thr                 1e-8        # Ry，SCF 收敛阈值
mixing_gg0              1.5
mixing_ndim             20

out_chg                 1           # 输出自洽后的电荷密度
out_mat_hs2             1           # 输出实空间 HR、SR 矩阵（PYATB 必需）
out_mat_r               1           # 输出实空间偶极矩阵 rR（PYATB 必需）
```

### 关键参数说明

| 参数 | 本案例值 | 作用 |
|------|----------|------|
| `basis_type` | `lcao` | PYATB 要求使用 LCAO 基组 |
| `symmetry` | `0` | 关闭对称性，避免矩阵输出不完整 |
| `out_mat_hs2` | `1` | 输出 `data-HR-sparse_SPIN0.csr`（哈密顿量）和 `data-SR-sparse_SPIN0.csr`（重叠矩阵） |
| `out_mat_r` | `1` | 输出 `data-rR-sparse.csr`（位置偶极矩阵，计算速度矩阵元必需） |
| `ecutwfc` | `100` Ry | 控制实空间格点密度，影响矩阵精度 |
| `scf_thr` | `1e-8` | 收敛阈值，光学计算对收敛精度要求较高 |
| `mixing_gg0` | `1.5` | 倒空间混合修正，帮助 SiO₂ 这类绝缘体收敛 |

## 3.2 STRU 文件

STRU 文件包含晶体结构信息。本案例为立方 SiO₂ 超胞，晶格常数 7.12 Å，共 8 个 Si 原子和 16 个 O 原子。

```
# STRU
ATOMIC_SPECIES
Si 28.086  Si_ONCV_PBE-1.0.upf
O 15.999   O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_7au_100Ry_2s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261246257702          # 1 Å = 1.8897 Bohr，此值使后续坐标单位为 Å

LATTICE_VECTORS
7.1199998856 0.0 0.0        # 单位：LATTICE_CONSTANT（即 Å）
0.0 7.1199998856 0.0
0.0 0.0 7.1199998856

ATOMIC_POSITIONS
Cartesian                   # 笛卡尔坐标，单位：LATTICE_CONSTANT（Å）
Si
0.0                         # 初始磁矩
8                           # Si 原子数
0.000000000000 0.000000000000 0.000000000000 1 1 1
0.000000000000 3.559999943000 3.559999943000 1 1 1
3.559999943000 3.559999943000 0.000000000000 1 1 1
3.559999943000 0.000000000000 3.559999943000 1 1 1
5.339999914000 1.779999971000 5.339999914000 1 1 1
1.779999971000 1.779999971000 1.779999971000 1 1 1
1.779999971000 5.339999914000 5.339999914000 1 1 1
5.339999914000 5.339999914000 1.779999971000 1 1 1
O
0.0                         # 初始磁矩
16                          # O 原子数
0.889999986000 0.889999986000 0.889999986000 1 1 1
6.229999900000 2.669999957000 4.449999928000 1 1 1
2.669999957000 4.449999928000 6.229999900000 1 1 1
4.449999928000 6.229999900000 2.669999957000 1 1 1
0.889999986000 4.449999928000 4.449999928000 1 1 1
6.229999900000 6.229999900000 0.889999986000 1 1 1
2.669999957000 0.889999986000 2.669999957000 1 1 1
4.449999928000 2.669999957000 6.229999900000 1 1 1
4.449999928000 0.889999986000 4.449999928000 1 1 1
2.669999957000 2.669999957000 0.889999986000 1 1 1
6.229999900000 4.449999928000 2.669999957000 1 1 1
0.889999986000 6.229999900000 6.229999900000 1 1 1
4.449999928000 4.449999928000 0.889999986000 1 1 1
2.669999957000 6.229999900000 4.449999928000 1 1 1
6.229999900000 0.889999986000 6.229999900000 1 1 1
0.889999986000 2.669999957000 2.669999957000 1 1 1
```

### 轨道基组说明

本案例选用 `2s2p1d-7au` 轨道，命名规则为 `元素_泛函_截断半径_截断能_轨道类型.orb`：

- `Si_gga_7au_100Ry_2s2p1d.orb`：Si 元素，GGA 泛函，7 au 截断半径，100 Ry 截断能，2 个 s 轨道 + 2 个 p 轨道 + 1 个 d 轨道（共 13 个基函数）
- `O_gga_7au_100Ry_2s2p1d.orb`：O 元素，参数同上（共 13 个基函数）

24 个原子 × 13 个基函数 = **312 个 NAO 基函数**（与 SCF 输出中 `NBASE = 312` 一致）。

## 3.3 KPT 文件

```
# KPT
K_POINTS
0
Gamma
6 6 6 0 0 0
```

采用 6×6×6 Gamma 中心 Monkhorst-Pack 网格，对应 112 个不等价 k 点。这是 SCF 计算的 k 网格；介电函数计算时 PYATB 会单独使用更密的 k 网格（20×20×20）。

## 3.4 运行与检查输出

提交计算：

```bash
export OMP_NUM_THREADS=1 && mpirun -np 16 abacus
```

### SCF 收敛过程

正常收敛的输出如下（共 13 步收敛到 `1e-8` 以下）：

```
ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)
GE1    -7.851283e+03  0.000000e+00   2.125e-01  1.011e+01
GE2    -7.826465e+03  2.481833e+01   1.691e-01  9.909e+00
...
GE13   -7.835176e+03  -1.392110e-11  7.527e-10  9.462e+00
```

### 从输出中读取关键信息

SCF 完成后，需要从日志中读取两个参数，用于后续 PYATB Input 文件：

```bash
grep 'occupied bands' OUT.silica/running_scf.log
grep 'E_Fermi' OUT.silica/running_scf.log | tail -1
```

```
occupied bands = 64
E_Fermi    0.4070749075    5.5385382545
```

- **occupied bands = 64**：SiO₂ 单胞有 24 个原子（Si 有 4 个价电子，O 有 6 个），总价电子数 = 8×4 + 16×6 = 128，nspin=1 时占据 64 条能带。
- **E_Fermi = 5.5385382545 eV**：第二列为 eV 单位，**直接填入 PYATB Input 文件**。

### 输出文件列表

SCF 完成后，`OUT.silica/` 目录中包含 PYATB 所需的三个矩阵文件：

```
OUT.silica/
├── data-HR-sparse_SPIN0.csr    # 实空间哈密顿量矩阵（稀疏 CSR 格式）
├── data-SR-sparse_SPIN0.csr    # 重叠矩阵
├── data-rR-sparse.csr          # 位置偶极矩阵（x、y、z 分量）
├── running_scf.log             # 计算日志
└── ...
```
# 四、PYATB 介电函数计算

本章使用 PYATB 读取 ABACUS 输出的紧束缚矩阵，在密 k 网格上计算介电函数张量。

## 4.1 准备矩阵文件

将 ABACUS 输出的三个矩阵文件复制到 PYATB 工作目录：

```bash
cp OUT.silica/data* ./pyatb_OpticConductivity/
cd pyatb_OpticConductivity
```

目录结构：

```
pyatb_OpticConductivity/
├── Input                        # PYATB 控制文件（需手动填写）
├── data-HR-sparse_SPIN0.csr     # 从 ABACUS 复制
├── data-SR-sparse_SPIN0.csr     # 从 ABACUS 复制
└── data-rR-sparse.csr           # 从 ABACUS 复制
```

## 4.2 Input 文件详解

PYATB 的控制文件固定命名为 `Input`，由三个区块组成：

```
# Input
INPUT_PARAMETERS
{
nspin               1
package             ABACUS
fermi_energy        5.5385382545    # eV，从 SCF 日志中读取（E_Fermi 第二列）
fermi_energy_unit   eV
HR_route            data-HR-sparse_SPIN0.csr
SR_route            data-SR-sparse_SPIN0.csr
rR_route            data-rR-sparse.csr
HR_unit             Ry              # ABACUS 输出的哈密顿量单位为 Rydberg
rR_unit             Bohr            # 偶极矩阵单位为 Bohr
}

LATTICE
{
lattice_constant        1.8897261246257702
lattice_constant_unit   Bohr
lattice_vector
7.1199998856 0.0 0.0
0.0 7.1199998856 0.0
0.0 0.0 7.1199998856
}

OPTICAL_CONDUCTIVITY
{
occ_band    64      # 占据能带数，从 SCF 日志 "occupied bands" 读取
omega       0 30    # 计算能量范围：0 ~ 30 eV
domega      0.01    # 能量步长：0.01 eV，共 3000 个点
eta         0.1     # 展宽参数（eV），模拟有限寿命效应
grid        20 20 20  # 布里渊区积分网格
}
```

### 参数说明

**INPUT_PARAMETERS 区块**

| 参数 | 本案例值 | 说明 |
|------|----------|------|
| `nspin` | `1` | 与 ABACUS 的 nspin 保持一致 |
| `package` | `ABACUS` | 指定矩阵文件来源为 ABACUS |
| `fermi_energy` | `5.5385382545` | 费米能（eV），从 SCF 日志最后一行 E_Fermi 读取 |
| `HR_unit` | `Ry` | ABACUS 输出的 HR 单位固定为 Rydberg |
| `rR_unit` | `Bohr` | 偶极矩阵单位固定为 Bohr |

**LATTICE 区块**

晶格参数必须与 ABACUS STRU 文件完全一致：
- `lattice_constant`：取 STRU 中的值 `1.8897261246257702`（对应 1 Å 的 Bohr 换算因子）
- `lattice_vector`：与 STRU 中的 LATTICE_VECTORS 相同

**OPTICAL_CONDUCTIVITY 区块**

| 参数 | 本案例值 | 说明 |
|------|----------|------|
| `occ_band` | `64` | 占据能带数，需与 SCF 结果一致 |
| `omega` | `0 30` | 能量范围（eV）。上限覆盖感兴趣的光学跃迁 |
| `domega` | `0.01` | 能量步长（eV）。越小谱线越细，计算量正比增加 |
| `eta` | `0.1` | 展宽参数（eV），对应 Lorentz 展宽。越小峰形越尖锐，但需要更密 k 网格；太小会出现噪声 |
| `grid` | `20 20 20` | k 积分网格。光学性质通常需要比 SCF 更密的网格；对本案例 20³ 已足够，复杂材料可能需要 40³ 以上 |

> **关于 `eta` 的选择**：`eta` 同时影响谱线宽度和 k 网格需求。减小 `eta` 需同步加密 `grid`，否则积分不收敛，谱线出现虚假震荡。

## 4.3 运行与输出

```bash
export OMP_NUM_THREADS=1 && mpirun -np 16 pyatb
```

计算过程中不会在屏幕上输出进度信息，可以追踪日志文件：

```bash
tail -f Out/running.log
```

计算完成后，进入输出目录：

```bash
cd Out/Optical_Conductivity
ls
```

```
dielectric_function_imag_part.dat    # 介电函数虚部 ε₂
dielectric_function_real_part.dat    # 介电函数实部 ε₁
optical_conductivity.dat             # 光电导率
```

### 输出文件格式

以 `dielectric_function_imag_part.dat` 为例，文件格式如下：

```
# omega(eV)  xx    xy    xz    yx    yy    yz    zx    zy    zz
  0.00000   0.000000e+00  -6.324391e-15  ...
  0.01000   9.564624e-06  -6.324391e-15  ...
  0.02000   1.912933e-05  ...
```

- 第一列：光子能量 $\omega$（eV）
- 后 9 列：介电函数张量 $\epsilon_2^{\alpha\beta}$（$\alpha,\beta=x,y,z$）的各分量，顺序为 xx、xy、xz、yx、yy、yz、zx、zy、zz

对于立方 SiO₂，三个对角分量 xx = yy = zz，非对角分量接近 0（~10⁻¹⁴ 量级），体现材料的各向同性。

低能量端（< 光学带隙）的 $\epsilon_2$ 应接近 0；$\epsilon_1(0)$（零频极限）应接近材料的静态介电常数，本案例约为 1.90（参考下一章实部数据）。

> **验证**：SCF 收敛精度不足时，$\epsilon_2$ 在低频段可能出现非零"尾巴"，需检查 `scf_thr` 设置。
# 五、后处理：从介电函数到光学性质

PYATB 输出了介电函数张量的全部 9 个分量。本章用 Python 读取数据，先可视化介电函数，再计算其余光学性质。

## 5.1 介电函数的读取与可视化

以下脚本读取虚部和实部数据，分别绘制各分量随能量的变化：

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

RealPartData = 'Out/Optical_Conductivity/dielectric_function_real_part.dat'
ImagPartData = 'Out/Optical_Conductivity/dielectric_function_imag_part.dat'

# 读取数据
with open(ImagPartData, 'r') as f:
    labels_imag = f.readline().strip().split()[1:]
data_imag = np.loadtxt(ImagPartData, skiprows=1)
omega = data_imag[:, 0]
imag_parts = data_imag[:, 1:]

with open(RealPartData, 'r') as f:
    labels_real = f.readline().strip().split()[1:]
data_real = np.loadtxt(RealPartData, skiprows=1)
real_parts = data_real[:, 1:]

# 绘图
markers = ['o', '*', '^', 'd', 's', 'p', '+', 'x', 'h']
fig = plt.figure(constrained_layout=True)
gs = GridSpec(1, 2, figure=fig)

ax_imag = fig.add_subplot(gs[0, 0])
for i, (imag_part, marker) in enumerate(zip(imag_parts.T, markers), start=1):
    ax_imag.plot(omega, imag_part, linestyle='-', marker=marker,
                 markersize=5, label=labels_imag[i])
ax_imag.legend()
ax_imag.set_title('Dielectric Function Imaginary Part')
ax_imag.set_xlabel('Energy (eV)')
ax_imag.set_ylabel(r'$\epsilon_2$')

ax_real = fig.add_subplot(gs[0, 1])
for i, (real_part, marker) in enumerate(zip(real_parts.T, markers), start=1):
    ax_real.plot(omega, real_part, linestyle='-', marker=marker,
                 markersize=5, fillstyle='none', mew=0.4, label=labels_real[i])
ax_real.legend()
ax_real.set_title('Dielectric Function Real Part')
ax_real.set_xlabel('Energy (eV)')
ax_real.set_ylabel(r'$\epsilon_1$')

plt.show()
```

由于 SiO₂ 是各向同性立方结构，三个对角分量（xx、yy、zz）完全重合，非对角分量为零。

实部 $\epsilon_1(0) \approx 1.90$，对应折射率 $n(0) \approx 1.38$。虚部 $\epsilon_2$ 在约 10 eV 以下（带隙以下）接近 0，之后随能量上升出现明显吸收峰。

## 5.2 各光学性质的计算

取各向同性情况下的平均值（(xx+yy+zz)/3）进行后续计算：

```python
def read_and_average(filename):
    """读取介电函数文件，返回能量轴和各向同性平均值"""
    data = np.loadtxt(filename, skiprows=1)
    energy = data[:, 0]
    # 取 xx(列1)、yy(列5)、zz(列9) 的平均值
    xx = data[:, 1]
    yy = data[:, 5]
    zz = data[:, 9]
    return energy, (xx + yy + zz) / 3

energy, eps2 = read_and_average(ImagPartData)    # 虚部
_,     eps1 = read_and_average(RealPartData)     # 实部

# ── 计算各光学量 ──────────────────────────────────────────
eps_abs = np.sqrt(eps1**2 + eps2**2)             # |ε|

# 折射率
n = np.sqrt((eps_abs + eps1) / 2)

# 消光系数
kappa = np.sqrt((eps_abs - eps1) / 2)

# 吸收系数（单位：cm⁻¹）
# ω = E/ħ，ħ ≈ 4.136e-15 eV·s，c = 3e10 cm/s
hbar = 4.135667696e-15   # eV·s
c    = 3.0e10            # cm/s（注意：此处需用 cm/s 而非 m/s，以得到 cm⁻¹ 单位）
omega_rad = energy / hbar                         # s⁻¹
alpha = np.sqrt(2) * omega_rad / c * np.sqrt(eps_abs - eps1)

# 能量损失函数
L = eps2 / (eps1**2 + eps2**2)

# ── 绘图 ──────────────────────────────────────────────────
fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))

axs[0, 0].plot(energy, n, 'b-')
axs[0, 0].set_ylabel('Refraction Index $n$')
axs[0, 0].set_xlabel('Energy (eV)')
axs[0, 0].grid(True)

axs[0, 1].plot(energy, kappa, 'm-')
axs[0, 1].set_ylabel('Extinction Coefficient $\\kappa$')
axs[0, 1].set_xlabel('Energy (eV)')
axs[0, 1].grid(True)

axs[1, 0].plot(energy, alpha, 'y-')
axs[1, 0].set_ylabel(r'Absorption Coefficient (cm$^{-1}$)')
axs[1, 0].set_xlabel('Energy (eV)')
axs[1, 0].set_yscale('log')
axs[1, 0].grid(True)

axs[1, 1].plot(energy, L, 'r-')
axs[1, 1].set_ylabel('Energy Loss Function $L$')
axs[1, 1].set_xlabel('Energy (eV)')
axs[1, 1].grid(True)

plt.suptitle('Optical Properties of SiO₂ vs Energy')
plt.tight_layout()
plt.show()
```

> **注意**：代码中各向同性平均的列索引需与输出文件列顺序对应：
> 列 1 = xx，列 5 = yy，列 9 = zz（从 0 开始计数）。

## 5.3 结果物理解读

**介电函数**

- $\epsilon_1(0) \approx 1.90$，对应静态介电常数，与 SiO₂ 的实验值（~2.1）接近（PBE 低估带隙导致轻微偏差）
- $\epsilon_2$ 在约 9 eV 以下接近零，反映 SiO₂ 宽带隙绝缘体特征（实验带隙 ~8.9 eV，PBE 计算约低估 1~2 eV）
- $\epsilon_1$ 在 $\epsilon_2$ 峰值附近出现过零点，对应强吸收区

**折射率与消光系数**

- 低能区（可见光及近紫外）折射率约 1.38，接近实验值（~1.46，偏差来自 PBE 低估带隙）
- 高能区（>10 eV，深紫外）消光系数 $\kappa$ 明显上升，对应强烈光吸收

**吸收系数**

- 吸收系数以对数坐标显示，在带隙以下接近 0（通常 < 10² cm⁻¹）
- 超过光学带隙后急剧上升，可达 10⁶ cm⁻¹ 量级
- 吸收边（$\alpha$ 开始上升的能量）可用于估读 DFT 计算的光学带隙

**能量损失函数**

- $L(\omega)$ 的峰值对应体等离激元频率（plasmon resonance）
- SiO₂ 的等离激元峰位于 ~20 eV 左右，是其典型特征

> **提高精度的建议**：
> - 增大 PYATB 的 `grid`（如 40×40×40）可改善 k 积分收敛性
> - 减小 `eta`（如 0.05 eV）可获得更尖锐的谱线，但需配合更密的 `grid`
> - 若需准确带隙，改用 HSE06 泛函或 G₀W₀ 修正
# 附录

## 参考资料

1. PYATB 官方文档（光学电导率模块）：https://pyatb.github.io/pyatb/functions/optical_conductivity.html
2. PYATB GitHub 仓库：https://github.com/pyatb/pyatb
3. ABACUS 官方文档：https://abacus.deepmodeling.com
4. 光学电导率 Wikipedia：https://en.wikipedia.org/wiki/Optical_conductivity
5. 介电函数计算技术笔记：https://www.openmx-square.org/tech_notes/Dielectric_Function_YTL.pdf

**PYATB 引用文献**：
> Gan Jin, Hongsheng Pang, Yuyang Ji, Zujian Dai, and Lixin He, *PYATB: An efficient Python package for electronic structure calculations using ab initio tight-binding model*, Comput. Phys. Commun. **291**, 108844 (2023).

## 常见错误与排查

| 错误现象 | 可能原因 | 解决方法 |
|----------|----------|----------|
| `data-rR-sparse.csr` 文件不存在 | ABACUS 未设置 `out_mat_r = 1` | 在 INPUT 中添加 `out_mat_r 1` 后重新计算 |
| `data-HR-sparse_SPIN0.csr` 文件不存在 | ABACUS 未设置 `out_mat_hs2 = 1` | 在 INPUT 中添加 `out_mat_hs2 1` 后重新计算 |
| PYATB 报错：费米能不匹配 | `fermi_energy` 填写了 Ry 单位的值 | 从 running_scf.log 读取 eV 值（E_Fermi 第二列） |
| 矩阵文件维度不一致 | ABACUS 对称性设置导致 k 点不完整 | 在 ABACUS INPUT 中设置 `symmetry 0` |
| $\epsilon_2$ 在低能区非零 | SCF 未收敛 | 收紧 `scf_thr`（如 `1e-9`）重新计算 |
| PYATB 运行无输出 | 正常现象 | 追踪 `Out/running.log` 查看进度 |

## 进阶学习方向

- **杂化泛函（HSE06）**：得到更准确的带隙和光学谱，代价是计算时间增加约 10 倍
- **非线性光学性质**：PYATB 还支持二阶非线性光学响应（SHG、二阶霍尔效应）的计算
- **k 点收敛性测试**：对于金属或带隙很小的材料，`grid` 需要仔细测试（50³ 甚至更密）
- **自旋极化体系**：设置 `nspin = 2` 时，需分别处理两个自旋通道的 HR 矩阵
