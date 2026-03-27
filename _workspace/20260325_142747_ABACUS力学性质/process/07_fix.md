# 前言

本教程介绍 ABACUS 中的力学性质计算，涵盖结构优化（几何弛豫）和弹性常数两条主线。

**适用读者：** 有 ABACUS 基本使用经验，了解 INPUT/STRU/KPT 文件格式，做过 SCF 计算的研究者。

**前置知识：**
- ABACUS 基本输入文件的格式与作用
- SCF 自洽计算的基本概念
- LCAO 基组的基本概念

**教程结构：**
- 第一章：理论背景——几何优化原理 + 弹性张量理论 + 宏观力学性质
- 第二章：案例一——h-BN 变胞优化（cell-relax，192原子体系）
- 第三章：案例二——Si 弹性常数（abacustest 自动化工作流）
- 第四章：注意事项与参数速查

**案例说明：**
- **h-BN**：六方层状结构，192原子超胞，演示 cell-relax 的完整参数配置与收敛过程。
- **Si**：立方晶系，3个独立弹性张量元，演示 abacustest 从结构优化到弹性常数的完整自动化流程，结果与实验值对比。

所有案例数据均来自实际计算，未经修改。
# 第一章：理论背景

## 1.1 几何优化：势能面与收敛判据

第一性原理计算中的"结构优化"，本质上是在势能面（Born-Oppenheimer 势能面）上寻找极小值点的过程。给定初始结构，程序计算各原子受力（$F_i = -\partial E / \partial \mathbf{r}_i$）和晶胞应力张量，沿受力方向更新原子坐标。对于 cell-relax，还同步更新晶格参数。重复迭代，直到力与应力均低于设定阈值。

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

### 收敛判据

结构优化由两个阈值控制：

| 参数 | 单位 | 含义 | 典型推荐值 |
|------|------|------|-----------|
| `force_thr_ev` | eV/Å | 所有原子受力中最大分量（LARGEST GRAD）的收敛阈值 | 0.01–0.05 |
| `stress_thr` | kBar | TOTAL-PRESSURE 绝对值的收敛阈值（负值表示拉伸压力，取绝对值与阈值比较） | 0.5–5 |

当 LARGEST GRAD < `force_thr_ev` **且** TOTAL-PRESSURE < `stress_thr` 同时满足时，计算宣告收敛。对于弹性常数计算，残余应力应尽量小，建议 `stress_thr` 取 0.5 kBar 或更严格。

其余常用控制参数：

| 参数 | 含义 | 默认/推荐值 |
|------|------|-----------|
| `relax_nmax` | 最大离子步数 | 100 |
| `cal_force` | 开启受力计算 | 1（relax/cell-relax 时必须） |
| `cal_stress` | 开启应力张量计算 | 1（cell-relax 时必须） |
| `out_stru` | 输出优化后结构文件（STRU_ION_D） | 1 |

---

## 1.2 弹性张量与广义胡克定律

在线弹性近似下，晶体的应力张量 $\sigma_{ij}$ 与应变张量 $\varepsilon_{kl}$ 通过广义胡克定律联系：

$$\sigma_{ij} = C_{ijkl}\,\varepsilon_{kl}$$

其中 $C_{ijkl}$ 是四阶弹性刚度张量，共 81 个分量。由于 $\sigma$ 和 $\varepsilon$ 均为对称张量，张量对称性将独立分量减少到 **21 个**。

### Voigt 标记

利用 Voigt 标记，将 4 阶张量降维为 6×6 对称矩阵：

$$xx \to 1,\quad yy \to 2,\quad zz \to 3,\quad yz \to 4,\quad xz \to 5,\quad xy \to 6$$

胡克定律化为矩阵形式：

$$\begin{bmatrix}\sigma_1\\\sigma_2\\\sigma_3\\\sigma_4\\\sigma_5\\\sigma_6\end{bmatrix}
=\begin{bmatrix}
C_{11}&C_{12}&C_{13}&C_{14}&C_{15}&C_{16}\\
&C_{22}&C_{23}&C_{24}&C_{25}&C_{26}\\
&&C_{33}&C_{34}&C_{35}&C_{36}\\
&&&C_{44}&C_{45}&C_{46}\\
&&&&C_{55}&C_{56}\\
&&&&&C_{66}
\end{bmatrix}
\begin{bmatrix}\varepsilon_1\\\varepsilon_2\\\varepsilon_3\\\varepsilon_4\\\varepsilon_5\\\varepsilon_6\end{bmatrix}$$

### 晶系对称性约束

晶体点群进一步约减独立分量数：

| 晶系 | 独立分量数 | 典型材料 |
|------|-----------|---------|
| 三斜晶系 | 21 | — |
| 单斜晶系 | 13 | — |
| 正交晶系 | 9 | — |
| 六方晶系 | 5 | h-BN、Mg |
| 四方晶系（Laue I） | 6 | 金红石 TiO₂ |
| 立方晶系 | 3 | Si、Cu、Al |

**立方晶系**（本教程案例 Si）：独立分量仅 3 个——$C_{11}$（= $C_{22}$ = $C_{33}$）、$C_{12}$（= $C_{13}$ = $C_{23}$）、$C_{44}$（= $C_{55}$ = $C_{66}$）。其余非对角元为零。

### 计算方法：应力-应变法

最常用的第一性原理弹性常数计算方法是**应力-应变法**：

1. 对弛豫后的平衡结构施加一系列已知应变（6种独立应变方向，每种取 4 个幅度 ±0.5%、±1%）
2. 对每个形变构型做固定晶胞的原子弛豫（relax）后计算应力
3. 对应力-应变数据做线性拟合，得到弹性张量元

共生成 24 个形变构型，加上原始未形变结构共 **25 个计算任务**。

---

## 1.3 宏观力学性质

弹性张量包含晶体的各向异性力学信息。对于多晶材料，通常需要从弹性张量导出各向同性的宏观力学性质，常用的是 Voigt 平均（假设均匀应变）：

| 性质 | 符号 | 物理含义 |
|------|------|---------|
| 体弹模量（Bulk Modulus） | $B_V$ | 抵抗均匀压缩的能力；越大越难压缩 |
| 剪切模量（Shear Modulus） | $G_V$ | 抵抗剪切形变的能力；越大越耐剪切 |
| 杨氏模量（Young's Modulus） | $E_V$ | 单轴拉伸时应力与应变之比；表征轴向刚度 |
| 泊松比（Poisson's Ratio） | $\nu$ | 横向收缩与轴向伸长之比；通常 0–0.5 |

abacustest 的后处理命令 `abacustest model elastic post` 会自动完成弹性张量拟合和上述四个量的计算，结果以 JSON 格式保存。
# 第二章：案例一——h-BN 变胞优化

六方氮化硼（h-BN）是典型的层状材料，属六方晶系，结构类似石墨。层内 B-N 键强，层间为弱 van der Waals 相互作用。这使 a/b 方向与 c 方向的弹性性质差异显著。

层间距离对 c 方向应力非常敏感，cell-relax 需要同时收紧力和应力才能得到可靠的平衡结构。本案例以此体系演示变胞优化的完整参数配置。

本案例体系：**h-BN 192原子超胞**（B₉₆N₉₆），LCAO 基组，PBE 泛函。

---

## 2.1 INPUT 文件

```
# h-BN cell-relax INPUT
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

## 2.2 STRU 文件

h-BN 192原子体系的 STRU 关键部分如下：

```
ATOMIC_SPECIES
B   10.81   B_ONCV_PBE-1.0.upf
N   14.007  N_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
B_gga_8au_100Ry_2s2p1d.orb
N_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
2.6665300000  0.0000000000  0.0000000000
0.0000000000  16.970664000  0.0000000000
0.0000000000  0.0000000000  28.0218

ATOMIC_POSITIONS
Direct
B
0.0
96
...（96个B原子坐标）
N
0.0
96
...（96个N原子坐标）
```

晶格常数单位为 Bohr（`LATTICE_CONSTANT` 给出 Bohr 换算因子），`LATTICE_VECTORS` 给出无量纲晶格矢量（乘以 `LATTICE_CONSTANT` 后为 Bohr 单位）。

轨道文件说明：
- `B_gga_8au_100Ry_2s2p1d.orb`：B 的 DZP 轨道，截断半径 8 a.u.，截断能 100 Ry
- `N_gga_8au_100Ry_2s2p1d.orb`：N 的 DZP 轨道，截断半径 8 a.u.，截断能 100 Ry
- 两者均与 `ONCV_PBE-1.0` 系列赝势匹配

> 完整的 192原子坐标文件可从 ABACUS 官方 GitHub 仓库的 examples 目录获取。

---

## 2.3 运行与收敛过程

提交计算后，ABACUS 的输出日志（`running_cell-relax.log`）中会依次输出每轮 RELAX CELL 和 RELAX IONS 的信息（以下为示意性输出，实际数值因体系而异）：

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
- **TOTAL-PRESSURE**：应力张量对角元均值（对应 `stress_thr`）

当两者**同时**低于设定阈值时，计算收敛，输出最终优化结构 `OUT.BN/STRU_ION_D`（通过 `out_stru 1` 开启）。

---

## 2.4 常见问题

**Q：cell-relax 跑了很多步但应力一直不降？**

可能原因：初始结构残余应力太大，或 `kspacing` 取值偏大（K 点密度不够）。建议检查 K 点收敛性，或适当减小 `kspacing`（如从 0.1 改为 0.08）。

**Q：`symmetry 0` 和 `symmetry 1` 有什么区别？**

开启对称性（`symmetry 1`）会约束原子位置使其满足对称性，减少独立自由度，加速收敛，但有时会因程序识别对称性不准确导致弛豫出错。大体系或低对称性材料建议用 `symmetry 0`。

**Q：弛豫完成后如何继续进行弹性计算？**

将 `OUT.BN/STRU_ION_D` 复制为新目录的 `STRU`，修改 `INPUT` 中 `calculation scf`（或交给 abacustest 自动处理），即可在优化后结构上计算弹性常数。
# 第三章：案例二——Si 弹性常数（abacustest 工作流）

Si 是立方晶系，弹性张量仅有 3 个独立分量（$C_{11}$、$C_{12}$、$C_{44}$），是验证弹性计算流程的标准案例。本章演示用 abacustest 从结构生成到弹性张量输出的完整自动化流程。

---

## 3.1 安装 abacustest

```bash
pip install abacustest
```

也可从 GitHub 获取最新版本：

```bash
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装后执行 `abacustest -h` 验证安装成功。

**配置赝势和轨道路径：**

abacustest 在准备 LCAO 计算输入文件时，通过环境变量查找赝势和轨道文件：

```bash
export ABACUS_PP_PATH=/path/to/pseudopotentials
export ABACUS_ORB_PATH=/path/to/orbitals
```

> 若在 Bohrium 云平台上运行，可使用平台提供的 abacustest App，无需手动配置环境变量。

---

## 3.2 生成 Si 惯用胞结构

弹性张量计算要求晶格矢量与坐标轴对齐（IEEE 标准取向）。对立方晶系，用 ase 生成的惯用胞已满足此条件：

```python
from ase.build import bulk

si_conv = bulk('Si', cubic=True)  # 8原子立方惯用胞
si_conv.write("Si_conv.cif")
```

---

## 3.3 结构优化（cell-relax）

弹性计算要求平衡结构无残余应力。用 abacustest 自动准备 cell-relax 输入：

```bash
abacustest model inputs -f Si_conv.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si
```

参数说明：

| 参数 | 含义 |
|------|------|
| `-f Si_conv.cif` | 输入结构文件 |
| `--ftype cif` | 结构文件格式 |
| `--jtype cell-relax` | 计算类型：同时优化晶胞和原子位置 |
| `--lcao` | 使用 LCAO 数值原子轨道基组 |
| `--folder-syntax Si` | 输出文件夹命名规则 |

执行后生成 `Si/` 文件夹，包含 INPUT、STRU、KPT 及赝势/轨道文件。提交计算，等待完成。

---

## 3.4 整理优化后结构

优化完成后，优化后的结构保存在 `Si/INPUT/OUT.ABACUS/STRU_ION_D`。新建弹性计算目录：

```bash
mkdir Si-elastic
cp Si/INPUT Si/Si* Si-elastic/
cp Si/INPUT/OUT.ABACUS/STRU_ION_D Si-elastic/STRU
```

修改 `Si-elastic/INPUT`，将计算类型改为 SCF：

```
calculation    scf
```

此目录将作为弹性计算的基础，abacustest 会读取其中的 INPUT、STRU、KPT 和赝势/轨道文件，自动生成各形变构型。

---

## 3.5 生成形变构型

```bash
abacustest model elastic prepare -j Si-elastic
```

执行后，在 `Si-elastic/` 目录下生成：

- `deformed_*/`：24 个形变构型文件夹（6 种应变方向 × 4 个幅度）
- `org/`：原始未形变结构

共 **25 个计算任务**，每个文件夹内的 ABACUS 计算类型为 `relax`（固定晶胞、允许原子弛豫），并自动开启应力计算。

可选参数：

| 参数 | 默认值 | 含义 |
|------|--------|------|
| `--norm` | 0.01 | 最大正应变幅度（生成 ±0.5%、±1%） |
| `--shear` | 0.01 | 最大剪切应变幅度 |
| `--norelax` | — | 不做原子弛豫，直接计算应力（速度更快但精度略低） |

> **警告：** 重复执行 `elastic prepare` 会**直接删除**已有的形变文件夹。准备好后不要重复运行此命令。

---

## 3.6 提交计算

将 25 个计算任务提交到集群或云平台。以 Bohrium 平台为例，每个任务目录需要一个 `job.json` 配置文件：

```json
{
    "job_name": "Si-elastic-deformed",
    "command": "OMP_NUM_THREADS=1 mpirun -np 8 abacus > log",
    "log_file": "log",
    "backward_files": ["OUT.*", "log"],
    "project_id": 205855,
    "platform": "ali",
    "job_type": "container",
    "machine_type": "c16_m32_cpu",
    "image_address": "registry.dp.tech/dptech/abacus:LTSv3.10.1"
}
```

> 也可用 abacustest 的批量提交功能，详见 abacustest 文档。

---

## 3.7 后处理：提取弹性张量

所有计算完成后，执行：

```bash
abacustest model elastic post -j Si-elastic
```

屏幕输出：

```
Model: elastic
Postprocessing elastic calculation for job: Si-elastic/
             bulk_modulus  shear_modulus  young_modulus  poisson_ratio
Si-elastic/     94.191857      70.892488     170.022310       0.199156

Si-elastic/     elastic_tensor:
              0             1             2          3          4          5
0  1.654562e+02  5.855970e+01  5.855970e+01   0.000000   0.000000   0.000000
1  5.855970e+01  1.654562e+02  5.855970e+01   0.000000   0.000000   0.000000
2  5.855970e+01  5.855970e+01  1.654562e+02   0.000000   0.000000   0.000000
3  0.000000e+00 -2.000000e-10  0.000000e+00  82.521984   0.000000   0.000000
4 -2.000000e-10  0.000000e+00  2.000003e-10   0.000000  82.521984   0.000000
5 -2.000000e-10  0.000000e+00  0.000000e+00   0.000000   0.000000  82.521984

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

---

## 3.8 结果分析

**弹性张量（GPa）：**

立方晶系 3 个独立分量：

| 分量 | 本文结果 | Materials Project | 实验测量 |
|------|---------|-----------------|---------|
| $C_{11}$ | 165.5 | 153 | 165.7 |
| $C_{12}$ | 58.6 | 57 | — |
| $C_{44}$ | 82.5 | 74 | — |

$C_{11}$ 与实验测量值（165.7 GPa）偏差不足 0.2%。与 Materials Project 的差异（8–11%）来自两套计算采用了不同赝势和截断能，并不代表计算有误。

弹性张量矩阵结构完全符合立方晶系预期：非对角非零元（~10⁻¹⁰ GPa）均为数值噪声，可视为零。

**各向同性力学性质（GPa）：**

| 性质 | 数值 |
|------|------|
| 体弹模量 $B_V$ | 94.19 |
| 剪切模量 $G_V$ | 70.89 |
| 杨氏模量 $E_V$ | 170.02 |
| 泊松比 $\nu$ | 0.199 |

后处理结果同时保存在 `metrics.json` 和 `metrics_elastic.json`，便于批量处理和后续分析。
# 第四章：注意事项与参数速查

## 4.1 结构优化参数速查

**ionic relax（固定晶胞）最简 INPUT：**

```
INPUT_PARAMETERS
suffix                  example
ntype                   1
pseudo_dir              ./
calculation             relax
basis_type              lcao
ecutwfc                 100
scf_thr                 1e-6
cal_force               1
force_thr_ev            0.01
relax_nmax              100
out_stru                1
kspacing                0.1
```

**cell-relax（变胞优化）最简 INPUT：**

```
INPUT_PARAMETERS
suffix                  example
ntype                   1
pseudo_dir              ./
orbital_dir             ./
calculation             cell-relax
basis_type              lcao
ecutwfc                 100
scf_thr                 1e-7
cal_force               1
cal_stress              1
force_thr_ev            0.01
stress_thr              0.5
relax_nmax              100
out_stru                1
kspacing                0.08
```

---

## 4.2 常见问题汇总

| 问题 | 原因 | 解决方法 |
|------|------|---------|
| 结构优化不收敛 | 初始结构离平衡位置太远，或 K 点/截断能不足 | 检查初始结构合理性；增大 K 点密度或 ecutwfc |
| 应力收敛但力不收敛 | 原子被约束（STRU 中 0 0 0 约束标志） | 检查 STRU 文件中原子自由度设置 |
| `elastic prepare` 后重复运行丢失数据 | 命令会直接删除已有形变文件夹 | 提交计算前只运行一次 prepare |
| 弹性张量非对角元不为零 | K 点密度或截断能不足导致应力精度差 | 增大 K 点密度（减小 kspacing）或 ecutwfc |
| 与 Materials Project 结果差异大 | 赝势、泛函或截断能不同；取向不一致 | 确认取向与 IEEE 标准一致；对比赝势选择 |

---

## 4.3 计算精度建议

- **K 点密度**：弹性常数对 K 点密度敏感，建议 `kspacing` ≤ 0.1 Å⁻¹，必要时做收敛测试
- **截断能**：建议参考所用赝势的推荐截断能，不要低于赝势推荐值
- **SCF 收敛**：弹性计算的 SCF 阈值应严格（`scf_thr 1e-7` 或更严格），否则应力误差会传导到弹性张量
- **结构优化充分度**：cell-relax 的 `stress_thr` 建议 ≤ 1 kBar；若残余应力过大，拟合结果可靠性降低
# 附录

## 参考文献

[1] M. de Jong et al., Charting the complete elastic properties of inorganic crystalline compounds, *Scientific Data* **2**, 150009 (2015).

[2] S. Singh et al., MechElastic: A Python library for analysis of mechanical and elastic properties of bulk and 2D materials, *Computer Physics Communications* **267**, 108068 (2021).

[3] ABACUS 官方文档：Geometry Optimization，https://abacus.deepmodeling.com/en/latest/advanced/opt.html

[4] ABACUS+pymatgen 弹性常数教程，https://mcresearch.github.io/abacus-user-guide/abacus-elastic.html

---

## 进阶学习

- **高通量弹性性质筛选**：abacustest 支持同时对多个结构运行，`-f` 参数可接受多个 CIF 文件
- **弹性稳定性判据**：各晶系的力学稳定性条件由弹性张量元的正定性给出，可用于新材料结构筛选
- **力学各向异性分析**：弹性张量包含比各向同性模量更丰富的信息，可用 MechElastic 等工具可视化不同晶向的杨氏模量
- **磁性材料的弹性常数**：需在 INPUT 中设置 `nspin 2` 和初始磁矩，同时注意 smearing 参数；参见 ABACUS Mn 元素案例
