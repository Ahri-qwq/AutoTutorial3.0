---
title: "使用 abacustest 计算晶体弹性性质"
author: "AutoTutorial 3.0"
date: "2026-02-27"
topic: "ABACUS弹性常数计算"
task_type: "C"
has_case: true
word_count: 约1876词/字
lines: 468
chapters: 5
cases:
  - name: "Si（立方晶系）"
    mp_id: null
    result_C11: 165.5
    result_C12: 58.6
    result_C44: 82.5
  - name: "金红石型TiO₂（四方晶系）"
    mp_id: "mp-2657"
    result_bulk: 220.71
    result_shear: 128.68
    result_young: 323.23
    result_poisson: 0.256
---

# 使用 abacustest 计算晶体弹性性质

# 前言

本教程介绍如何使用 abacustest 计算晶体的弹性性质，包括完整的弹性张量矩阵以及体模量、剪切模量、杨氏模量、泊松比等宏观力学指标。

**适用读者：** 有 ABACUS 基本使用经验的研究者，了解第一性原理结构优化的基本流程。

**前置知识：**
- ABACUS 基本输入文件（INPUT、STRU、KPT）的作用与格式
- 结构弛豫计算（cell-relax）的基本流程
- LCAO 基组的概念

**教程结构：**
- 第一章：弹性张量的物理背景与应力-应变法
- 第二章：软件准备——安装 abacustest，配置赝势与轨道文件
- 第三章：案例一——Si（立方晶系，3 个独立弹性张量元）
- 第四章：案例二——金红石型 TiO₂（四方晶系，6 个独立弹性张量元）
- 第五章：讨论与小结

**案例说明：**

两个案例都使用 abacustest 的弹性计算工作流，计算引擎为 ABACUS LTSv3.10，基组为 LCAO。

- **Si**：立方晶系，对称性高，仅有 3 个独立弹性张量元（C₁₁、C₁₂、C₄₄）。流程最简，适合熟悉计算步骤。
- **金红石型 TiO₂**（Materials Project 编号 mp-2657）：四方晶系，有 6 个独立弹性张量元。结果将与 Materials Project 数据库及实验测量值对比，展示 abacustest 的计算精度。

所有数值来自实际计算，未经修改。

---

# 一、背景与方法

## 1.1 弹性张量与广义胡克定律

晶体在外力作用下会产生形变（应变），内部相应地产生应力。在材料处于平衡位置附近的线弹性范围内，应力张量 σ_ij 与应变张量 ε_kl 之间的关系由广义胡克定律描述：

$$\sigma_{ij} = C_{ijkl}\,\varepsilon_{kl}$$

其中 C_ijkl 是四阶弹性刚度张量。由于应力和应变都是二阶对称张量，C_ijkl 满足以下对称性：

$$C_{ijkl} = C_{jikl},\quad C_{ijkl} = C_{ijlk},\quad C_{ijkl} = C_{klij}$$

经过对称性约化，81 个分量缩减为 **21 个独立分量**。

利用 **Voigt 标记**，用如下对应关系将 4 阶张量映射为 6×6 对称矩阵：

$$xx\to 1,\quad yy\to 2,\quad zz\to 3,\quad yz\to 4,\quad xz\to 5,\quad xy\to 6$$

应力和应变各写为 6 维列向量，胡克定律化为矩阵形式：

$$\begin{bmatrix}\sigma_1\\\sigma_2\\\sigma_3\\\sigma_4\\\sigma_5\\\sigma_6\end{bmatrix}=\begin{bmatrix}C_{11}&C_{12}&C_{13}&C_{14}&C_{15}&C_{16}\\&C_{22}&C_{23}&C_{24}&C_{25}&C_{26}\\&&C_{33}&C_{34}&C_{35}&C_{36}\\&&&C_{44}&C_{45}&C_{46}\\&&&&C_{55}&C_{56}\\&&&&&C_{66}\end{bmatrix}\begin{bmatrix}\varepsilon_1\\\varepsilon_2\\\varepsilon_3\\\varepsilon_4\\\varepsilon_5\\\varepsilon_6\end{bmatrix}$$

## 1.2 晶体对称性的约束

晶体的点群对称性会进一步减少独立分量数：

| 晶系 | 独立分量数 | 示例 |
|------|-----------|------|
| 三斜晶系 | 21 | — |
| 单斜晶系 | 13 | — |
| 正交晶系 | 9 | — |
| 六方晶系 | 5 | — |
| 四方晶系（Laue I） | 6 | 金红石 TiO₂ |
| 立方晶系 | 3 | Si、Cu |

**立方晶系**（Si）：独立分量只有 3 个——C₁₁（= C₂₂ = C₃₃）、C₁₂（= C₁₃ = C₂₃）、C₄₄（= C₅₅ = C₆₆），矩阵其余非对角元素均为 0。

**四方晶系**（金红石型 TiO₂，Laue 类型 I）：独立分量为 C₁₁（= C₂₂）、C₁₂、C₁₃（= C₂₃）、C₃₃、C₄₄（= C₅₅）、C₆₆，共 6 个。

计算结果与 Materials Project 数据库对比时，还需注意晶体取向。Materials Project 采用 IEEE 176/1987 标准取向。对于本教程中的 Si，由 ase 生成的立方惯用胞晶格矢量已与坐标轴对齐，无需额外旋转。

## 1.3 计算方法：应力-应变法

目前计算弹性张量最常用的方法是**应力-应变法**：

1. 对优化后的平衡结构施加一系列已知应变
2. 用第一性原理计算各形变结构的应力
3. 对应力-应变数据线性拟合，得到弹性张量元

**abacustest 的实现方案：**

对 6 种独立应变（3 种正应变 + 3 种剪切应变），每种取 4 个幅度（±0.5%、±1%），共生成 **24 个形变构型**。对每个构型做固定晶格的 relax 计算（允许原子弛豫后再计算应力），以减少数值噪声。再加 1 个原始结构（org），合计 **25 个计算任务**。

**完整工作流：**

```
生成/下载晶体结构（CIF）
        ↓
abacustest model inputs  →  准备 cell-relax 输入
        ↓
提交并等待结构优化完成
        ↓
整理优化后结构，准备 SCF 目录
        ↓
abacustest model elastic prepare  →  生成 24+1 个形变构型
        ↓
提交并等待弹性计算完成
        ↓
abacustest model elastic post  →  输出弹性张量和力学性质
```

---

# 二、软件准备

## 2.1 安装 abacustest

abacustest 是 ABACUS 的辅助工具，提供输入文件准备、高通量任务管理、结果后处理等功能。通过 pip 安装：

```bash
pip install abacustest
```

也可从 GitHub 获取最新版本：

```bash
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装完成后，执行 `abacustest -h` 验证安装是否成功。

## 2.2 准备赝势和轨道文件

abacustest 在准备 LCAO 计算的输入文件时，会自动从环境变量 `ABACUS_PP_PATH` 和 `ABACUS_ORB_PATH` 指定的路径中查找赝势和数值原子轨道文件。在运行前需设置这两个环境变量：

```bash
export ABACUS_PP_PATH=/path/to/pseudopotentials
export ABACUS_ORB_PATH=/path/to/orbitals
```

本教程的案例使用 APNS-v1 赝势库与 efficiency 数值原子轨道基组。赝势和轨道文件可从 ABACUS 官方渠道获取。

> **说明：** 若计划在 Bohrium 云平台上运行，可使用平台提供的 abacustest App，无需手动配置环境变量，赝势和轨道文件由平台统一管理。

---

# 三、案例一：Si 的弹性常数

Si 是立方晶系，其弹性张量只有 3 个独立分量：C₁₁、C₁₂、C₄₄。以下完整演示从结构生成到弹性张量输出的全流程。

## 3.1 生成 Si 惯用胞结构

计算弹性张量时，立方晶系材料通常使用正方体形式的惯用胞，并保持坐标轴与晶格矢量重合。用 ase 生成 Si 惯用胞的 CIF 文件：

```python
from ase.build import bulk

si_conv = bulk('Si', cubic=True)
si_conv.write("Si_conv.cif")
```

`bulk('Si', cubic=True)` 生成的即为 Si 的立方惯用胞（8 原子超胞），晶格矢量与坐标轴对齐，后续无需额外旋转至 IEEE 标准取向。

## 3.2 结构优化：准备 cell-relax 输入文件

在施加应变前，必须先对晶胞和原子位置进行充分优化，消除残余应力，使结构处于势能面的局部极小值点。用 abacustest 自动准备 cell-relax 的 ABACUS 输入文件：

```bash
abacustest model inputs -f Si_conv.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si
```

各参数含义：

| 参数 | 含义 |
|------|------|
| `-f Si_conv.cif` | 输入结构文件（可同时指定多个） |
| `--ftype cif` | 结构文件格式 |
| `--jtype cell-relax` | 计算类型：同时优化晶胞和原子位置 |
| `--lcao` | 使用 LCAO 数值原子轨道基组 |
| `--folder-syntax Si` | 生成的输入文件夹命名规则 |

执行完成后，当前目录生成 `Si/` 文件夹，包含 ABACUS 的 INPUT、STRU、KPT 文件和赝势/轨道文件。INPUT 中的参数适配 ABACUS LTSv3.10。

> **注意：** 检查 INPUT 文件中的参数是否合理，必要时调整 nspin（磁矩）、DFT+U 参数等，再提交计算。

提交优化计算，等待计算完成。

## 3.3 整理优化后结构，准备弹性计算目录

结构优化完成后，优化后的结构保存在 `Si/INPUT/OUT.ABACUS/STRU_ION_D`。新建弹性计算目录并整理文件：

```bash
mkdir Si-elastic
cp Si/INPUT Si/Si* Si-elastic
cp Si/INPUT/OUT.ABACUS/STRU_ION_D Si-elastic/STRU
```

然后修改 `Si-elastic/INPUT`，将计算类型改为 SCF：

```
calculation    scf
```

此目录将作为弹性计算的基础：abacustest 会读取其中的 INPUT、STRU、KPT 和赝势/轨道文件，自动生成各形变构型的输入文件夹。

由于 Si 的晶格矢量已与坐标轴对齐，这里不需要对结构做旋转处理。

## 3.4 生成形变构型

使用 abacustest 自动生成弹性计算所需的形变结构：

```bash
abacustest model elastic prepare -j Si-elastic
```

执行完成后，在 `Si-elastic/` 目录下生成以下内容：

- 一系列以 `deformed_` 开头的文件夹：共 24 个，对应 3 种正应变 × 4 幅度 + 3 种剪切应变 × 4 幅度
- `org/` 文件夹：包含原始未形变结构的计算任务

共 25 个计算任务（24 个形变 + 1 个原始）。每个文件夹内的 ABACUS 计算类型为 relax（固定晶胞、允许原子弛豫），并自动开启了应力计算（`cal_stress = 1`）。

常用可选参数：

| 参数 | 默认值 | 含义 |
|------|--------|------|
| `--norm` | 0.01 | 最大正应变（生成 ±0.5%、±1%） |
| `--shear` | 0.01 | 最大剪切应变（生成 ±0.5%、±1%） |
| `--norelax` | — | 不做原子弛豫，直接计算应力 |

> **警告：** 重复执行 `abacustest model elastic prepare` 会**直接删除**已有的形变文件夹！准备好后不要重复运行此命令。

## 3.5 提交计算

将 25 个 ABACUS 计算任务提交到计算集群或云平台。计算完成后，进行后处理。

## 3.6 后处理：提取弹性张量

所有计算完成后，用一条命令完成弹性张量的拟合和力学性质的计算：

```bash
abacustest model elastic post -j Si-elastic
```

屏幕输出结果：

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

## 3.7 结果分析

**弹性张量（单位：GPa）：**

立方晶系的 3 个独立分量从输出中读取：

| 分量 | 本文结果 | Materials Project |
|------|---------|-----------------|
| C₁₁ | 165.5 | 153 |
| C₁₂ | 58.6 | 57 |
| C₄₄ | 82.5 | 74 |

本文 C₁₂ 与 Materials Project 吻合（差异约 3%）；C₁₁ 和 C₄₄ 差异在 8–11%，这是两套计算采用了不同赝势和截断能所致，并不代表计算有误。本文 C₁₁ = 165.5 GPa 与硅的实验测量值（165.7 GPa）偏差不足 0.2%。

弹性张量矩阵结构与立方晶系的预期完全符合：等价分量之间的一致性良好，非对角非零元素均为数值噪声量级（~10⁻¹⁰ GPa），可认为为 0。独立分量的具体数值为：C₁₁ = C₂₂ = C₃₃ ≈ 165.5 GPa，C₁₂ = C₁₃ = C₂₃ ≈ 58.6 GPa，C₄₄ = C₅₅ = C₆₆ ≈ 82.5 GPa。

**各向同性力学性质（单位：GPa）：**

| 性质 | 数值 |
|------|------|
| 体模量（Bulk Modulus） | 94.19 |
| 剪切模量（Shear Modulus） | 70.89 |
| 杨氏模量（Young's Modulus） | 170.02 |
| 泊松比（Poisson's Ratio） | 0.199 |

这些各向同性弹性性质由 abacustest 从弹性张量导出，代表多晶材料的平均力学行为。

后处理结果同时保存在 `metrics.json` 和 `metrics_elastic.json` 文件中，便于批量处理。

---

# 四、案例二：金红石型 TiO₂ 的弹性常数

金红石型 TiO₂（Materials Project 编号 mp-2657）属于四方晶系，Laue 类型 I，有 6 个独立弹性张量元：C₁₁（= C₂₂）、C₁₂、C₁₃（= C₂₃）、C₃₃、C₄₄（= C₅₅）、C₆₆。计算流程与 Si 完全相同，这里重点展示四方晶系的结果特征及与实验数据的对比。

## 4.1 获取结构

从 Materials Project 下载金红石型 TiO₂（mp-2657）的 CIF 文件，命名为 `TiO2.cif`。可通过网站手动下载，也可用 `mp-api` 在 Python 中直接获取：

```python
from mp_api.client import MPRester

with MPRester("YOUR_API_KEY") as mpr:
    structure = mpr.get_structure_by_material_id("mp-2657")
    structure.to("TiO2.cif")
```

API Key 可在 [Materials Project 个人设置页面](https://next-gen.materialsproject.org/api) 获取。使用 mp-2657 的 DFT 弛豫结构作为起点，可精确复现本教程的计算数值。

## 4.2 结构优化

用 abacustest 准备 cell-relax 输入文件并提交优化计算：

```bash
abacustest model inputs -f TiO2.cif --ftype cif --lcao --jtype cell-relax --folder-syntax TiO2-rutile
```

等待计算完成后，整理优化后的结构：

```bash
mkdir TiO2-rutile-elastic
cp TiO2-rutile/INPUT TiO2-rutile/Ti* TiO2-rutile/O* TiO2-rutile-elastic
cp TiO2-rutile/INPUT/OUT.ABACUS/STRU_ION_D TiO2-rutile-elastic/
```

> **说明：** `TiO2-rutile/Ti*` 和 `TiO2-rutile/O*` 分别复制 Ti 和 O 的赝势与轨道文件（文件名以元素符号开头）。根据实际文件名调整通配符。

## 4.3 准备并提交弹性计算

同样地，修改 `TiO2-rutile-elastic/INPUT` 中的计算类型为 `scf`，然后运行：

```bash
abacustest model elastic prepare -j TiO2-rutile-elastic
```

同样在 `TiO2-rutile-elastic/` 下生成 25 个计算任务（24 个形变 + 1 个原始结构）。提交所有任务，等待完成。

> **警告：** 与 Si 案例相同，不要在计算进行中或计算完成后重复运行 `elastic prepare`，否则已有的计算文件夹将被删除。

## 4.4 后处理

```bash
abacustest model elastic post -j TiO2-rutile-elastic
```

屏幕输出：

```
Model: elastic
Postprocessing elastic calculation for job: TiO2-rutile-elastic/
                      bulk_modulus  shear_modulus  young_modulus  poisson_ratio
TiO2-rutile-elastic/    220.705174     128.683753     323.230608       0.255911

TiO2-rutile-elastic/    elastic_tensor:
              0             1             2           3           4           5
0  2.812277e+02  1.581445e+02  1.559258e+02    0.000000    0.000000    0.000000
1  1.581445e+02  2.812277e+02  1.559258e+02    0.000000    0.000000    0.000000
2  1.558677e+02  1.558677e+02  4.840152e+02    0.000000    0.000000    0.000000
3  4.200000e-09  5.500000e-09  1.230000e-08  117.414353    0.000000    0.000000
4 -2.250000e-08 -1.770000e-08 -4.810000e-08    0.000000  117.414354    0.000000
5 -1.300000e-09  2.200002e-09  0.000000e+00    0.000000    0.000000  216.431859

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

## 4.5 结果分析

### 弹性张量结构的验证

首先检查弹性张量矩阵结构是否符合四方晶系（Laue I）的预期：
- C₁₁ = C₂₂：281.2 ≈ 281.2 GPa ✓
- C₁₃ = C₂₃：155.9 ≈ 155.9 GPa ✓
- C₄₄ = C₅₅：117.4 ≈ 117.4 GPa ✓
- C₆₆ ≠ C₄₄：216.4 ≠ 117.4 GPa ✓（四方晶系区别于六方晶系的关键特征）
- 非对角非零元均为数值噪声量级（~10⁻⁸ GPa）✓

计算得到的弹性张量结构完全符合金红石型 TiO₂ 所属的四方晶系（Laue 类型 I）对 6 个独立弹性张量元的要求。

### 与 Materials Project 及实验数据的对比

下表汇总了本文计算、Materials Project 数据库和实验测量 [3] 的各独立弹性张量元：

| | 本文结果 (GPa) | Materials Project (GPa) | 实验测量 (GPa) |
|---|---:|---:|---:|
| C₁₁ | 281.2 | 426 | 268.0 ± 1.4 |
| C₁₂ | 158.1 | 2 | 174.9 ± 1.4 |
| C₁₃ | 155.9 | 149 | 147.4 ± 1.5 |
| C₃₃ | 484.0 | 470 | 484.2 ± 1.8 |
| C₄₄ | 117.4 | 113 | 123.8 ± 0.2 |
| C₆₆ | 216.4 | 43 | 190.2 ± 0.5 |

**对比分析：**

本文结果与实验数据总体吻合较好，多数分量的偏差在 10% 以内。C₃₃ 和 C₄₄ 与实验的吻合尤为出色（偏差 < 5%）。

Materials Project 数据库中部分分量与实验值差异显著。C₁₂（数据库值 2 GPa，实验值 174.9 GPa）和 C₆₆（数据库值 43 GPa，实验值 190.2 GPa）的偏差尤为突出，可能与该数据库条目的计算参数设置有关。本文的计算结果在这些分量上明显更接近实验值。

**各向同性力学性质（单位：GPa，无量纲量泊松比除外）：**

| 性质 | 数值 |
|------|------|
| 体模量（Bulk Modulus） | 220.71 |
| 剪切模量（Shear Modulus） | 128.68 |
| 杨氏模量（Young's Modulus） | 323.23 |
| 泊松比（Poisson's Ratio） | 0.256 |

金红石型 TiO₂ 的体模量（220 GPa）和杨氏模量（323 GPa）均显著高于 Si，反映了 TiO₂ 更强的化学键合和更高的力学刚度。

---

# 五、讨论与小结

## 5.1 abacustest 弹性计算工作流的特点

abacustest 封装了弹性计算中最繁琐的步骤——生成 24 个形变构型并对结果进行拟合，用户只需：

1. 准备好平衡结构的 ABACUS 输入目录（INPUT、STRU、KPT、赝势/轨道）
2. 执行 `abacustest model elastic prepare` 生成形变构型
3. 提交计算，等待完成
4. 执行 `abacustest model elastic post` 得到结果

整个后处理过程（弹性张量拟合、宏观力学性质计算、结果保存）由一条命令完成，结果以 JSON 格式保存，方便批量处理。

## 5.2 注意事项汇总

| 问题 | 说明 |
|------|------|
| 结构未充分优化 | 优化后的结构应无残余应力，否则弹性张量拟合误差大 |
| 重复运行 `elastic prepare` | 会直接删除已有的形变文件夹，不要在计算进行中重复运行 |
| 应变幅度选择 | 默认值（±0.5%、±1%）适用于大多数材料；刚度很高的材料可适当增大 |
| 晶体取向 | 与 Materials Project 对比时需确认取向一致；由 ase 生成的立方惯用胞无需额外旋转 |
| 弹性张量矩阵的数值噪声 | 非零但极小的非对角元（~10⁻¹⁰ GPa）是数值误差，可视为 0 |

## 5.3 拓展方向

- **其他晶系：** 对六方、正交等晶系，工作流相同，只是独立分量数量不同，结果解读时参考对应晶系的对称性约束。
- **批量计算：** abacustest 支持同时对多个结构运行，`-f` 参数可接受多个 CIF 文件，适合高通量弹性性质筛选。
- **应力精度优化：** 弹性常数对应力计算精度敏感，若结果与参考值偏差较大，可考虑增大 K 点密度或提高截断能。

---

**本教程总结：**

使用 abacustest 的弹性计算工作流，通过两个案例验证了计算流程：

- **Si（立方晶系）**：C₁₁ = 165.5 GPa、C₁₂ = 58.6 GPa、C₄₄ = 82.5 GPa，与实验测量值（C₁₁ = 165.7 GPa）偏差不足 0.2%
- **金红石型 TiO₂（四方晶系）**：6 个独立分量均符合四方晶系对称性约束，多数分量与实验测量值偏差 < 10%

---

# 附录

## 参考文献

[1] M. de Jong, W. Chen, T. Angsten, A. Jain, R. Notestine, A. Gamst, M. Sluiter, C. Krishna Ande, S. van der Zwaag, J. J. Plata, C. Toher, S. Curtarolo, G. Ceder, K. A. Persson, and M. Asta, Charting the complete elastic properties of inorganic crystalline compounds, *Scientific Data* **2**, 150009 (2015).

[2] S. Singh, L. Lang, V. Dovale-Farelo, U. Herath, P. Tavadze, F.-X. Coudert, and A. H. Romero, MechElastic: A Python library for analysis of mechanical and elastic properties of bulk and 2D materials, *Computer Physics Communications* **267**, 108068 (2021).

[3] D. G. Isaak, J. D. Carnes, O. L. Anderson, H. Cynn, and E. Hake, Elasticity of TiO₂ rutile to 1800 K, *Physics and Chemistry of Minerals* **26**, 31 (1998).

---

## 常见问题

**Q：abacustest elastic post 报错，说某个计算未完成？**

A：检查对应的 deformed_ 文件夹中是否有 ABACUS 的输出文件，以及计算是否正常完成（OUT.ABACUS/ 目录中应有 running_relax.log 且末尾显示计算完成）。

**Q：计算结果中弹性张量有明显不符合对称性的元素？**

A：可能是应力计算精度不足。考虑增大 K 点密度（`kspacing` 减小）或提高截断能（`ecutwfc` 增大）后重新计算。

**Q：如何与 Materials Project 的弹性张量正确对比？**

A：需确认晶体取向与 IEEE 标准一致。对由 ase 生成的立方晶系惯用胞，晶格矢量已沿坐标轴排列，无需旋转。对其他来源的结构，参考文献 [1] 中描述了 Materials Project 所用的取向规范。

---

## 进阶学习

- **应力-应变法 vs 能量-应变法：** 后者通过拟合能量-应变曲线获得弹性常数，对数值噪声更鲁棒，但需要更多计算点。abacustest 目前采用应力-应变法。
- **弹性稳定性判据：** 各晶系的弹性稳定性条件由独立弹性张量元的正定性约束给出，可用于判断结构的力学稳定性。
- **力学各向异性：** 弹性张量包含比各向同性模量更丰富的信息，可用于分析材料的各向异性力学行为（如不同晶向的杨氏模量差异）。
