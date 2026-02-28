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

```
Si-elastic/
├── deform-norm_xx-0.005/    # 正应变 xx，幅度 +0.5%
├── deform-norm_xx+0.005/    # 正应变 xx，幅度 -0.5%
├── deform-norm_yy-0.005/
├── ...                       # 3 种正应变 × 4 幅度 = 12 个文件夹
├── deform-shear_yz-0.005/   # 剪切应变 yz，幅度 ±0.5%
├── ...                       # 3 种剪切应变 × 4 幅度 = 12 个文件夹
└── org/                      # 原始未形变结构
```

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
Si-elastic/     90.705919      65.134208     157.664088       0.210302

Si-elastic/     elastic_tensor:
            0             1             2          3          4          5
0  155.464972  5.832639e+01  5.832639e+01   0.000000   0.000000   0.000000
1   58.326393  1.554650e+02  5.832639e+01   0.000000   0.000000   0.000000
2   58.326393  5.832639e+01  1.554650e+02   0.000000   0.000000   0.000000
3    0.000000  0.000000e+00  0.000000e+00  76.177486   0.000000   0.000000
4    0.000000 -2.000000e-10  0.000000e+00   0.000000  76.177486   0.000000
5    0.000000  0.000000e+00 -2.000000e-10   0.000000   0.000000  76.177486

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

## 3.7 结果分析

**弹性张量（单位：GPa）：**

立方晶系的 3 个独立分量从输出中读取：

| 分量 | 本文结果 | Materials Project |
|------|---------|-----------------|
| C₁₁ | 155.5 | 153 |
| C₁₂ | 58.2 | 57 |
| C₄₄ | 76.2 | 74 |

计算结果与 Materials Project 数据库吻合，差异在 2% 以内。

弹性张量矩阵结构也与立方晶系的预期完全符合：C₁₁ = C₂₂ = C₃₃ ≈ 155.5 GPa，C₁₂ = C₁₃ = C₂₃ ≈ 58.3 GPa，C₄₄ = C₅₅ = C₆₆ ≈ 76.2 GPa，非对角非零元素均为数值噪声量级（~10⁻¹⁰ 量级），可认为为 0。

**各向同性力学性质（单位：GPa）：**

| 性质 | 数值 |
|------|------|
| 体模量（Bulk Modulus） | 90.71 |
| 剪切模量（Shear Modulus） | 65.13 |
| 杨氏模量（Young's Modulus） | 157.66 |
| 泊松比（Poisson's Ratio） | 0.210 |

这些各向同性弹性性质由 abacustest 通过 Voigt-Reuss-Hill 平均方法从弹性张量导出，代表多晶材料的平均力学行为。

后处理结果同时保存在 `metrics.json` 和 `metrics_elastic.json` 文件中，便于批量处理。
