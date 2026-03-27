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
