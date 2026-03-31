# 第三章：Al 弹性常数计算实战

本章以铝（Al）的 FCC 结构为例，详细演示弹性常数计算的完整流程。Al 是典型的金属材料，属于立方晶系，只有 3 个独立弹性常数。

## 3.1 案例介绍

### 材料体系

**材料：** 铝（Aluminum, Al）
**晶体结构：** FCC（面心立方）
**空间群：** Fm-3m (225)
**晶格常数：** a ≈ 4.05 Å（实验值）
**原子数：** 4 个 Al 原子（FCC 原胞）

### 计算参数

**基组类型：** 平面波（pw）
**截断能：** 70 Ry
**k 点网格：** 6×6×6
**赝势：** Al.upf（ONCV 模守恒赝势）
**展宽方法：** Gaussian，σ = 0.015 Ry（金属体系）

### 两步计算流程

**步骤 1：晶胞优化（cell-relax）**
- 同时优化晶格参数和原子位置
- 得到平衡结构

**步骤 2：应变构型计算（relax）**
- 基于平衡结构生成 24 个应变构型
- 对每个构型计算应力
- 拟合弹性常数

## 3.2 步骤 1：晶胞优化（cell-relax）

### 3.2.1 INPUT 文件

```
INPUT_PARAMETERS

# System variables
calculation         cell-relax
symmetry            1
kspacing            0.14
precision           double

# Plane wave related variables
ecutwfc             70.0
pw_diag_nmax        20
pw_diag_ndim        2

# Electronic structure
basis_type          pw
ks_solver           dav_subspace
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-08

# Geometry relaxation
relax_method        cg
relax_nmax          60
cal_force           1
force_thr_ev        0.01
cal_stress          1
stress_thr          0.5
fixed_axes          None
```

**关键参数说明：**

- `calculation = cell-relax`：同时优化晶格和原子位置
- `kspacing = 0.14`：自动生成 k 点网格（单位：1/Bohr），对于 Al 约等于 6×6×6
- `ecutwfc = 70.0`：平面波截断能 70 Ry，对于 Al 已足够收敛
- `scf_thr = 1e-08`：高精度收敛阈值，确保应力计算准确
- `smearing_method = gauss`：Gaussian 展宽，适合金属
- `smearing_sigma = 0.015`：展宽宽度 0.015 Ry
- `mixing_type = broyden`：Broyden 混合方法，收敛性好
- `force_thr_ev = 0.01`：力收敛阈值 0.01 eV/Å
- `stress_thr = 0.5`：应力收敛阈值 0.5 kBar

### 3.2.2 STRU 文件

```
ATOMIC_SPECIES
Al 26.981500 Al.upf

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    4.05000000000     0.00000000000     0.00000000000
    0.00000000000     4.05000000000     0.00000000000
    0.00000000000     0.00000000000     4.05000000000

ATOMIC_POSITIONS
Cartesian

Al
0.000000
4
    0.00000000000     0.00000000000     0.00000000000 1 1 1
    0.00000000000     2.02500000000     2.02500000000 1 1 1
    2.02500000000     0.00000000000     2.02500000000 1 1 1
    2.02500000000     2.02500000000     0.00000000000 1 1 1
```

**文件说明：**

- `LATTICE_CONSTANT = 1.889726`：转换因子，使得 LATTICE_VECTORS 单位为 Å
- 晶格矢量：立方晶格，a = 4.05 Å（初始值，接近实验值）
- 坐标系：Cartesian（笛卡尔坐标）
- 4 个 Al 原子：FCC 结构的原胞
- `1 1 1`：允许原子在三个方向上移动

### 3.2.3 KPT 文件

```
K_POINTS
0
Gamma
6 6 6 0 0 0
```

**说明：**
- 使用 Gamma 中心的 Monkhorst-Pack 网格
- 6×6×6 k 点网格，对于 Al 金属已足够收敛

> **注意：** 也可以使用 INPUT 中的 `kspacing = 0.14` 自动生成 k 点，此时 KPT 文件内容会被忽略。

### 3.2.4 运行计算

```bash
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

计算完成后，优化后的结构保存在 `OUT.ABACUS/STRU_ION_D` 文件中。

### 3.2.5 检查收敛

查看输出日志，确认优化收敛：

```bash
grep "TOTAL-PRESSURE" OUT.ABACUS/running_cell-relax.log | tail -5
grep "LARGEST GRAD" OUT.ABACUS/running_cell-relax.log | tail -5
```

收敛判据：
- `LARGEST GRAD < 0.01 eV/Å`（力收敛）
- `TOTAL-PRESSURE < 0.5 kBar`（应力收敛）

优化后的晶格常数约为 **4.04 Å**，与实验值 4.05 Å 接近（PBE 泛函通常略微低估）。

## 3.3 步骤 2：应变构型计算（relax）

### 3.3.1 准备优化后的结构

将优化后的结构复制到新的工作目录：

```bash
mkdir elastic_calc
cd elastic_calc
cp ../OUT.ABACUS/STRU_ION_D ./STRU
```

### 3.3.2 修改 INPUT 文件

将 `calculation` 改为 `relax`（固定晶格，允许原子弛豫）：

**关键变化：**
- `calculation = relax`（不再是 cell-relax）
- 移除 `kspacing`，使用固定 k 点网格（KPT 文件）

其他参数（ecutwfc、scf_thr、smearing_method 等）与步骤 1（3.2.1）保持一致

### 3.3.3 STRU 文件（优化后）

优化后的 STRU 文件（`STRU_ION_D`）内容示例：

```
ATOMIC_SPECIES
Al 26.981500 Al.upf

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    4.04160318060     0.00000000000    -0.00000000000
    0.00000000000     4.04160318060    -0.00000000000
    0.00000000000    -0.00000000000     4.04160318060

ATOMIC_POSITIONS
Direct

Al
0.000000
4
    0.00000000000     0.00000000000     0.00000000000 1 1 1
    0.00000000000     0.50000000000     0.50000000000 1 1 1
    0.50000000000     0.00000000000     0.50000000000 1 1 1
    0.50000000000     0.50000000000     0.00000000000 1 1 1
```

**变化说明：**
- 晶格常数：从 4.05 Å 优化为 4.04160 Å
- 坐标系：从 Cartesian 变为 Direct（分数坐标）
- 原子位置：优化后的分数坐标

### 3.3.4 KPT 文件

```
K_POINTS
0
Gamma
6 6 6 0 0 0
```

保持与步骤 1 相同的 k 点设置。

### 3.3.5 生成应变构型

使用 pymatgen 脚本生成 24 个应变构型。这部分将在第四章详细介绍，这里先说明原理：

**应变构型生成原理：**
1. 读取优化后的 STRU 文件
2. 对晶格矢量施加 6 种应变模式
3. 每种应变施加 4 种大小：±0.5%, ±1.0%
4. 生成 24 个 task 文件夹（task.000 到 task.023）

**每个 task 文件夹包含：**
- INPUT：与上述相同
- KPT：与上述相同
- STRU：应变后的结构
- strain.json：记录应变大小

### 3.3.6 批量计算应力

对 24 个应变构型分别运行 ABACUS 计算：

```bash
for i in {000..023}; do
    cd task.$i
    OMP_NUM_THREADS=1 mpirun -np 8 abacus
    cd ..
done
```

每个计算完成后，应力数据保存在 `OUT.ABACUS/running_relax.log` 中。

### 3.3.7 提取应力数据

从每个 task 的输出中提取应力张量：

```bash
grep "TOTAL-STRESS" task.*/OUT.ABACUS/running_relax.log
```

应力输出格式示例：

```
TOTAL-STRESS (KBAR):
-0.12   0.00   0.00
 0.00  -0.12   0.00
 0.00   0.00  -0.12
```

单位为 kBar，需要转换为 GPa（1 kBar = 0.1 GPa）。

## 3.4 结果解读

### 3.4.1 弹性常数矩阵

通过 pymatgen 拟合得到的弹性常数矩阵（单位：GPa）：

```
C_11  C_12  C_12   0     0     0
C_12  C_11  C_12   0     0     0
C_12  C_12  C_11   0     0     0
 0     0     0    C_44   0     0
 0     0     0     0    C_44   0
 0     0     0     0     0    C_44
```

典型计算结果（Al，PBE 泛函）：
- C_11 ≈ 114 GPa
- C_12 ≈ 62 GPa
- C_44 ≈ 32 GPa

### 3.4.2 弹性模量

从弹性常数计算得到的弹性模量：

| 弹性模量 | 计算公式 | 典型值（Al） |
|---------|---------|-------------|
| 体弹模量 B | (C_11 + 2C_12) / 3 | ~79 GPa |
| 剪切模量 G | (C_11 - C_12 + 3C_44) / 5 | ~28 GPa |
| 杨氏模量 E | 9BG / (3B + G) | ~72 GPa |
| 泊松比 ν | (3B - 2G) / (6B + 2G) | ~0.34 |

### 3.4.3 与实验值对比

| 物理量 | 计算值（PBE） | 实验值 | 偏差 |
|--------|--------------|--------|------|
| C_11 | ~114 GPa | 108 GPa | +5.6% |
| C_12 | ~62 GPa | 62 GPa | 0% |
| C_44 | ~32 GPa | 28 GPa | +14% |
| 体弹模量 B | ~79 GPa | 76 GPa | +3.9% |
| 杨氏模量 E | ~72 GPa | 70 GPa | +2.9% |

PBE 泛函对 Al 的弹性常数预测较为准确，偏差在 15% 以内。

