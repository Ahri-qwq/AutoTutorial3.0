# 第五章：结构优化实战

前几章分别介绍了三个输入文件。本章以一个具体例子，展示三个文件如何配合完成一次结构优化计算。

## 5.1 结构优化的物理意义

实验制备或数据库中的晶体结构并非总处于 DFT 意义下的能量极小值。结构优化（几何弛豫）通过调整原子位置和晶格参数，使体系达到理论上的稳定构型，是大多数材料计算工作的第一步。

**relax 与 cell-relax 的区别：**

| 类型 | 优化对象 | 需要额外设置 |
|------|----------|--------------|
| `relax` | 原子位置（晶格固定）| `force_thr_ev` |
| `cell-relax` | 原子位置 + 晶格矢量 | `force_thr_ev` + `stress_thr` |

通常建议先做 `cell-relax` 同时弛豫晶格和原子，再视情况做 `relax` 精细优化原子位置。

## 5.2 完整案例：Si 晶体 cell-relax

以金刚石结构硅（8 原子超胞）为例，演示一次完整的结构优化。

### INPUT 文件

```
INPUT_PARAMETERS

suffix          Si
ntype           1
pseudo_dir      ./
calculation     cell-relax
basis_type      pw

ecutwfc         50
scf_thr         1e-6
scf_nmax        100

smearing_method gauss
smearing_sigma  0.01

symmetry        1

force_thr_ev    0.01
stress_thr      0.5
relax_nmax      100
out_stru        1
```

### KPT 文件

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

### STRU 文件

```
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
5.50000   0.00000   0.00000
0.00000   5.50000   0.00000
0.00000   0.00000   5.50000

ATOMIC_POSITIONS
Direct
Si
0.0
8
0.000  0.000  0.000  1  1  1
0.250  0.250  0.250  1  1  1
0.500  0.500  0.000  1  1  1
0.750  0.750  0.250  1  1  1
0.500  0.000  0.500  1  1  1
0.750  0.250  0.750  1  1  1
0.000  0.500  0.500  1  1  1
0.250  0.750  0.750  1  1  1
```

注意：初始晶格常数故意设为 5.50 Å（实验值约 5.43 Å），目的是观察 cell-relax 对晶格的优化效果。

### 运行计算

```bash
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

## 5.3 判断优化是否收敛

每完成一步离子迭代，输出中会显示：

```
TOTAL-PRESSURE: -2.070e+00 KBAR
ETOT DIFF (eV)       : -1.234e-03
LARGEST GRAD (eV/A)  : 8.521e-03
```

**收敛判据（cell-relax 需同时满足）：**

| 量 | 含义 | 收敛条件 |
|----|------|----------|
| `LARGEST GRAD (eV/A)` | 所有原子中最大受力 | < `force_thr_ev`（本例 0.01 eV/Å）|
| `TOTAL-PRESSURE (KBAR)` | 晶胞总压强（应力迹的三分之一）| < `stress_thr`（本例 0.5 kBar）|

两个条件同时满足时，计算收敛并自动停止。`relax` 只需满足力的判据。

## 5.4 查看优化结果

**优化轨迹：** 设置 `out_stru 1` 后，每步优化的结构保存在 `OUT.Si/STRU_ION_D`，可用于追踪优化过程。

**最终结构：** 优化收敛后，最终的晶格参数和原子位置会出现在 `running_cell-relax.log`（或屏幕输出）中：

```
 Cell Parameters:  (Angstrom)
  a = 5.43020   b = 5.43020   c = 5.43020
  alpha = 90.000   beta = 90.000   gamma = 90.000
```

优化后的晶格常数（约 5.43 Å）与实验值吻合。最终结构也以 STRU 格式保存在 `OUT.Si/` 下，可直接用于后续计算。
