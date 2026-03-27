# 第二章：INPUT 文件详解

INPUT 文件是 ABACUS 计算的核心控制文件，决定了计算类型、使用的基组、精度要求以及输出内容。

## 2.1 文件格式与基本规则

INPUT 文件必须以 `INPUT_PARAMETERS` 作为第一行（有效内容），之前的内容会被忽略。参数格式为：

```
关键字   值
```

注意：**不使用等号**（区别于 VASP 的 INCAR）。

其他规则：
- 以 `#` 或 `/` 开头的行为注释，整行忽略
- 同一参数出现多次时，**取最后一次**的值（区别于 VASP 取第一次）
- 布尔值可用 `1`/`0`、`True`/`False`、`T`/`F`，大小写均可
- 文件名固定为 `INPUT`，不可更改

一个最简 INPUT 示例（Si 的 SCF 计算）：

```
INPUT_PARAMETERS
ntype           1
pseudo_dir      ./
ecutwfc         50
basis_type      pw
calculation     scf
```

## 2.2 通用参数

这些参数在几乎所有计算中都需要设置：

### suffix

```
suffix          ABACUS
```

输出文件夹的后缀名。计算完成后会生成 `OUT.ABACUS/` 文件夹。建议设为体系名称，如 `suffix Si`，则输出到 `OUT.Si/`。

### ntype

```
ntype           2
```

体系中元素的种类数。有多少种元素就填几。**必须与 STRU 文件中的元素种数一致。**

### pseudo_dir / orbital_dir

```
pseudo_dir      ./
orbital_dir     ./
```

赝势文件和轨道文件所在的目录。可以使用相对路径或绝对路径。`./` 表示当前目录。使用 LCAO 基组时，`orbital_dir` 必须设置。

### calculation

```
calculation     scf
```

最重要的参数之一，决定计算类型：

| 值 | 含义 |
|----|------|
| `scf` | 自洽电子结构计算，获取总能量和基态电荷密度 |
| `nscf` | 非自洽计算，用于能带/DOS（需已有电荷密度）|
| `relax` | 固定晶胞，优化原子位置 |
| `cell-relax` | 同时优化晶胞参数和原子位置 |
| `md` | 分子动力学模拟 |

### basis_type

```
basis_type      pw
```

基组类型：
- `pw`：平面波，通用性强，精度基准
- `lcao`：数值原子轨道，速度更快，需额外提供轨道文件

## 2.3 精度控制参数

### ecutwfc

```
ecutwfc         50
```

平面波截断能，单位 **Rydberg（Ry）**。控制平面波基组的大小——截断能越大，基组越完备，计算越精确，但计算量也越大。

- 默认值：50 Ry
- 典型范围：50–100 Ry（轻元素偏低，重元素或含 d/f 轨道偏高）
- **即使使用 LCAO 基组，也需要设置** `ecutwfc`（局部赝势和力的计算仍用平面波）

### scf_thr

```
scf_thr         1e-6
```

SCF 自洽迭代的收敛阈值，单位 **Ry**。当相邻两步迭代的电荷密度差小于该值时，认为收敛。

- 常规计算：`1e-6`（默认）
- 高精度计算（如力、应力）：`1e-7` 或更小
- 快速预判：`1e-4`（不推荐用于最终结果）

### scf_nmax

```
scf_nmax        100
```

SCF 迭代的最大步数，默认 100。若体系较难收敛（如磁性金属），可适当增大至 200–300。

## 2.4 常用功能参数

### smearing_method / smearing_sigma

展宽（smearing）方法，用于处理费米面附近的电子占据，改善 SCF 收敛性：

```
smearing_method     gauss
smearing_sigma      0.01
```

`smearing_method` 可选值：

| 值 | 适用场景 |
|----|----------|
| `gauss` | 绝缘体、半导体（默认）|
| `mp` | 金属（推荐）|
| `mp2` | 金属（高阶 Methfessel-Paxton）|
| `fd` | 有限温度模拟（Fermi-Dirac）|
| `fixed` | 绝缘体，固定占据数 |

`smearing_sigma` 为展宽宽度，单位 Ry，默认 0.015 Ry：
- 金属：0.005–0.02 Ry
- 绝缘体：0.01 Ry（实际影响不大）

### symmetry

```
symmetry        1
```

是否使用晶体对称性：
- `1`：打开对称性（加速计算，默认）
- `0`：关闭晶体对称性，保留时间反演对称性
- `-1`：关闭所有对称性

做能带计算（nscf）时通常设为 `0`。

### cal_force / cal_stress

```
cal_force       1
cal_stress      1
```

是否计算原子受力和晶胞应力。做结构优化时必须开启（relax/cell-relax 会自动开启）。

### nspin

```
nspin           1
```

自旋极化设置：
- `1`：无自旋极化（默认）
- `2`：共线自旋极化（铁磁/反铁磁体系）
- `4`：非共线磁矩（需设置 `noncolin 1`）

## 2.5 结构优化专用参数

使用 `calculation relax` 或 `calculation cell-relax` 时需要额外设置：

### force_thr_ev

```
force_thr_ev    0.01
```

原子受力的收敛阈值，单位 **eV/Å**。所有原子中的最大受力小于该值时认为力收敛。

- 常规优化：0.01–0.05 eV/Å
- 高精度：0.005 eV/Å 或更小

### stress_thr

```
stress_thr      0.5
```

晶胞应力的收敛阈值，单位 **kBar**（`cell-relax` 专用）。

- 常规优化：0.5–5 kBar
- 高精度：0.1 kBar

### relax_nmax

```
relax_nmax      100
```

结构优化的最大离子迭代步数，默认 50。复杂体系可增大至 100–200。

### out_stru

```
out_stru        1
```

是否在每个离子步后输出 STRU 文件，默认 0。设为 1 时，每步优化结果保存为 `OUT.suffix/STRU_ION_D`，便于追踪优化过程。

## 2.6 输出控制参数

### out_chg

```
out_chg         1
```

是否输出电荷密度文件（`SPIN1_CHG.cube` 等）。做能带计算时需先用 SCF 输出电荷密度，再用 NSCF 读入。

### init_chg

```
init_chg        atomic
```

初始电荷密度的来源：
- `atomic`：从原子电荷密度叠加（默认）
- `file`：从已有的 `SPIN1_CHG` 等文件读入（NSCF 时使用）

## 2.7 完整示例

### 示例1：Si 晶体 SCF 计算（PW 基组）

```
INPUT_PARAMETERS

# 通用参数
suffix          Si
ntype           1
pseudo_dir      ./
calculation     scf
basis_type      pw

# 精度控制
ecutwfc         50          # 单位：Ry
scf_thr         1e-6        # 单位：Ry
scf_nmax        100

# 展宽（Si 是半导体）
smearing_method gauss
smearing_sigma  0.01

# 对称性
symmetry        1
```

### 示例2：MgO 晶体结构优化（LCAO 基组）

```
INPUT_PARAMETERS

# 通用参数
suffix          MgO
ntype           2
pseudo_dir      ./
orbital_dir     ./
calculation     cell-relax
basis_type      lcao

# 精度控制
ecutwfc         100         # 单位：Ry
scf_thr         1e-6
scf_nmax        100

# 展宽（MgO 是绝缘体）
smearing_method gauss
smearing_sigma  0.01

# 对称性
symmetry        1

# 结构优化收敛判据
force_thr_ev    0.01        # 单位：eV/Å
stress_thr      0.5         # 单位：kBar
relax_nmax      100
out_stru        1
```

> **参数设置原则：** 参数不是越多越好。对于大多数场景，默认值已经是合理选择。只设置你明确知道含义的参数。完整的参数列表见官方文档：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html
