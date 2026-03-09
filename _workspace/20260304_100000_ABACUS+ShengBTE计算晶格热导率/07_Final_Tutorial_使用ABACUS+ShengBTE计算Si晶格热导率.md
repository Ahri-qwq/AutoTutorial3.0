---
title: "使用 ABACUS+ShengBTE 计算 Si 的晶格热导率"
author: "AutoTutorial 3.0"
date: "2026-03-04"
topic: "晶格热导率"
task_type: "C"
has_case: true
case_material: "金刚石结构 Si，LCAO 基组"
word_count: 631
chapters: 8
orbital_fix: "Si_gga_7au → Si_gga_8au_100Ry_2s2p1d.orb（自动修正）"
---

# 使用 ABACUS+ShengBTE 计算 Si 的晶格热导率

**作者：AutoTutorial 3.0**

**最后更新：2026-03-04**

---

# 第一章：概述

本教程以金刚石结构 Si 为例，演示从 ABACUS 第一性原理计算出发，
结合 Phonopy、thirdorder 和 ShengBTE，完整地计算晶格热导率。
计算采用数值原子轨道基组（LCAO），并在结尾给出平面波（PW）基组的参数差异。

整体流程如下：

```
  ABACUS (relax)
        │ 优化后 STRU
        ▼
  Phonopy -d                        thirdorder sow
        │ 超胞构型 ×1                      │ 微扰构型 ×40
        ▼                                 ▼
  ABACUS SCF                     ABACUS SCF ×40
        │ 原子受力                         │ 原子受力
        ▼                                 ▼
  phonopy -f                       aba2vasp.py
        │ FORCE_SETS                       │ vasprun.xml ×40
        ▼                                 ▼
  phonopy band.conf                thirdorder reap
        │ FORCE_CONSTANTS                  │
      au2si.py                            │
        │                                 │
        ▼                                 ▼
  FORCE_CONSTANTS_2ND         FORCE_CONSTANTS_3RD
                    │                 │
                    └────────┬────────┘
                             ▼
                         ShengBTE
                             │
                             ▼
                       κ(T)  [W/(m·K)]
```

体系信息：
- 材料：Si，金刚石结构（Fd-3m），2 原子原胞
- 基组：LCAO（主），PW（对比）
- 赝势：`Si_ONCV_PBE-1.0.upf`（模守恒，GGA-PBE）
- 轨道文件：`Si_gga_8au_100Ry_2s2p1d.orb`（DZP，8 au 截断，100 Ry 能量截断）

---

# 第二章：物理背景

## 2.1 声子与晶格热导率

固体中的热量主要通过声子传导。声子是晶格振动的量子，可以理解为原子集体运动的能量包。存在温度梯度时，高温端声子密度更高，声子向低温端扩散，形成热流。

晶格热导率 κ 可以写成各声子模式贡献之和：

$$\kappa = \frac{1}{V} \sum_{\mathbf{q},s} C_{\mathbf{q}s}\, v_{\mathbf{q}s}^2\, \tau_{\mathbf{q}s}$$

其中：
- $C_{\mathbf{q}s}$：声子模式 (q, s) 的热容（由声子频率决定）
- $v_{\mathbf{q}s}$：声子群速度（声子色散曲线的斜率）
- $\tau_{\mathbf{q}s}$：声子寿命（散射弛豫时间）

**计算 κ 需要同时知道声子频率（决定 C 和 v）和声子寿命（τ）。**
这两项分别由二阶和三阶力常数提供。

## 2.2 二阶与三阶力常数

**二阶力常数**（2nd IFC）描述谐振相互作用：一个原子偏离平衡位置时，对周围原子产生多大的回复力。
由二阶力常数构造动力学矩阵，对角化后得到声子色散关系 ω(q) 和群速度 v(q)。

**三阶力常数**（3rd IFC）描述非谐振相互作用：两个原子同时偏离平衡位置时，力的非线性部分。
它决定了声子-声子三声子散射过程，直接控制声子寿命：

$$\tau^{-1} \propto \left|\Phi^{(3)}_{ijk}\right|^2$$

其中 $\Phi^{(3)}_{ijk}$ 是三阶力常数张量的矩阵元。

简言之：**二阶力常数 → 声子频率和群速度；三阶力常数 → 声子散射和寿命。**

## 2.3 有限位移法

计算力常数的标准方法是**有限位移法**（Finite Displacement Method）：
将超胞中的原子沿特定方向位移一个小量（通常约 0.01 Å），用第一性原理计算所有原子受到的力，再通过数值差分提取力常数。

- **二阶力常数**：每次位移一个原子，Phonopy 负责生成构型和从受力中提取力常数
- **三阶力常数**：每次同时位移两个原子（情况更复杂），thirdorder 负责生成构型和提取力常数

由于超胞尺寸决定了力常数的截断范围，超胞越大，力常数描述越完整，热导率结果越收敛，但计算代价也越高。

## 2.4 ShengBTE 的角色

ShengBTE 是求解声子玻尔兹曼输运方程（BTE）的专用软件，承担最后一步计算。
它读取二阶和三阶力常数，在布里渊区密网格上：
1. 构造声子色散并计算群速度
2. 计算声子-声子散射矩阵元（三声子过程）
3. 考虑同位素散射等其他散射机制
4. 求解线性化 BTE，得到各温度下的热导率张量 κ(T)

ShengBTE 的三个必要输入文件：
- `FORCE_CONSTANTS_2ND`：单位 eV/Å²
- `FORCE_CONSTANTS_3RD`：三阶力常数
- `CONTROL`：体系晶体结构信息和计算参数

---

# 第三章：前置准备

## 3.1 案例文件

案例文件可从 Gitee 下载：

```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
```

进入 `examples/interface_ShengBTE/LCAO/` 目录，结构如下：

```
LCAO/
├── 2nd/              # 二阶力常数计算
│   ├── STRU          # relax 后的结构文件
│   ├── INPUT         # ABACUS 输入文件
│   ├── KPT           # k 点设置
│   ├── setting.conf  # Phonopy 超胞设置
│   ├── band.conf     # Phonopy 声子谱和力常数设置
│   └── au2si.py      # 单位转换脚本
├── 3rd/              # 三阶力常数计算
│   ├── POSCAR        # 由 STRU 转换而来（已提供）
│   ├── pos2stru.py   # POSCAR→STRU 转换脚本（调用 ASE）
│   ├── run_stru.sh   # 批量提交 SCF 的脚本
│   └── aba2vasp.py   # 将 ABACUS 受力封装为 vasprun.xml 的脚本
└── shengbte/         # ShengBTE 计算
    ├── CONTROL       # ShengBTE 参数文件
    ├── FORCE_CONSTANTS_2ND   # 参考文件（已提供）
    ├── FORCE_CONSTANTS_3RD   # 参考文件（已提供）
    └── Ref/          # 参考计算结果
```

## 3.2 软件依赖

| 软件 | 说明 |
|------|------|
| ABACUS ≥3.2.0 | 第一性原理计算，提供结构优化和 SCF 受力 |
| Phonopy ≥2.19.1 | 二阶力常数计算，已支持 ABACUS 接口 |
| ASE | 原子结构格式转换（pos2stru.py 内部调用） |
| thirdorder | 三阶力常数，依赖 VASP/QE 格式的结构和受力输入 |
| ShengBTE | BTE 求解，输出 κ(T) |

相关文档：
- Phonopy 接口：https://abacus.deepmodeling.com/en/latest/advanced/interface/phonopy.html
- ASE 接口：https://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html
- ShengBTE / thirdorder：https://bitbucket.org/sousaw/shengbte/src/master/

---

# 第四章：Step 1 — 计算二阶力常数

进入 `2nd` 目录。本步骤使用 ABACUS + Phonopy 计算二阶力常数，最终得到 `FORCE_CONSTANTS_2ND`。

## 4.1 结构优化

计算力常数前，先对体系进行结构优化（`calculation = relax`），确保原子处于受力平衡的平衡构型。
本例使用 2×2×2 k 点采样，ecut = 100 Ry（注意：实际研究中应测试 k 点收敛性）。

结构优化完成后，得到以下 STRU 文件：

```
ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
0 2.81594778072 2.81594778072 #latvec1
2.81594778072 0 2.81594778072 #latvec2
2.81594778072 2.81594778072 0 #latvec3

ATOMIC_POSITIONS
Direct # direct coordinate

Si #label
0 #magnetism
2 #number of atoms
0.875  0.875  0.875  m  0  0  0
0.125  0.125  0.125  m  0  0  0
```

说明：
- `LATTICE_CONSTANT = 1.88972612546`：单位为 Bohr，此值等于 1 Å（1 Å = 1.88972612546 Bohr），格矢量乘以此缩放因子得到实际格矢长度
- 实际格矢长度：1.88972612546 × 2.81594778072 ≈ 5.32 Bohr ≈ 2.82 Å
- `m 0 0 0` 表示原子坐标固定，为 relax 优化后的最终平衡位置
- Si 的质量（28.0855）在力常数计算中不起作用

## 4.2 Phonopy 产生超胞微扰构型

在 `2nd` 目录下执行：

```bash
phonopy setting.conf --abacus -d
```

其中 `setting.conf` 的内容为：

```
DIM = 2 2 2
ATOM_NAME = Si
```

`DIM = 2 2 2` 表示将 2 原子原胞扩展为 2×2×2 超胞（共 16 个原子）。
Phonopy 会生成需要计算受力的微扰超胞构型文件。

Si 金刚石结构对称性高，只需产生 **1 个**微扰构型 `STRU-001`。

## 4.3 ABACUS SCF 计算受力

对 `STRU-001` 进行单点 SCF 计算，获取原子受力。INPUT 中的关键设置：

```
calculation   scf
cal_force     1        # 必须开启受力输出
stru_file     STRU-001 # 指向 Phonopy 生成的超胞构型
```

**小技巧**：通过 `stru_file` 参数直接指定 Phonopy 生成的构型文件，无需手动重命名。

计算完成后，用以下命令从 ABACUS 输出中提取受力，生成 `FORCE_SETS` 文件：

```bash
phonopy -f OUT.DIA-50/running_scf.log
```

其中 `OUT.DIA-50` 是 ABACUS 的输出目录（目录名由 INPUT 中的 `suffix` 参数决定，
若 `suffix = DIA-50` 则输出目录为 `OUT.DIA-50`）。

## 4.4 计算声子谱和二阶力常数

执行：

```bash
phonopy -p band.conf --abacus
```

`band.conf` 的完整内容为：

```
ATOM_NAME = Si
DIM = 2 2 2
MESH = 8 8 8
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.5  0.625 0.25 0.625, 0.375 0.375 0.75  0.0 0.0 0.0  0.5 0.5 0.5
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
FORCE_CONSTANTS = WRITE
FULL_FORCE_CONSTANTS = .TRUE.
```

参数说明：
- `MESH = 8 8 8`：声子谱积分的 q 点网格
- `PRIMITIVE_AXES`：原始布里渊区的基矢，这里 Si FCC 取标准设置
- `BAND`：高对称路径点，覆盖 Γ→X→U/K→Γ→L 路径（FCC 标准声子谱路径，逗号表示路径断点）
- `FORCE_CONSTANTS = WRITE`：将力常数写入文件
- **`FULL_FORCE_CONSTANTS = .TRUE.`**：输出完整力常数矩阵（必须设置，见注意事项）

命令执行后，Phonopy 输出：
- `band.yaml`：声子谱数据（可用于绘图）
- `FORCE_CONSTANTS`：二阶力常数文件

如需导出 gnuplot 格式的声子谱：

```bash
phonopy-bandplot --gnuplot > pho.dat
```

## 4.5 单位转换：生成 FORCE_CONSTANTS_2ND

ShengBTE 要求 `FORCE_CONSTANTS_2ND` 的单位为 **eV/Å²**，
但 ABACUS 结合 Phonopy 输出的 `FORCE_CONSTANTS` 单位为 **eV/(Å·au)**，其中 1 au = 0.52918 Å。

运行单位转换脚本：

```bash
python au2si.py
```

脚本将 `FORCE_CONSTANTS` 乘以相应的单位换算系数，输出 `FORCE_CONSTANTS_2ND`，即 ShengBTE 所需的格式。

至此，第一步完成。`shengbte/` 目录中已提供参考用的 `FORCE_CONSTANTS_2ND` 文件，可用于对照验证。

---

# 第五章：Step 2 — 计算三阶力常数

进入 `3rd` 目录。本步骤使用 ABACUS + thirdorder 计算三阶力常数，最终得到 `FORCE_CONSTANTS_3RD`。

**背景说明**：thirdorder 目前只支持读取 VASP 和 QE 格式的输入/输出文件，
因此需要将 ABACUS 的结构文件转换为 POSCAR，将输出的受力封装为 vasprun.xml，再交给 thirdorder 处理。

## 5.1 产生微扰构型

首先，将 relax 后的 STRU 文件转换为 POSCAR 格式（`3rd/` 目录下已提供转换好的 POSCAR，
若需自行转换，可使用 ASE 完成）。

运行 thirdorder 产生微扰构型：

```bash
thirdorder_vasp.py sow 2 2 2 -2
```

参数说明：
- `sow`：播种模式，生成微扰超胞构型
- `2 2 2`：超胞扩展倍数（与二阶力常数计算保持一致）
- `-2`：近邻截断参数，负号表示取到第 2 近邻范围内的原子对

此命令为 Si 体系产生 **40 个** POSCAR 文件：`3RD.POSCAR.1`、`3RD.POSCAR.2`、……、`3RD.POSCAR.40`。

将这些 POSCAR 文件转换为 ABACUS 的 STRU 格式：

```bash
python pos2stru.py
```

> **注意**：转换时只能使用 `pos2stru.py`（内部调用 ASE），**不能使用 dpdata**。
> dpdata 在转换时会强制输出下三角晶格矩阵，等效于对晶格做了旋转，
> 导致 ABACUS 计算的受力方向同步旋转，三阶力常数因此完全错误。
> 这是一个容易忽视但影响致命的问题。

## 5.2 批量 SCF 计算

需要对 40 个微扰构型分别进行单点 SCF 计算，获取原子受力。
`run_stru.sh` 脚本会自动创建 `SCF-1/`、`SCF-2/`、……、`SCF-40/` 子目录并批量提交计算。

INPUT 中的关键参数：

```
calculation   scf
cal_force     1      # 必须输出受力
scf_thr       1e-8   # LCAO 基组：收敛阈值至少需要 1e-8
```

**`scf_thr` 对三阶力常数精度至关重要**：受力误差直接传递到力常数误差，进而影响声子寿命和热导率。
LCAO 基组至少需要 `1e-8`；若使用 PW 基组，需要至少 `1e-12`（原因详见第七章注意事项）。

40 个 SCF 任务彼此独立，建议并行提交到集群以节省时间。

## 5.3 封装受力为 vasprun.xml

40 个 SCF 计算完成后，执行：

```bash
python aba2vasp.py
```

该脚本读取每个 `SCF-*/` 中 ABACUS 输出的原子受力，将其封装为 thirdorder 能识别的 vasprun.xml 格式，放置在各 `SCF-*/` 目录中。

生成的 vasprun.xml 格式如下（以某一构型的受力为例）：

```xml
<modeling>
    <calculation>
        <varray name="forces">
            <v>1.865e-05 -0.04644196 -0.00153852</v>
            <v>-1.77e-05 -0.00037715 -0.00149635</v>
            <v>1.973e-05  0.002213   -0.00149461</v>
            <v>-1.976e-05 0.00065303 -0.0014804</v>
            ...（共 16 个原子的受力分量）
        </varray>
    </calculation>
</modeling>
```

每个 `<v>` 标签对应一个原子在 x、y、z 三个方向上的受力（单位：eV/Å）。

## 5.4 提取三阶力常数

执行以下命令，收集所有构型的受力并提取三阶力常数：

```bash
find SCF-* -name vasprun.xml | sort -n | thirdorder_vasp.py reap 2 2 2 -2
```

命令说明：
- `find SCF-* -name vasprun.xml | sort -n`：按编号顺序找到所有 vasprun.xml 文件
- `thirdorder_vasp.py reap`：收割模式，从受力数据中提取三阶力常数
- `2 2 2 -2`：与 `sow` 步骤使用相同的参数（超胞倍数和截断近邻数）

运行完成后生成 `FORCE_CONSTANTS_3RD`。
`shengbte/` 目录中已提供参考文件，可用于对照。

至此，两种力常数文件均已准备好，可以进入最后一步。

---

# 第六章：Step 3 — 运行 ShengBTE

进入 `shengbte/` 目录，确认以下三个文件已就位：
- `FORCE_CONSTANTS_2ND`（来自 Step 1）
- `FORCE_CONSTANTS_3RD`（来自 Step 2）
- `CONTROL`（ShengBTE 参数文件，已提供）

## 6.1 CONTROL 文件

`CONTROL` 是 ShengBTE 的核心参数文件，使用 Fortran namelist 格式，包含四个段：

```fortran
&allocations
    nelements=1
    natoms=2
    ngrid(:)=10 10 10
&end
&crystal
    lfactor=0.100000
    lattvec(:,1)=0 2.81594778072 2.81594778072
    lattvec(:,2)=2.81594778072 0 2.81594778072
    lattvec(:,3)=2.81594778072 2.81594778072 0
    elements="Si"
    types=1 1
    positions(:,1)=0.8750000000000000  0.8750000000000000  0.8750000000000000
    positions(:,2)=0.1250000000000000  0.1250000000000000  0.1250000000000000
    scell(:)=2 2 2
&end
&parameters
    !T=300,
    T_min=200
    T_max=500
    T_step=50
    scalebroad=1.0
&end
&flags
    !espresso=.true.
    nonanalytic=.true.,
    isotopes=.true.
&end
```

各段参数说明：

**&allocations**
| 参数 | 值 | 说明 |
|------|----|------|
| `nelements` | 1 | 元素种类数（Si 只有 1 种） |
| `natoms` | 2 | 原胞中原子数 |
| `ngrid(:)` | 10 10 10 | 布里渊区 q 点积分网格，越大精度越高 |

**&crystal**
| 参数 | 说明 |
|------|------|
| `lfactor` | 格矢量长度单位因子（0.1 nm = 1 Å） |
| `lattvec(:,i)` | 三个格矢量，单位为 lfactor × nm = Å，与 STRU 中 LATTICE_VECTORS 对应 |
| `elements` | 元素名称 |
| `types` | 每个原子的元素类型编号（两个 Si 均为类型 1） |
| `positions(:,i)` | 每个原子的分数坐标，与 STRU 中的原子坐标一致 |
| `scell(:)` | 计算力常数时使用的超胞大小，**必须与 Step 2 保持一致** |

**&parameters**
| 参数 | 值 | 说明 |
|------|----|------|
| `T_min` | 200 | 计算温度下限（K） |
| `T_max` | 500 | 计算温度上限（K） |
| `T_step` | 50 | 温度步长（K），本例计算 200、250、300、350、400、450、500 K |
| `scalebroad` | 1.0 | 声子线宽展宽因子，控制散射积分的展宽精度 |

（`!T=300` 行以感叹号开头表示注释，ShengBTE 不会读取该行。）

**&flags**
| 参数 | 值 | 说明 |
|------|----|------|
| `nonanalytic` | .true. | 开启非解析修正（对极性材料有较大影响，Si 中效果较小） |
| `isotopes` | .true. | 考虑天然同位素散射 |

## 6.2 运行 ShengBTE

```bash
mpirun -n 10 ShengBTE
```

ShengBTE 运行期间会输出声子对称性信息、q 点数量等日志，例如：

```
Info: symmetry group Fd-3m detected
Info: 48 symmetry operations
Info: Ntot = 1000
Info: Nlist = 47
Info: about to obtain the spectrum
Info: expecting Phonopy 2nd-order format
```

计算完成后，热导率结果保存在 `BTE.kappa_*` 系列文件中。

## 6.3 计算结果

本教程参数（2×2×2 超胞，2×2×2 k 点，ngrid=10×10×10）下，Si 的计算结果如下：

| 基组 | 300 K 热导率 | 实验值（300 K） |
|------|-------------|----------------|
| LCAO | ~100 W/(m·K) | ~150 W/(m·K) |
| PW   | ~100 W/(m·K) | ~150 W/(m·K) |

LCAO 和 PW 基组的结果高度一致，说明两种基组对 Si 体系均能给出可靠的力常数。

计算值低于实验值约 30%，原因是本教程为演示目的使用了较小的超胞（2×2×2）和较稀的 k 点（2×2×2），导致力常数截断误差和 k 点采样误差同时存在。实际科研中需要系统地测试这两个参数（见第八章）。

对于 PW 基组的计算，流程与 LCAO 完全相同，参数差异见第七章注意事项。

---

# 第七章：关键注意事项

整个计算流程中有几处容易出错，汇总如下。

## 7.1 必须设置 FULL_FORCE_CONSTANTS = .TRUE.

`band.conf` 中必须加入这一行：

```
FULL_FORCE_CONSTANTS = .TRUE.
```

不设置时，Phonopy 默认只输出精简形式的力常数文件（仅含不等价部分），ShengBTE 读取时会因格式不匹配而直接报错。这是一个不设置时不会有任何警告、但结果完全无法运行的问题。

## 7.2 二阶力常数的单位转换不能跳过

ABACUS 结合 Phonopy 输出的 `FORCE_CONSTANTS` 单位为 **eV/(Å·au)**，
而 ShengBTE 要求 `FORCE_CONSTANTS_2ND` 的单位为 **eV/Å²**（1 au = 0.52918 Å）。
如果跳过 `au2si.py` 的单位转换步骤，直接将 `FORCE_CONSTANTS` 重命名为 `FORCE_CONSTANTS_2ND`，
ShengBTE 的计算不会报错，但热导率结果的数量级完全错误。

## 7.3 结构格式转换不能用 dpdata

将 `3RD.POSCAR.*` 转换为 STRU 时，只能使用 `pos2stru.py`（调用 ASE），**不能用 dpdata**。

dpdata 在转换时会将晶格矩阵强制变换为下三角形式，这等效于对整个晶格做了一次旋转。旋转后，晶体结构的物理是等价的，但 ABACUS 计算出的原子受力方向也会对应旋转——而 thirdorder 在 `reap` 阶段并不知道这个旋转，最终提取出的三阶力常数方向完全错误。

## 7.4 SCF 收敛阈值须足够严格

计算三阶力常数时，受力精度直接影响力常数精度，进而影响声子寿命和热导率。

| 基组 | `scf_thr` 最低要求 |
|------|------------------|
| LCAO | 1e-8 |
| PW   | 1e-12 |

PW 基组对收敛阈值要求更严格（需要 1e-12），是因为平面波基组下哈密顿矩阵元的数值精度对 scf_thr 更敏感。如果不满足此要求，受力误差会使三阶力常数出现较大噪声，导致热导率结果不收敛。

## 7.5 CONTROL 中的 scell 须与 thirdorder 保持一致

CONTROL 文件中的 `scell(:)=2 2 2` 必须与 thirdorder 中使用的超胞参数（`thirdorder_vasp.py sow 2 2 2 -2` 中的 `2 2 2`）完全一致。
如果不一致，ShengBTE 在读取 `FORCE_CONSTANTS_3RD` 时会出现维度不匹配的错误，或者无声地给出错误结果。

---

# 第八章：进阶与展望

## 8.1 收敛性测试

本教程使用的参数（2×2×2 超胞、2×2×2 k 点、ngrid=10×10×10）是面向演示的最小参数集，计算速度快但精度有限（300 K 结果偏低约 30%）。实际研究中需要对以下三个参数分别测试收敛性：

- **超胞大小**：逐步增大（2×2×2 → 3×3×3 → 4×4×4），直到热导率不再随超胞增大而显著变化
- **k 点采样**：增大原胞 SCF 和 relax 的 k 点密度，确保电子结构收敛
- **ShengBTE 的 ngrid**：增大 q 点积分网格（如 15×15×15），确保热导率积分收敛

以 Si 为例，使用 4×4×4 超胞和更密的 k 点可将 300 K 热导率提高到接近实验值（~150 W/(m·K)）。

## 8.2 扩展到其他体系

本教程展示的 ABACUS+Phonopy+thirdorder+ShengBTE 工作流适用于大多数半导体和绝缘体材料。
对于有强非谐效应的材料（如热电材料 PbTe、SnSe），流程相同，但通常需要更大的超胞和更严格的 `scf_thr`。
对于有极性的离子晶体（如 GaN、BN），建议在 CONTROL 中开启非解析修正（`nonanalytic=.true.`）并提供 Born 有效电荷数据，以正确处理 LO-TO 劈裂对热导率的贡献。

---

# 附录

## 参考资料

- ABACUS 官方文档：https://abacus.deepmodeling.com/
- ABACUS+ShengBTE 官方教程：https://mcresearch.github.io/abacus-user-guide/abacus-shengbte.html
- 案例文件（Gitee）：https://gitee.com/mcresearch/abacus-user-guide/tree/master/examples/interface_ShengBTE
- Phonopy 文档：https://phonopy.github.io/phonopy/
- ShengBTE / thirdorder：https://bitbucket.org/sousaw/shengbte/src/master/

## 常见问题

**Q：ShengBTE 运行报错 "expecting Phonopy 2nd-order format"**
A：检查 FORCE_CONSTANTS_2ND 是否由 `au2si.py` 生成（不能直接用 Phonopy 的 FORCE_CONSTANTS 重命名）；另外确认 `band.conf` 中设置了 `FULL_FORCE_CONSTANTS = .TRUE.`。

**Q：thirdorder reap 步骤报错**
A：确认 `sort -n` 排序正确，所有 SCF 计算已完成，各 SCF-* 目录中均有 vasprun.xml。

**Q：LCAO 计算受力为零或异常**
A：检查 INPUT 中是否设置了 `cal_force = 1`；同时检查 `stru_file` 是否正确指向 Phonopy 生成的构型文件。

**Q：热导率与实验值偏差较大**
A：检查超胞大小和 k 点采样是否收敛；同时确认 `scf_thr` 满足最低要求（LCAO：1e-8，PW：1e-12）。

## 快速命令参考

```bash
# === Step 1：二阶力常数 ===
# 生成超胞微扰构型
phonopy setting.conf --abacus -d
# 提取受力（SCF 完成后）
phonopy -f OUT.DIA-50/running_scf.log
# 计算声子谱和力常数
phonopy -p band.conf --abacus
# 单位转换
python au2si.py

# === Step 2：三阶力常数 ===
# 生成 40 个微扰构型
thirdorder_vasp.py sow 2 2 2 -2
# POSCAR→STRU 转换
python pos2stru.py
# 批量提交 SCF（需修改 run_stru.sh 适配集群）
bash run_stru.sh
# 封装受力（所有 SCF 完成后）
python aba2vasp.py
# 提取三阶力常数
find SCF-* -name vasprun.xml | sort -n | thirdorder_vasp.py reap 2 2 2 -2

# === Step 3：ShengBTE ===
mpirun -n 10 ShengBTE
```
