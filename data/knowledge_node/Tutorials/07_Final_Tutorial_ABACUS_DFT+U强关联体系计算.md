---
title: "ABACUS 使用教程｜DFT+U 强关联体系计算"
author: "AutoTutorial 3.0"
date: "2026-02-28"
topic: "DFT+U强关联体系计算设置"
task_type: "C"
has_case: true
case_system: "反铁磁NiO（LCAO, SCF+DFT+U+omc）"
word_count: ~5000
chapters: 5
---

# ABACUS 使用教程｜DFT+U 强关联体系计算

作者：AutoTutorial 3.0

---

## 1. 引言

过渡金属氧化物（如 NiO、MnO、FeO）是凝聚态物理和材料科学中的经典研究对象。用常规的 LDA 或 GGA 泛函计算这类体系时，往往得到错误的基态：NiO 本应是绝缘体，用 GGA 却常给出金属态，能隙严重低估，磁矩也不准确。

这个失败的根源在于 GGA 无法正确描述 Ni 的 d 轨道中强烈的电子关联效应。DFT+U 方法通过在 d/f 轨道上叠加 Hubbard 模型修正，以很低的计算代价改善这一问题。

本教程以反铁磁 NiO 为算例，展示 ABACUS LCAO 基组下 DFT+U 的完整计算流程：
- 设置反铁磁初始磁矩
- 配置 DFT+U 输入参数
- 提取并验证总能量、能隙、磁矩
- 使用 occupation matrix control (omc) 功能

> **注意**：ABACUS 的 DFT+U 功能**仅支持 LCAO 基组**（`basis_type lcao`），平面波基组下不可用。
## 2. 体系与文件准备

### 2.1 算例体系：反铁磁 NiO

NiO 具有岩盐结构，原胞含 2 个 Ni 原子和 2 个 O 原子。在反铁磁基态下，两个 Ni 原子的磁矩大小相等、方向相反，体系总磁矩为零。

在 ABACUS 中，需要将**磁矩方向不同的同种原子定义为不同的 atomic species**，才能分别设置初始磁矩。本算例将两个 Ni 分别定义为 `Ni1`（磁矩 +2.0）和 `Ni2`（磁矩 -2.0），使用完全相同的赝势和轨道文件。

### 2.2 STRU 文件

STRU 文件中，`ATOMIC_SPECIES` 部分定义三类原子：

```
# STRU 文件（关键部分）

ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf    # 第一种Ni：赝势文件
Ni2 58.693 Ni_ONCV_PBE-1.0.upf    # 第二种Ni：同一个赝势文件
O   15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Ni_gga_9au_100Ry_4s2p2d1f.orb     # Ni1轨道文件（双d基组）
Ni_gga_9au_100Ry_4s2p2d1f.orb     # Ni2共用同一轨道文件
O_gga_7au_100Ry_2s2p1d.orb
```

在 `ATOMIC_POSITIONS` 部分，每种 species 下方跟随初始磁矩值：

```
ATOMIC_POSITIONS
Direct

Ni1
2.0              # 初始磁矩（Bohr mag），spin-up
1
0.000000000  0.000000000  0.000000000  0  0  0

Ni2
-2.0             # 初始磁矩（Bohr mag），spin-down，方向相反
1
0.500000000  0.500000000  0.000000000  0  0  0

O
0.0
2
0.500000000  0.000000000  0.500000000  0  0  0
0.000000000  0.500000000  0.500000000  0  0  0
```

> 上述坐标为示意（岩盐结构的典型直角坐标），完整 STRU 文件以算例数据集为准。

> **技巧**：把不同磁矩的同种原子当作不同 species，是 ABACUS 设置共线反铁磁的标准做法。两种 Ni 可以共用同一套赝势和轨道文件，不影响计算。

### 2.3 INPUT 文件

完整的 INPUT 文件如下。DFT+U 相关参数集中在文件末尾。

```
# INPUT 文件

INPUT_PARAMETERS
#Parameters (General)
suffix                  NiO
calculation             scf
esolver_type            ksdft
symmetry                0

#Parameters (Basis)
basis_type              lcao
ecutwfc                 100

#Parameters (Spin)
nspin                   2

#Parameters (Iteration)
scf_thr                 1e-7
scf_nmax                200

#Parameters (Smearing)
smearing_method         gauss    # 典型设置，具体值以算例为准
smearing_sigma          0.01

#Parameters (Mixing)
mixing_type             broyden  # 典型设置
mixing_beta             0.4

#Parameters (Output)
out_bandgap             1        # 输出能隙信息
out_mul                 1        # 输出 Mulliken 电荷分析（含原子磁矩）
out_chg                 1        # 输出电荷密度（同时输出 onsite.dm）

#Parameters (DFT+U)
dft_plus_u              1        # DFT+U 总开关：1=开，0=关
orbital_corr            2 2 -1   # 各 species 施加+U的l量子数：2=d轨道，-1=不施加
hubbard_u               5.0 5.0 0.0  # 各 species 的U值（eV）
```

#### DFT+U 三个核心参数详解

**`dft_plus_u`**
- 总控制开关
- `1`：开启 DFT+U；`0`：关闭（其他+U参数被忽略）

**`orbital_corr`**
- 数组，长度等于 `ntype`（atomic species 总数，本例为 3）
- 每个数字指定对应 species 施加+U修正的角量子数：`0`=s，`1`=p，`2`=d，`3`=f，`-1`=不施加
- 本例：`2 2 -1` 表示对 Ni1 的 d 轨道和 Ni2 的 d 轨道各施加+U，O 不施加

**`hubbard_u`**
- 数组，长度与 `orbital_corr` 相同
- 单位：eV
- 本例：对 Ni1 和 Ni2 各施加 5.0 eV，O 为 0.0 eV

> **注意事项**：
> - `orbital_corr` 和 `hubbard_u` 的数组长度必须等于 `ntype`，缺少或多余都会报错
> - U 值的选取取决于体系，NiO 的 Ni d 轨道常用 4~6 eV，本例采用 5 eV
> - `symmetry` 需要设为 0，否则可能因对称性操作影响磁结构

### 2.4 赝势与轨道文件

| 文件 | 说明 |
|------|------|
| `Ni_ONCV_PBE-1.0.upf` | Ni 的 ONCV 模守恒赝势（PBE泛函），Ni1 和 Ni2 共用 |
| `O_ONCV_PBE-1.0.upf` | O 的 ONCV 模守恒赝势 |
| `Ni_gga_9au_100Ry_4s2p2d1f.orb` | Ni 的数值原子轨道（4s2p2d1f，双d基组） |
| `O_gga_7au_100Ry_2s2p1d.orb` | O 的数值原子轨道 |

轨道文件中，Ni 只有一个 d 基组（zeta=0），这在后续分析 occupation matrix 时会用到。

---

## 3. 运行计算

### 3.1 运行 SCF 计算

在 `ABACUS_DFT+U/` 目录下执行：

```bash
# 进入工作目录
cd ABACUS_DFT+U

# 使用 8 个 MPI 进程，单线程运行
OMP_NUM_THREADS=1 mpirun -n 8 abacus
```

根据机器配置调整 MPI 进程数。计算在普通工作站上通常几分钟内完成。

### 3.2 检查计算收敛

计算完成后，检查 `OUT.NiO/running_scf.log`。每步 SCF 迭代会输出一行，包含总磁矩（TMAG）、绝对磁矩（AMAG）、总能量（ETOT）和收敛指标（EDIFF、DRHO）：

```
 ITER   TMAG      AMAG      ETOT(eV)          EDIFF(eV)     DRHO       TIME(s)
 GE1    ...       ...       ...              0.000e+00      ...        ...
 GE2    ...       ...       ...              ...            ...        ...
 ...
```

关键检查项：
- `ETOT` 趋于稳定：相邻步差值小于 `scf_thr`（本例为 1e-7 eV）
- `TMAG` 始终为 0：确认体系保持反铁磁态（总磁矩为零）
- `AMAG` 收敛到约 3.35：绝对磁矩反映两个 Ni 的磁矩大小

若 SCF 不收敛，可尝试减小 `mixing_beta`（如 0.2）或增大 `scf_nmax`。

---

## 4. 结果分析

计算完成后，所有输出文件位于 `OUT.NiO/` 目录下。

### 4.1 总能量

在 `running_scf.log` 中搜索 `FINAL_ETOT`：

```bash
grep "FINAL_ETOT" OUT.NiO/running_scf.log
```

预期输出：

```
 --------------------------------------------
 !FINAL_ETOT_IS -9255.7279034240546025 eV
 --------------------------------------------
```

### 4.2 能隙

INPUT 中设置了 `out_bandgap 1`，在 `running_scf.log` 中搜索 `E_bandgap`，取最后一次出现的值：

```bash
grep "E_bandgap" OUT.NiO/running_scf.log | tail -1
```

预期输出：

```
E_bandgap               +0.205369322748               +2.794192983776
```

两列分别是 spin-up 和 spin-down 的能隙（单位：eV）。NiO 是绝缘体，DFT+U 给出了非零能隙。若使用纯 GGA，NiO 通常给出 0 能隙（金属态），这正是 DFT+U 修正的意义所在。

> spin-up 能隙（0.21 eV）和 spin-down 能隙（2.79 eV）不对称，反映了两种自旋通道中能带结构的差异。

### 4.3 磁矩

**总磁矩与绝对磁矩**

搜索 `absolute magnetism`，取最后一次出现的值：

```bash
grep -A1 "total magnetism" OUT.NiO/running_scf.log | tail -4
```

预期输出：

```
      total magnetism (Bohr mag/cell) = 0.00000000
   absolute magnetism (Bohr mag/cell) = 3.35321634
```

- 总磁矩为 0：验证了体系是反铁磁态
- 绝对磁矩 3.35：两个 Ni 原子磁矩的绝对值之和

**原子磁矩（Mulliken 分析）**

INPUT 中设置了 `out_mul 1`，生成 `OUT.NiO/Mulliken.txt`。搜索 `Magnetism`：

```bash
grep "Total Magnetism" OUT.NiO/Mulliken.txt
```

预期输出：

```
Total Magnetism on atom  Ni1           1.8268646
Total Magnetism on atom  Ni2          -1.8268646
Total Magnetism on atom  O      -3.6718263e-13
Total Magnetism on atom  O       1.7330755e-13
```

Ni1 和 Ni2 的磁矩大小相等（约 1.83 μ_B）、方向相反，O 原子的磁矩几乎为零。这是反铁磁 NiO 的预期结果。

### 4.4 Occupation Matrix 的读取

DFT+U 计算会在 `running_scf.log` 中的每一步输出 occupation matrix。搜索以 `L(S)DA+U` 开头的块：

```bash
grep -A 50 "L(S)DA+U" OUT.NiO/running_scf.log | head -60
```

该块的结构如下：

```
L(S)DA+U:
atom_type=0  L=2  chi=0    U=5ev    # Ni1：d轨道，U=5 eV
atom_type=1  L=2  chi=0    U=5ev    # Ni2：d轨道，U=5 eV

atoms  0          # 第一个Ni原子（Ni1）
L  2              # d 轨道（l=2）
zeta  0           # 第一个（也是唯一一个）d基组
spin  0           # spin-up 的 occupation matrix（5×5）
  0.9xx  ...
  ...             # 5行5列
spin  1           # spin-down 的 occupation matrix
  0.0xx  ...
  ...

atoms  1          # 第二个Ni原子（Ni2）
L  2
zeta  0
spin  0
  ...
spin  1
  ...
```

occupation matrix 的物理含义：每个元素 n_{mm'} 表示 d 轨道磁量子数 m 和 m' 之间的占据，对角元代表各 d 分量的占据数（0~1）。Ni1 的 spin-up 分量占据较满，Ni2 的 spin-down 分量占据较满，体现了两者的磁矩相反。

**`onsite.dm` 文件**

计算结束后（因设置了 `out_chg 1`），在 `OUT.NiO/` 中生成 `onsite.dm` 文件，保存最后一步的 occupation matrix，格式与上述 log 中的输出相同。该文件在下一节 omc 演示中使用。

---

## 5. 进阶：Occupation Matrix Control (omc)

### 5.1 功能说明

在 DFT+U 计算中，SCF 迭代过程中 occupation matrix 会随着电子密度一起更新。对于复杂体系，不同的初始磁矩可能导致收敛到不同的局域极小（"多解问题"）。occupation matrix control（omc）提供了一种手动控制 occupation matrix 初始值的方法。

INPUT 中的 `omc` 参数有三种取值：

| omc 值 | 行为 |
|--------|------|
| 0 | 标准 DFT+U，不使用外部 occupation matrix，从原子初始值开始迭代 |
| 1 | 第一步读入 `initial_onsite.dm`，后续步骤照常更新 |
| 2 | 始终使用 `initial_onsite.dm` 中的 occupation matrix，全程固定，不更新 |

`initial_onsite.dm` 的文件格式与 `OUT.NiO/onsite.dm` 完全相同。

### 5.2 演示：使用 omc=2 固定 occupation matrix

#### 步骤 1：复制收敛后的 onsite.dm

将上一步计算生成的 `onsite.dm` 复制到工作目录，重命名为 `initial_onsite.dm`：

```bash
cp OUT.NiO/onsite.dm ./initial_onsite.dm
```

#### 步骤 2：修改 INPUT 文件

在 INPUT 中添加 `omc` 参数：

```
#Parameters (DFT+U)
dft_plus_u              1
orbital_corr            2 2 -1
hubbard_u               5.0 5.0 0.0
omc                     2        # 固定使用 initial_onsite.dm
```

#### 步骤 3：重新运行计算

```bash
OMP_NUM_THREADS=1 mpirun -n 8 abacus
```

#### 步骤 4：验证结果

使用 omc=2 的计算结果应与 omc=0 的标准 DFT+U 结果完全一致：

```
 !FINAL_ETOT_IS -9255.7279034240546025 eV   # 与标准计算相同
```

在 `running_scf.log` 中，每一步的 occupation matrix 应保持不变（等于 `initial_onsite.dm` 中的值）。

### 5.3 omc 的实际应用场景

- **多解问题**：当体系可能存在多个磁态（如铁磁、反铁磁、亚铁磁），可以通过指定不同的 `initial_onsite.dm` 引导收敛到目标磁态
- **重启计算**：对已知合理磁态的体系，用 omc=1 给出一个好的初始 occupation matrix，加快 SCF 收敛
- **固定轨道占据**：某些场景需要约束特定轨道的占据（如研究激发态），omc=2 可实现全程固定

### 5.4 补充：Yukawa Potential 自动确定 U 值

如果不想手动设置 U 值，ABACUS 还提供了通过 Yukawa potential 自动计算 U 值的方法：

```
yukawa_potential         1        # 开启 Yukawa potential 方法
# yukawa_lambda          X.X     # 可选：手动指定 screening length
```

开启后，程序在 SCF 过程中通过将电子相互作用近似为 Yukawa potential，自动计算各轨道的 U 值。更多细节参见 Qu et al. (2022)。

---

## 附录

### A. DFT+U 参数速查表

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `dft_plus_u` | int | 0 | DFT+U 总开关（1=开，0=关） |
| `orbital_corr` | int 数组 | — | 各 species 施加+U的l量子数（-1=不施加；长度=ntype） |
| `hubbard_u` | float 数组 | — | 各 species 的U值，单位 eV（长度=ntype） |
| `omc` | int | 0 | occupation matrix control（0=标准，1=读入后更新，2=固定） |
| `yukawa_potential` | int | 0 | 自动计算U值开关（1=开） |
| `yukawa_lambda` | float | — | Yukawa screening length（手动设置，可选） |

### B. 常用 grep 命令速查

```bash
# 总能量
grep "FINAL_ETOT" OUT.NiO/running_scf.log

# 能隙（最后一步）
grep "E_bandgap" OUT.NiO/running_scf.log | tail -1

# 总磁矩和绝对磁矩（最后一步）
grep "magnetism" OUT.NiO/running_scf.log | tail -2

# 原子磁矩（Mulliken）
grep "Total Magnetism" OUT.NiO/Mulliken.txt

# occupation matrix（每SCF步）
grep -n "L(S)DA+U" OUT.NiO/running_scf.log
```

### C. 常见问题

**Q：计算报错 "orbital_corr length does not match ntype"**
A：`orbital_corr` 和 `hubbard_u` 的数组长度必须等于 STRU 中 `ATOMIC_SPECIES` 的种类数（ntype）。本例 ntype=3（Ni1/Ni2/O），所以两个数组都需要 3 个数值。

**Q：SCF 不收敛，磁矩振荡**
A：反铁磁体系 SCF 有时不稳定。尝试：减小 `mixing_beta`（0.2~0.3）、设置 `mixing_gg0`（约 1.0）或增大 `scf_nmax`。

**Q：想算铁磁态，如何设置？**
A：在 STRU 中将所有 Ni 原子设为同一 species（或两种 species 但磁矩同号），`orbital_corr` 和 `hubbard_u` 也只需对 Ni 的 species 设置。

**Q：只有一种 Ni 的情况如何设置 orbital_corr？**
A：如果 STRU 中只有 Ni 和 O（ntype=2），则 `orbital_corr 2 -1`，`hubbard_u 5.0 0.0`。

### D. 参考文献

[1] Qu X, Xu P, Jiang H, He L, Ren X. DFT+U within the framework of linear combination of numerical atomic orbitals. *The Journal of Chemical Physics*. 2022;156(23):234104. [doi:10.1063/5.0090122](https://doi.org/10.1063/5.0090122)
