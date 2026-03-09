---
title: "ABACUS 随机密度泛函理论（SDFT）使用教程"
---

# ABACUS 随机密度泛函理论（SDFT）使用教程

## 教程目标

本教程介绍 ABACUS 中随机波函数密度泛函理论（SDFT）和混合随机密度泛函理论（MDFT）的使用方法，以三个官方算例为主线，覆盖电子自洽迭代（SCF）、分子动力学（MD）和态密度（DOS）三种计算场景。

## 适用读者

- 需要模拟高温高压体系（温稠密物质，WDM）的研究者
- 已有 ABACUS KSDFT 使用经验、需要迁移到 SDFT 的用户

## 前置知识

- ABACUS 基本输入文件（INPUT / STRU / KPT）的准备
- DFT 基本概念（赝势、平面波基组、自洽迭代）

## 教程结构

| 章节 | 内容 |
|------|------|
| 第一章 | SDFT/MDFT 背景与适用场景 |
| 第二章 | 算例下载与文件准备 |
| 第三章 | 案例一：SCF 计算（Si，MDFT，T=8.16 eV） |
| 第四章 | 案例二：分子动力学（Al，纯SDFT，T=100 eV） |
| 第五章 | 案例三：态密度计算（Si，MDFT+DOS，T=8.16 eV） |
| 第六章 | 精度控制与参数选取 |
| 附录 | 参考文献、常见问题 |

## 一、引言

### 1.1 KSDFT 在高温下的局限

Kohn-Sham 密度泛函理论（KSDFT）通过对角化哈密顿矩阵求解电子波函数，计算量正比于电子数的三次方（O(N³)）。在常温下，只有费米能级附近的少量轨道被占据，计算尚可承受。

当体系温度升高到数十至上千 eV 时（1 eV ≈ 11604.5 K），费米-狄拉克分布的展宽随温度增大，大量高能轨道获得可观的占据数，需要纳入计算的波函数数目急剧增加。温稠密物质（Warm Dense Matter，WDM）广泛存在于行星内部、惯性约束聚变等场景，使用传统 KSDFT 计算这类体系几乎无法实现。

### 1.2 SDFT 的核心思路

随机波函数密度泛函理论（Stochastic DFT，SDFT）用一组随机产生的波函数轨道代替对角化，通过随机轨道对体系的迹（Trace）进行统计估计，从而完全绕开哈密顿矩阵对角化这一计算瓶颈。

SDFT 有两个显著特点：

1. **计算量与温度成反比**：温度越高，切比雪夫展开所需阶数越少，计算反而更快。这与 KSDFT 在高温时变慢的趋势完全相反。
2. **并行效率高**：不同随机轨道的操作互相独立，可实现近线性的并行扩展。

ABACUS 中 SDFT 的算法细节见：Qianrui Liu and Mohan Chen, *Phys. Rev. B* **106**, 125132 (2022)。

### 1.3 MDFT：混合加速策略

纯 SDFT 的随机误差随轨道数增加缓慢收敛。混合随机密度泛函理论（Mixed stochastic-deterministic DFT，MDFT）引入少量低能 Kohn-Sham 轨道，将随机误差的来源限制在高能部分，收敛速度大幅加快。

在 ABACUS 中，通过 `nbands` 和 `nbands_sto` 两个参数控制两种轨道的数目：
- `nbands = 0`，`nbands_sto > 0`：纯 SDFT
- `nbands > 0`，`nbands_sto > 0`：MDFT

### 1.4 本教程的三个算例

| 算例 | 体系 | 温度 | 模式 | 计算类型 |
|------|------|------|------|----------|
| pw_Si2 | 2原子Si（金刚石）| 0.6 Ry（≈8.16 eV）| MDFT | SCF |
| pw_md_Al | 16原子Al | 7.35 Ry（≈100 eV）| 纯SDFT | MD |
| 186_PW_SDOS_10D10S | 1原子Si | 0.6 Ry（≈8.16 eV）| MDFT | SCF+DOS |

三个算例从简到繁，分别展示最基础的 SCF、高温 MD 和态密度的完整使用流程。

## 二、准备工作

### 2.1 软件版本

SDFT 功能从 ABACUS 2.3.0 版本开始支持，DOS 计算和多 K 点采样等完整功能从 3.2.0 版本起可用。本教程基于 ABACUS 3.2.0 版本的算例。

### 2.2 算例下载

官方算例包可从 Gitee 获取：

```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
```

或从 GitHub 获取：

```bash
git clone https://github.com/MCresearch/abacus-user-guide.git
```

下载后进入随机波函数算例目录：

```bash
cd abacus-user-guide/examples/stochastic
```

目录中包含三个文件夹：

```
stochastic/
├── pw_Si2/                  # 案例一：Si SCF
├── pw_md_Al/                # 案例二：Al MD
└── 186_PW_SDOS_10D10S/      # 案例三：Si DOS
```

每个文件夹均包含 `INPUT`、`STRU`、`KPT` 三个输入文件，与常规 KSDFT 计算的文件结构完全相同。赝势文件位于 `../../PP_ORB/` 目录（相对于各算例文件夹）。

### 2.3 赝势说明

案例一使用 `Si.pz-vbc.UPF`，包含硅的 4 个价电子。

当电子温度特别高（如超过 100 eV）时，内壳层电子可能被热激发，常规赝势的可移植性会下降。此时需要选用包含更多内壳层电子的赝势，甚至需要定制生成新的赝势。ABACUS 目前支持模守恒赝势（NCPP）。

## 三、案例一：SCF 计算（Si，MDFT，T ≈ 8.16 eV）

### 3.1 场景说明

`pw_Si2` 是一个 2 个原子的金刚石结构硅体系，电子温度设为 0.6 Ry（约 8.16 eV）。本例采用 MDFT 模式：4 条 KS 轨道（覆盖费米面以下的低能部分）加 64 条随机轨道，执行电子自洽迭代（SCF）计算。

这是 SDFT 功能最基础的入口案例，所涉及的核心参数在此一并详解。

### 3.2 INPUT 文件

```
INPUT_PARAMETERS
#Parameters (General)
calculation    scf
esolver_type   sdft
pseudo_dir     ../../PP_ORB
nbands         4
nbands_sto     64
nche_sto       100
method_sto     1
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  0.6
```

### 3.3 关键参数详解

**`esolver_type`**
选择系统总能量的求解方式。默认值为 `ksdft`，设为 `sdft` 后程序将根据 `nbands` 和 `nbands_sto` 的取值自动判断使用纯 SDFT 还是 MDFT。

**`nbands`**
确定性 KS 轨道数（deterministic orbitals），通过严格对角化矩阵计算得到。一般取费米面以下所有低能轨道的数目，效率最高。本例中取 4，对应 Si 的 4 个价电子轨道。

- `nbands = 0` 且 `nbands_sto > 0`：纯 SDFT
- `nbands > 0` 且 `nbands_sto > 0`：MDFT（本例）

**`nbands_sto`**
随机波函数轨道数（stochastic orbitals）。取值越大，随机误差越小，但计算量相应增加。本例取 64。

如何判断 64 是否足够？见第六章的收敛测试方法。

**`nche_sto`**
切比雪夫展开阶数。对哈密顿量进行切比雪夫多项式展开以评估费米-狄拉克函数，阶数越高精度越高、效率越低。

经验规律：
- 与温度**成反比**：温度越高，所需阶数越少
- 与 `ecutwfc`（平面波截断能）**成正比**：截断能越大，所需阶数越多

判断标准：查看输出文件 `running_scf.log`，找到 `Chebyshev Precision` 一行，确保其值小于 `1e-8`。本例 T=8.16 eV，取 100 阶。

**`method_sto`**
SDFT 内部计算方法的选择：
- `1`：内存消耗较小，速度稍慢
- `2`：速度更快，但需要更多内存（默认值）

本例取 `1`，适合内存有限的情形。

**`smearing_method`**
展宽方式。SDFT/MDFT 当前**只支持** `fd`（Fermi-Dirac）展宽，不能使用高斯或 Methfessel-Paxton 等方式。

**`smearing_sigma`**
电子温度，单位为 **Ry**（里德伯，1 Ry ≈ 13.6 eV）。本例 0.6 Ry ≈ 8.16 eV。注意这是电子温度，在 MD 计算中与离子温度（`md_tfirst`）是两个独立参数。

**常规参数**

| 参数 | 值 | 说明 |
|------|----|------|
| `ecutwfc` | 50 Ry | 平面波截断能 |
| `scf_nmax` | 20 | 最大 SCF 迭代步数 |
| `symmetry` | 1 | 开启晶体对称性 |

### 3.4 STRU 和 KPT

STRU 文件描述金刚石结构 Si 的晶胞和原子位置，格式与常规 KSDFT 计算完全相同。KPT 文件设置布里渊区 K 点采样，ABACUS 的 SDFT 支持多 K 点，在某些性质（如 DOS）计算中要注意 K 点收敛性。

### 3.5 运行方法

```bash
cd pw_Si2
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

计算完成后，结果位于 `OUT.ABACUS/` 目录，总能量和收敛信息可在 `running_scf.log` 中查看。

### 3.6 注意事项

1. **K 点收敛**：SDFT 支持多 K 点采样，高精度计算前应测试 K 点数目对结果的影响。

3. **`nbands_sto = 0` 的特殊行为**：ABACUS 3.2.2 版本以后，若将 `nbands_sto` 设为 0，程序会自动切换为 KSDFT 模式计算。

## 四、案例二：分子动力学模拟（Al，纯SDFT，T ≈ 100 eV）

### 4.1 场景说明

`pw_md_Al` 是一个 16 个铝原子的体系，电子温度设为 7.35 Ry（约 100 eV）。本例使用**纯 SDFT**（`nbands = 0`，没有任何 KS 轨道），执行分子动力学（MD）模拟。

极端高温（100 eV）是 SDFT 最能发挥优势的场景：此时 KSDFT 需要数百条以上的波函数，而 SDFT 只需 64 条随机轨道，且切比雪夫展开仅需 20 阶（远少于低温时的 100 阶），计算效率大幅提升。

### 4.2 INPUT 文件

```
INPUT_PARAMETERS
#Parameters (General)
calculation    md
esolver_type   sdft
pseudo_dir     ../../PP_ORB
nbands         0
nbands_sto     64
nche_sto       20
method_sto     2
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
scf_thr        1e-6
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  7.34986072
#Parameters (MD)
md_tfirst      1160400
md_dt          0.2
md_nstep       10
```

### 4.3 关键参数详解

**`calculation = md`**
设为 `md` 启动分子动力学模拟。在每个 MD 步中，程序先执行 SDFT 的 SCF 自洽迭代计算出受力，再根据受力更新原子位置。

**`nbands = 0`**
本例不使用任何 KS 轨道，全部依赖随机轨道，是纯 SDFT 模式。与案例一（`nbands = 4`，MDFT）相比，适用于极高温场景——此时费米面附近已无明显的低能 KS 轨道可供利用。

**`nche_sto = 20`**
温度 100 eV 时，费米-狄拉克函数变得非常平滑，切比雪夫展开只需 20 阶即可达到足够精度。相比案例一（T=8.16 eV，100 阶），高温大幅降低了展开阶数，这也是 SDFT 在高温下效率优异的根本原因。

**`method_sto = 2`**
高温 MD 计算通常需要更快的速度，选用方法 2（速度更快，内存更多）。

**`smearing_sigma = 7.34986072`**
电子温度，单位 Ry，约等于 100 eV（7.35 Ry × 13.6 eV/Ry ≈ 99.96 eV）。

**MD 参数**

| 参数 | 值 | 说明 |
|------|----|------|
| `md_tfirst` | 1160400 K | 离子初始温度（1160400 K ÷ 11604.5 K/eV ≈ 100 eV） |
| `md_dt` | 0.2 fs | MD 时间步长 |
| `md_nstep` | 10 | MD 模拟总步数 |

注意 `md_tfirst` 是**离子温度**（单位 K），而 `smearing_sigma` 是**电子温度**（单位 Ry）。高温 WDM 模拟中，两者通常取相同的等效温度。

### 4.4 纯 SDFT 还是 MDFT？

| 场景 | 推荐模式 | 原因 |
|------|----------|------|
| 温度 < 10 eV | MDFT（`nbands > 0`） | 低能 KS 轨道可显著减少随机误差 |
| 温度 ≥ 10 eV | 纯 SDFT（`nbands = 0`）或 MDFT | SDFT 效率已足够，MDFT 收益递减 |
| 极端高温（> 50 eV）| 纯 SDFT | KS 轨道对角化成本高，不值得 |

### 4.5 运行方法

```bash
cd pw_md_Al
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

MD 轨迹和各步能量、受力、压强等信息输出在 `OUT.ABACUS/` 目录的 `running_md.log` 和 `MD_dump` 文件中。

### 4.6 注意事项

1. **`scf_thr = 1e-6`**：MD 中每一步的 SCF 不需要收敛到静态计算的精度，`1e-6` 已足够。
2. **高温赝势**：100 eV 下赝势可移植性需重点评估，参见第二章 2.3。
3. **SDFT 随机误差在 MD 中的影响**：每步 SCF 的受力含有随机误差，会引入额外的随机"热噪声"。通过足够多的 `nbands_sto` 可将误差控制在可接受范围。

## 五、案例三：态密度计算（Si，MDFT+DOS，T ≈ 8.16 eV）

### 5.1 场景说明

`186_PW_SDOS_10D10S` 是一个 1 个 Si 原子的体系，电子温度 0.6 Ry（约 8.16 eV），采用 MDFT（10 条 KS 轨道 + 10 条随机轨道）计算 SCF 并输出态密度（DOS）。算例名中的 `10D10S` 正是指 10 条确定性轨道和 10 条随机轨道。

与前两个算例相比，本例的核心增量是：在 SDFT/MDFT 的 SCF 基础上，通过一组专用参数控制态密度的能量范围、展宽和精度，最终在 `OUT.autotest/` 目录下生成 `DOS1_smearing.dat` 文件。

### 5.2 INPUT 文件

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix          autotest
calculation     scf
esolver_type    sdft
method_sto      2
nbands          10
nbands_sto      10
nche_sto        120
emax_sto        0
emin_sto        0
seed_sto        20000
pseudo_dir      ../../PP_ORB
symmetry        1
kpar            1
bndpar          2
#Parameters (2.Iteration)
ecutwfc         20
scf_thr         1e-6
scf_nmax        20
#Parameters (3.Basis)
basis_type      pw
#Parameters (4.Smearing)
smearing_method fd
smearing_sigma  0.6
#Parameters (5.Mixing)
mixing_type     broyden
mixing_beta     0.4
out_dos         1
dos_emin_ev     -20
dos_emax_ev     100
dos_edelta_ev   0.1
dos_sigma       4
dos_nche        240
npart_sto       2
```

### 5.3 SDFT/MDFT 控制参数（本例特有项）

**`nbands = 10` / `nbands_sto = 10`**
MDFT 模式：10 条低能 KS 轨道 + 10 条随机轨道。算例名称中的 `10D10S` 即由此得名（10 Deterministic + 10 Stochastic）。

**`nche_sto = 120`**
同样 T=8.16 eV，本例用 120 阶（多于案例一的 100 阶）。DOS 计算对精度要求更高，建议适当增加展开阶数，确认 `Chebyshev Precision < 1e-8`。

**`emax_sto = 0` / `emin_sto = 0`**
随机轨道的能量范围上下限。设为 0 时程序自动确定，通常无需手动设置。

**`seed_sto = 20000`**
固定随机种子，确保计算结果可复现。收敛测试时应使用不同的 `seed_sto` 值（如 1、2、3……），以评估随机误差。

**`bndpar = 2`**
将所有 MPI 进程分成 2 组，随机轨道平均分配到每组，提高并行效率。实际使用时优先用 `kpar` 进行 K 点并行，再测试不同的 `bndpar` 取值。`bndpar` 并不是越大越好。

### 5.4 DOS 专用参数

| 参数 | 值 | 说明 |
|------|----|------|
| `out_dos` | 1 | 开关：设为 1 才会输出 DOS |
| `dos_emin_ev` | -20 eV | DOS 能量范围下限 |
| `dos_emax_ev` | 100 eV | DOS 能量范围上限 |
| `dos_edelta_ev` | 0.1 eV | DOS 能量间隔（分辨率） |
| `dos_sigma` | 4 eV | 高斯展宽宽度（WDM 体系用较大展宽） |
| `dos_nche` | 240 | DOS 后处理切比雪夫展开阶数（独立于 `nche_sto`） |
| `npart_sto` | 2 | 内存控制：使用正常内存的 1/2 |

**`dos_nche`** 与 `nche_sto` 是两个独立参数：前者专用于 DOS 的后处理步骤，通常需要比 `nche_sto` 取更大的值以获得高质量的态密度曲线。本例取 240，是 `nche_sto`（120）的两倍。

**`npart_sto`** 用于处理 `method_sto = 2` 时 DOS 后处理内存不足的问题。设为 `n` 时，程序将内存使用量控制在正常的 1/n，以较慢的速度换取内存的节省。默认值为 1（不拆分），建议在内存受限时设为 2 或更大。

### 5.5 其他参数

**`mixing_type = broyden` / `mixing_beta = 0.4`**
Broyden 混合方法，`mixing_beta` 控制每步更新的比例。SDFT 的 SCF 收敛有时比 KSDFT 稍慢，可适当降低 `mixing_beta`（如从 0.7 降至 0.4）来改善收敛。

**`ecutwfc = 20`**
本例截断能取 20 Ry，低于案例一（50 Ry）。

### 5.6 运行与输出

```bash
cd 186_PW_SDOS_10D10S
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

输出文件位于 `OUT.autotest/`：

- `running_scf.log`：SCF 收敛信息，含 `Chebyshev Precision` 检查项
- `DOS1_smearing.dat`：态密度数据文件，两列格式（能量(eV) / DOS强度）

### 5.7 注意事项

1. **DOS 能量范围选取**：`dos_emin_ev` 和 `dos_emax_ev` 应覆盖感兴趣的能量区间。对高温 WDM，费米面以上数十 eV 仍有可观态密度，上限取 100 eV 是合理的。

2. **`dos_sigma` 的选取**：展宽过小会导致 DOS 曲线有数值噪声（源于随机误差），展宽过大会掩盖能量细节。WDM 体系通常用 1~5 eV 的展宽。

3. **`kpar` 与 `bndpar` 的优先级**：`kpar`（K 点并行）的并行效率通常高于 `bndpar`。有多个 K 点时，优先使用 `kpar`，再调整 `bndpar`。

## 六、精度控制与参数选取

SDFT 引入了随机误差，参数选取策略与 KSDFT 有所不同。本章给出三个关键参数的收敛测试方法和选取建议。

### 6.1 `nbands_sto` 收敛测试

随机轨道数决定随机误差的大小。确定合理的 `nbands_sto` 步骤如下：

1. 选定初始值（如 `nbands_sto = 32`）
2. 准备 10 组不同随机种子：在 INPUT 中修改 `seed_sto`（如分别设为 1、2、3……10），每组独立运行一次 SDFT 计算
3. 取 10 组总能量的平均值和标准差
4. 若标准差 / |平均值| > 1e-4（万分之一），则增大 `nbands_sto`，重复上述步骤
5. 直到能量相对误差 < 1e-4 为止

**经验参考值：**

| 温度 | 体系大小 | 建议 `nbands_sto` |
|------|----------|-------------------|
| < 10 eV（MDFT）| 小体系（~10 原子）| 10~64 |
| < 10 eV（SDFT）| 小体系 | 64~256 |
| > 50 eV | 任意 | 32~64（高温误差本身较小） |

### 6.2 `nche_sto` 选取

切比雪夫展开阶数直接影响费米-狄拉克函数的近似精度。

**判断标准：** 运行计算后，在 `running_scf.log` 中搜索 `Chebyshev Precision`，该值应小于 `1e-8`。若超过此阈值，增大 `nche_sto` 后重新计算。

**经验规律：**

$$\text{nche\_sto} \propto \frac{\text{ecutwfc}}{T}$$

| 温度 | `ecutwfc = 50 Ry` 时参考阶数 |
|------|------------------------------|
| 100 eV（7.35 Ry） | ~20 |
| 8.16 eV（0.6 Ry） | ~100 |
| 1 eV（0.073 Ry） | ~800 |

温度越低，所需阶数越高，计算越慢——这正是 SDFT 适合高温计算的内在原因。

### 6.3 `ecutwfc` 收敛测试

SDFT 的 `ecutwfc` 测试原理与 KSDFT 相同，但因存在随机误差，需要统计处理：

1. 固定好 `nbands_sto`
2. 对每个 `ecutwfc` 值（从基准值起每次增加 10 Ry），使用 10 个不同随机种子计算，取平均能量
3. 相邻两个 `ecutwfc` 对应的平均能量之差 < 目标精度（如万分之一）时，认为收敛

WDM 体系的总能量一般很高（绝对值大），建议使用**相对误差**（而非绝对值）作为收敛判断标准。

### 6.4 并行策略

| 并行方式 | 参数 | 效率 | 建议 |
|----------|------|------|------|
| K 点并行 | `kpar` | 高 | 优先使用，设为 K 点总数的因子 |
| 随机轨道并行 | `bndpar` | 中 | 在 `kpar` 充分后再调整 |

`bndpar` 并非越大越好，建议在目标机器上实测 `bndpar = 1, 2, 4` 的效率，选最优值。

### 6.5 参数速查表

| 参数 | 影响 | 增大效果 | 减小效果 |
|------|------|----------|----------|
| `nbands_sto` | 随机误差 | 误差降低，效率降低 | 误差增大，效率提高 |
| `nbands`（MDFT）| 系统误差 | 误差降低（饱和后不变），效率降低 | — |
| `nche_sto` | 展开精度 | 精度提高，效率降低 | 精度降低（慎用） |
| `ecutwfc` | 基组完备性 | 精度提高，效率降低 | 精度降低（慎用） |
| `method_sto` | 速度/内存 | 2比1更快但更耗内存 | — |

## 附录

### 参考文献

1. Qianrui Liu and Mohan Chen, "Plane-wave-based stochastic-deterministic density functional theory for extended systems," *Phys. Rev. B* **106**, 125132 (2022). DOI: 10.1103/PhysRevB.106.125132

2. ABACUS 官方文档 — SDFT 参数说明：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#electronic-structure-sdft

3. 算例下载（Gitee）：https://gitee.com/mcresearch/abacus-user-guide/tree/master/examples/stochastic

### 常见问题

**Q：Chebyshev Precision 始终大于 1e-8，怎么办？**
A：增大 `nche_sto`，每次增加 20~50，直到 Chebyshev Precision 满足要求。若 `ecutwfc` 较大，所需阶数会显著增加。

**Q：SCF 不收敛怎么办？**
A：尝试降低 `mixing_beta`（如从 0.7 改为 0.4），或换用 `mixing_type = broyden`。SDFT 的 SCF 收敛比 KSDFT 稍慢属于正常现象。

**Q：DOS 计算内存溢出？**
A：设置 `npart_sto = 2` 或更大，程序将分批处理，内存用量减半（但速度相应降低）。

**Q：`nbands_sto` 增大到多少都有随机噪声？**
A：高温（>50 eV）时体系的电子分布非常弥散，仅靠随机轨道也难以彻底消除噪声。此时可同时增大 `nbands_sto` 和使用 K 点并行（多 K 点对不同 K 处的随机轨道取平均，额外降低误差）。

### 进阶学习

- **电子电导与热导**：ABACUS 支持在 SDFT 框架下计算电子电导率和热导率（`cal_cond = 1`），详见 ABACUS 官方教程"SDFT 计算电子电导热导"。
- **OFDFT**：对于极高温（>1000 eV）或超大体系，无轨道密度泛函理论（OF-DFT）是另一种选择，ABACUS 也提供相应功能。
- **自定义赝势**：高温下需要自制赝势时，可参考 ABACUS 文档中模守恒赝势的生成指南。
