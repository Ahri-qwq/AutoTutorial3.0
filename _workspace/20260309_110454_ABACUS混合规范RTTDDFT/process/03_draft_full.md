# 前言

## 教程目标

本教程介绍如何在 ABACUS 中使用混合规范（hybrid gauge）进行实时含时密度泛函理论（RT-TDDFT）计算。重点解决以下问题：

- 为什么在周期性体系中用速度规范会出错？
- 混合规范如何修正这个问题？只需改动哪个参数？
- 如何用混合规范计算 Si 的介电函数和高次谐波谱（HHG）？

## 适用读者

- 熟悉 ABACUS 基本操作（能跑 SCF），想进一步做 RT-TDDFT 计算的用户
- 对线性光学响应或非线性光学（HHG）感兴趣的研究者
- 已有 RT-TDDFT 基础，想了解混合规范优势的用户

## 前置知识

- ABACUS 基本使用（INPUT/STRU/KPT 文件格式）
- 密度泛函理论基本概念
- 了解 RT-TDDFT 原理（可参考 ABACUS RT-TDDFT 基础教程）

## 教程结构

1. 引言：RT-TDDFT 的用途与周期体系的挑战
2. 从 SCF 到 TDDFT：输入文件准备，以及速度规范为何在 NAO 下不够准确
3. 混合规范：一个参数解决问题，以及物理思想简介
4. 案例一：HCO 分子（非周期体系验证，三种规范对比）
5. 案例二：Si 单胞——介电函数与高次谐波谱
6. 常见问题与注意事项
7. 附录

## 关于本文

混合规范 RT-TDDFT 方法由中国科学技术大学赵昊天与何力新教授提出，已在 ABACUS 中实现，相关成果发表于 *Journal of Chemical Theory and Computation*（2025, 21, 3335–3341）。本文基于该工作的理论和测试案例编写，所有计算设置均以 ABACUS LCAO 基组为前提。
# 一、引言

RT-TDDFT 是研究材料在激发态和外场作用下电子动力学的第一性原理方法。与线性响应 TDDFT 不同，RT-TDDFT 直接在时域中传播含时 Kohn-Sham 方程，能处理强场、非线性、超快等线性响应方法无法触及的问题。常见应用包括：

- **光吸收谱**：通过弱脉冲激发，提取材料的线性光学响应
- **介电函数**：周期性体系的光学性质
- **高次谐波谱（HHG）**：强场下的非线性响应，谐波阶次可达 7 次乃至更高
- **超快载流子动力学**：光激发后的电子弛豫过程

ABACUS 的 RT-TDDFT 基于数值原子轨道（NAO）基组实现，支持 LCAO 框架下的实时传播。这一选择带来了效率优势，但也引出了一个问题：**如何在 NAO 基组下正确地引入外加电场？**

这个问题的答案就是规范选择。传统方法中，速度规范在理论上适用于周期性体系，但在 NAO 基组下会引入系统性误差，导致计算结果不可靠——尤其是在小基组和强场条件下。混合规范通过引入矢势诱导的局域相位修正，有效解决了这一问题，并同时保持了较高的计算效率。

本教程将带你从 SCF 出发，逐步搭建周期性体系的 RT-TDDFT 计算，重点展示如何用 `td_stype = 2` 切换到混合规范，以及这一改动对结果的影响。
# 二、从基态 SCF 到 TDDFT 演化

RT-TDDFT 计算的起点是基态 SCF。先跑 SCF，再在 SCF 基态上开始时间演化。本章以 Si 原胞为例，介绍完整的输入文件准备，并讨论规范选择的背景。

## 2.1 必要输入文件

ABACUS 的 RT-TDDFT 计算需要四类输入文件：

| 文件 | 作用 |
|------|------|
| `INPUT` | 控制计算类型、基组、时间步、电场参数等 |
| `STRU` | 晶体结构（晶格矢量 + 原子位置） |
| `KPT` | 布里渊区 k 点采样 |
| 赝势 (`.upf`) 和轨道 (`.orb`) 文件 | 元素的赝势和数值原子轨道 |

### Si 原胞结构（STRU）

Si 金刚石结构原胞含 2 个原子，采用面心立方原胞表示：

```
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
   0.5000000000     0.5000000000     0.0000000000
   0.0000000000     0.5000000000     0.5000000000
   0.5000000000     0.0000000000     0.5000000000

ATOMIC_POSITIONS
Direct
Si
0
2
   0.0000000000     0.0000000000     0.0000000000  m  1  1  1
   0.2500000000     0.2500000000     0.2500000000  m  1  1  1
```

晶格常数 10.2 Bohr ≈ 5.40 Å，接近 Si 的实验值（5.43 Å）。LCAO 基组选用 `Si_gga_8au_60Ry_2s2p1d.orb`（截断半径 8 Bohr，60 Ry 截断能，2s2p1d 构型）。

### k 点设置（KPT）

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

RT-TDDFT 对 k 点采样的要求与静态计算类似，周期体系中需合理采样布里渊区。4×4×4 的 Gamma 中心网格适用于 Si 原胞的基础计算。

### 基态 SCF 的 INPUT

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_scf
calculation      scf
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Accuracy)
ecutwfc          60
scf_thr          1e-8
scf_nmax         100

#Parameters (3.Basis)
basis_type       lcao
gamma_only       0

#Parameters (4.Smearing)
smearing_method  gauss
smearing_sigma   0.01
```

> **注意：** `gamma_only` 必须设为 `0`。RT-TDDFT 使用传播子迭代计算复数波函数，不兼容 `gamma_only = 1`（后者只存储实数）。这个设置在 SCF 阶段就必须确定，以保证 TDDFT 阶段能读取正确格式的电荷密度。

## 2.2 为什么周期体系必须用速度规范——以及它的局限性

在 RT-TDDFT 中，加入外加电场的方式取决于**规范选择**。从电动力学的基本关系出发：

$$\mathbf{E} = -\nabla V - \frac{\partial \mathbf{A}}{\partial t}$$

有两种常用选择：

**长度规范**：将电场归为纯标势，加入 $V \propto \mathbf{E} \cdot \mathbf{r}$。哈密顿量形式简单，计算快，但 $\mathbf{r}$ 在周期性体系边界处不连续——只能用锯齿状势来衔接，对周期性体系不适用。

**速度规范**：将电场归为矢势 $\mathbf{A}(t) \propto \int_0^t \mathbf{E}(t') dt'$，保持哈密顿量的周期性。这是处理周期性体系的标准方法。

**速度规范在 NAO 下的问题**

理论上，长度规范和速度规范通过如下规范变换完全等价：

$$\psi^{\text{velocity}}(\mathbf{r}, t) = e^{i\mathbf{A}(t)\cdot\mathbf{r}} \psi^{\text{length}}(\mathbf{r}, t)$$

这个相位因子 $e^{i\mathbf{A}(t)\cdot\mathbf{r}}$ 描述了矢势引起的空间相位变化。在完备基组下，两种规范的结果完全一致。

但在 NAO 基组下，波函数被展开为原子轨道的线性组合：

$$\psi(\mathbf{r}, t) = \sum_{\mu} c_\mu(t) \phi_\mu(\mathbf{r})$$

这里的原子轨道 $\phi_\mu(\mathbf{r})$ 在规范变换时不会自动获得 $e^{i\mathbf{A}(t)\cdot\mathbf{r}}$ 相位。结果是：直接使用速度规范的 NAO 展开，**忽略了矢势在每个原子轨道内部引起的相位变化**。

这一遗漏导致计算的电流和能量出现系统性误差。在响应电流的计算中，低频区域会出现明显的虚假发散；在介电函数等关键物理量上，精度随基组减小而快速降低。

增大基组可以缓解误差，但计算成本急剧上升——这正是 NAO 基组的速度优势被抵消的地方。**混合规范直接从根源上修正这一误差。**
# 三、混合规范：一个参数解决问题

## 3.1 物理思想

混合规范的核心想法很直接：既然速度规范的问题在于原子轨道缺少矢势诱导的相位，那就显式地把这个相位加回来。

具体做法是在每个原子轨道前引入一个依赖于矢量势的相位因子：

$$\tilde{\phi}_\mu(\mathbf{r}, t) = e^{i\mathbf{A}(t)\cdot\mathbf{R}_\mu} \phi_\mu(\mathbf{r})$$

其中 $\mathbf{R}_\mu$ 是第 $\mu$ 个原子轨道的中心位置，$\mathbf{A}(t)$ 是时变矢势。波函数展开变为：

$$\psi(\mathbf{r}, t) = \sum_{\mu} c_\mu(t) \, e^{i\mathbf{A}(t)\cdot\mathbf{R}_\mu} \phi_\mu(\mathbf{r})$$

在这个修正下，哈密顿量矩阵元和重叠矩阵均需相应修改。由于混合规范同时包含矢量势 $\mathbf{A}$ 和标量场，它本质上是一种新的规范选择——相比长度规范，混合规范下的外场项保持周期性，不破坏布洛赫定理；相比速度规范，它补偿了 NAO 基组不完备带来的相位误差。

混合规范和长度规范之间满足严格的规范变换关系，理论上结果完全等价。在 ABACUS 已有测试中，混合规范与长度规范（非周期体系）的结果高度吻合，而速度规范则存在明显偏差。

> 深入的理论推导参见原始论文：赵昊天, 何力新, *J. Chem. Theory Comput.* **2025**, 21, 3335–3341 (DOI: 10.1021/acs.jctc.5c00111)

## 3.2 在 ABACUS 中启用混合规范

启用混合规范只需设置一个参数：

```
td_stype   2     # 0=长度规范  1=速度规范  2=混合规范
```

以下是一个完整的 TDDFT INPUT 示例（Si 周期体系，弱场线性响应）：

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_tddft
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         4000          # 总时间步数 = 4000 × 0.02 fs = 80 fs
md_dt            0.02          # 时间步长 (fs)，建议 0.02–0.05
md_tfirst        0             # 离子初始温度，0 表示固定离子

#Parameters (4.外场设置)
td_vext          1             # 开启外加电场
td_stype         2             # 混合规范（重点参数）

td_tstart        1             # 第 1 步开始施加电场
td_tend          4000          # 与 md_nstep 相同（全程施加）

td_vext_dire     1             # 沿 x 方向施加电场

td_ttype         0             # 高斯型脉冲
td_gauss_amp     0.001         # 振幅 (V/Å)，弱场取 ~0.001
td_gauss_t0      50            # 脉冲中心位于第 50 步
td_gauss_sigma   0.5           # 高斯宽度 (fs)，宽带覆盖
td_gauss_freq    0.0           # 中心频率 0 (宽带激发)
td_gauss_phase   0.0           # 初相位

#Parameters (5.输出)
out_current      1             # 输出电流（周期体系/混合规范用此项）
out_efield       1             # 输出外加电场时程
```

> **与速度规范的区别：** 将 `td_stype 1` 改为 `td_stype 2`，其余所有参数保持不变。

## 3.3 外场参数速查

| 参数 | 说明 | 典型值 |
|------|------|--------|
| `td_vext` | 外场开关 | `1`（开） |
| `td_stype` | 规范选择 | `0`/`1`/`2`（长度/速度/混合） |
| `td_tstart` | 开始步数 | `1` |
| `td_tend` | 结束步数 | 通常等于 `md_nstep` |
| `td_vext_dire` | 方向 | `1`=x, `2`=y, `3`=z |
| `td_ttype` | 电场形状 | `0`=高斯, `1`=梯形, `2`=三角, `3`=阶跃 |
| `td_gauss_amp` | 振幅 (V/Å) | 线性响应 ~0.001；强场 0.01–0.1 |
| `td_gauss_t0` | 中心步数 | 脉冲中心位置 |
| `td_gauss_sigma` | 宽度 (fs) | 0.2–1.0；越小带宽越大 |
| `td_gauss_freq` | 中心频率 | `0` 为宽带；特定值则为窄带 |
| `out_current` | 电流输出 | `1`（周期体系用） |
| `out_dipole` | 偶极矩输出 | `1`（非周期体系用） |
| `out_efield` | 外场输出 | `1` |

**规范选择建议**

| 体系类型 | 推荐规范 | `td_stype` |
|---------|---------|------------|
| 非周期（分子、团簇） | 长度规范或混合规范 | `0` 或 `2` |
| 周期（晶体） | **混合规范**（推荐） | `2` |
| 周期（旧方法兼容） | 速度规范 | `1` |

在 ABACUS LCAO 框架下，**混合规范适用于所有体系**，精度优于速度规范，与长度规范等价（非周期）。对周期体系，混合规范是目前最可靠的选择。
# 四、案例验证：HCO 分子（非周期体系）

本章用 HCO 分子（甲酰基，非周期体系）验证混合规范的正确性，并直观展示速度规范的系统性误差。

## 4.1 体系设置

HCO 分子放置于含真空层的超胞中，沿所有方向留出足够真空（≥10 Å），确保周期性边界条件不影响分子响应。

**INPUT（混合规范，HCO 分子）**

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           HCO_hybrid
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         2000          # 40 fs 演化
md_dt            0.02
md_tfirst        0

#Parameters (4.外场设置)
td_vext          1
td_stype         2             # 混合规范
td_tstart        1
td_tend          2000

td_vext_dire     1             # x 方向
td_ttype         0
td_gauss_amp     0.01
td_gauss_t0      300
td_gauss_sigma   0.2
td_gauss_freq    3.66          # eV/hbar
td_gauss_phase   0.0

#Parameters (5.输出)
out_dipole       1             # 非周期体系：输出电偶极矩
out_efield       1
```

> **注意**：对非周期体系，应输出电偶极矩（`out_dipole 1`），而非电流。分子体系中也可使用长度规范（`td_stype 0`），此时 `out_dipole` 同样适用。

非周期体系的 KPT 文件只需 Gamma 点：
```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

## 4.2 结果：三种规范的对比

施加高斯脉冲激发后，从输出的 `SPIN1_DIPOLE` 文件读取各时刻的电偶极矩 $D(t)$，配合 `efield_0.dat` 中的外加电场 $E(t)$，通过傅里叶变换得到吸收谱：

$$\alpha_i(\omega) = \frac{\int D_i(t) e^{i\omega t} dt}{\int E_i(t) e^{i\omega t} dt}$$

其中 $i$ 为激发方向，$\alpha_i(\omega)$ 的虚部对应吸收强度。

后处理使用 ABACUS 自带的 `tools/plot-tools` 脚本：

```bash
python tools/plot-tools/plot_absorption.py \
    --dipole SPIN1_DIPOLE \
    --efield efield_0.dat \
    --output absorption_HCO.dat
```

**三种规范的结果对比（HCO 分子）**

| 规范 | `td_stype` | 吸收谱精度 | 低频行为 |
|------|-----------|-----------|---------|
| 长度规范 | `0` | 参考基准 | 正常 |
| 混合规范 | `2` | 与长度规范完全一致 | 正常 |
| 速度规范 | `1` | 明显偏差 | 低频发散 |

赵昊天等的计算结果显示，对 HCO 分子，长度规范与混合规范的响应电流和吸收谱结果高度吻合，验证了混合规范在非周期体系中的正确性。速度规范的响应电流与前两者存在明显偏差，对应的吸收谱在低频区出现显著的虚假发散，这正是 NAO 基组下忽略轨道内部相位变化导致的系统性误差的体现。

这个对比实验说明：**混合规范可以完全替代长度规范用于非周期体系，同时为周期体系提供准确的结果**，而速度规范在 NAO 基组下无论对哪类体系都存在不可忽视的误差。

> 参考算例：https://gitee.com/mcresearch/abacus-user-guide/tree/master/examples/tddft
# 五、核心案例：Si 单胞（周期体系）

Si 是验证周期体系 RT-TDDFT 方法的标准体系。本章包含两个计算任务：
1. 弱场激发 → 线性响应 → 介电函数
2. 强场激发 → 非线性响应 → 高次谐波谱（HHG）

## 5.1 弱场：介电函数（线性响应）

### 输入文件准备

使用第二章中的 STRU 和 KPT，仅修改 INPUT 如下：

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_linear
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         4000          # 80 fs，足够分辨率
md_dt            0.02
md_tfirst        0

#Parameters (4.外场：宽带弱场激发)
td_vext          1
td_stype         2             # 混合规范

td_tstart        1
td_tend          4000

td_vext_dire     1             # x 方向
td_ttype         0
td_gauss_amp     0.001         # 弱场振幅 (V/Å)，线性响应区间
td_gauss_t0      50
td_gauss_sigma   0.5           # 较宽的高斯，覆盖宽频域
td_gauss_freq    0.0
td_gauss_phase   0.0

#Parameters (5.输出)
out_current      1             # 周期体系：输出电流
out_efield       1
```

确保场强足够弱（`td_gauss_amp ≤ 0.001 V/Å`），使体系响应处于线性区间。

### 运行

先跑 SCF，再跑 TDDFT：
```bash
# 第一步：SCF
cp -r scf_input/ tddft_linear/
cd tddft_linear/
# 替换 INPUT 为 TDDFT 版本
mpirun -n 8 abacus
```

### 后处理：电流 → 介电函数

TDDFT 演化结束后，输出的 `SPIN1_CURRENT` 文件包含各时刻的电流密度 $J(t)$。

电流 $J(t)$ 与电偶极矩的时间导数成正比：$J(t) \propto \dot{D}(t)$，因此傅里叶变换时：

$$D(\omega) \propto \frac{J(\omega)}{i\omega}$$

介电函数的虚部（对应光吸收）由此得到：

$$\varepsilon_2(\omega) \propto \frac{\text{Im}[J(\omega)]}{\omega \cdot E(\omega)}$$

使用 ABACUS 的后处理脚本：

```bash
python tools/plot-tools/plot_current.py \
    --current SPIN1_CURRENT \
    --efield  efield_0.dat \
    --output  dielectric_Si.dat
```

### 结果与分析

Si 的介电函数有两个主要特征峰，分别位于约 3.4 eV 和 4.3 eV，对应 E₁ 和 E₂ 带间跃迁。

| 规范 | 介电函数精度 | 说明 |
|------|------------|------|
| 混合规范（`td_stype = 2`） | 准确 | 峰位、峰高与参考值吻合 |
| 速度规范（`td_stype = 1`） | 偏低且展宽 | 小基组下系统性低估 |

赵昊天等的计算表明，速度规范在 Si 原胞计算中存在明显的系统性误差，通过经验修正方案可部分缓解，但混合规范无需任何修正即可稳定地给出准确且物理一致的结果。

---

## 5.2 强场：高次谐波谱（HHG，非线性响应）

当场强增大到非线性区间（0.001–0.01 V/Å 量级），体系响应不再线性，电流信号中出现基频的高次谐波成分。

### 强场 INPUT

将线性响应 INPUT 中的场强参数修改如下（其余不变）：

```
#Parameters (4.外场：强场激发)
td_vext          1
td_stype         2             # 混合规范

td_tstart        1
td_tend          4000

td_vext_dire     1
td_ttype         0
td_gauss_amp     0.01          # 强场振幅，比线性响应大 10 倍
td_gauss_t0      50
td_gauss_sigma   0.5
td_gauss_freq    0.372         # 约 1.0 eV，基频激光频率
td_gauss_phase   0.0
```

> **参数说明**：`td_gauss_freq` 设置激光的中心频率（单位 eV/ℏ，具体单位请参考 ABACUS 在线文档）。高次谐波谱通过对 Si 单胞施加一系列不同场强（0.001、0.003、0.005、0.01 V/Å）的高斯脉冲来测试。

### 后处理：电流 → 高次谐波谱

HHG 谱通过对电流信号 $J(t)$ 作傅里叶变换得到：

$$I_{\text{HHG}}(\omega) \propto |\omega^2 \cdot J(\omega)|^2$$

对 $J(\omega)$ 取绝对值平方（强度谱），在频率轴上以基频 $\omega_0$ 的倍数标注，即可得到 HHG 谱。

```bash
python tools/plot-tools/plot_hhg.py \
    --current SPIN1_CURRENT \
    --efield  efield_0.dat \
    --fundamental_freq 0.372 \
    --output  hhg_Si.dat
```

### 5.3 场强从弱到强：谐波阶次的演化

赵昊天等对 Si 单胞进行了系统的强场测试，将高斯脉冲场强从弱到强递增：

| 场强 (V/Å) | 响应类型 | 出现的谐波阶次 |
|-----------|---------|--------------|
| ~0.001 | 线性响应 | 仅 1 次（基频） |
| ~0.003 | 弱非线性 | 1、3 次 |
| ~0.005 | 非线性 | 1、3、5 次 |
| ~0.010 | 强非线性 | 1、3、5、7 次 |

实验规律清晰：随场强增大，HHG 谱中依次出现 1、3、5、7 次谐波信号（奇次谐波，源于 Si 的反演对称性）。这一结果与 HHG 的理论预期完全吻合，也证明混合规范在强场和弱场条件下均能给出可靠结果。

> **混合规范的计算效率优势**：相比速度规范，混合规范在处理非局域势时更高效，在相同基组下精度更高。这意味着可以使用较小的基组（如 2s2p1d）完成高精度的 HHG 计算，而速度规范为达到同等精度则需要更大的基组（计算成本显著增加）。
# 六、常见问题与注意事项

## 6.1 计算无法启动或报错

**问题：** 程序报错提示 `gamma_only` 不兼容。

**原因：** `gamma_only = 1` 只存储实数波函数，rt-TDDFT 需要复数传播子，不兼容。

**解决：** 在 INPUT 中设置 `gamma_only 0`，且从 SCF 阶段就必须设置，以保证输出的电荷密度格式正确。

---

**问题：** 报错提示 `basis_type` 必须为 `lcao`。

**原因：** ABACUS 的 rt-TDDFT 目前仅支持 LCAO（数值原子轨道）基组，不支持平面波基组。

**解决：** 设置 `basis_type lcao`，并在 STRU 文件中提供 `NUMERICAL_ORBITAL` 轨道文件。

---

## 6.2 演化时间不够长，吸收谱分辨率低

**问题：** 吸收谱峰型很宽，分辨不清相邻峰。

**原因：** 频率分辨率 $\Delta\omega \approx 1/T_{\text{total}}$，演化时间 $T_{\text{total}}$ 越短，分辨率越差。1000 步 × 0.02 fs = 20 fs 的演化时间只能分辨 ~0.2 eV。

**解决：** 增大 `md_nstep`，使总演化时间达到所需分辨率。经验建议：
- 分辨 ~0.1 eV 的峰：演化时间 ≥ 40 fs
- 分辨 ~0.05 eV 的峰：演化时间 ≥ 80 fs

---

## 6.3 速度规范在低频出现发散

**问题：** 使用 `td_stype = 1`（速度规范）时，吸收谱或介电函数在低频区出现明显上翘或发散。

**原因：** 周期体系中，电流 $J(\omega)$ 除以 $\omega$ 得到等效偶极，低频时分母趋近于零，NAO 基组下速度规范的系统性误差在低频被放大。

**解决：** 切换到混合规范 `td_stype = 2`。这是根本解决方案，不需要调整其他参数。

---

## 6.4 电场设置相关

**同时施加多个电场方向**

如需同时激发多个方向（例如计算全张量介电函数），可以在 `td_vext_dire` 中指定多个方向：

```
td_vext_dire    1 2 3    # 同时沿 x、y、z 方向施加
```

相应的 `td_gauss_amp`、`td_gauss_t0` 等参数也需提供等数量的值。

**电场频率单位**

`td_gauss_freq` 的单位为 eV/ℏ，换算关系：1 eV/ℏ ≈ 1.52 × 10¹⁵ rad/s。对于宽带激发（线性响应），设置 `td_gauss_freq 0.0` 并配合小 sigma 即可。

---

## 6.5 td_tend 的设置

`td_tend` 设置电场结束的步数，之后标势清零，矢势固定在 `td_tend` 时刻的值。对于吸收谱计算，可以令 `td_tend = md_nstep`，全程施加电场（脉冲型电场会自然衰减到零）。也可以提前关闭以节约后续无场演化的资源——但需确保脉冲已经完整通过（通常要求 `td_tend` > `td_gauss_t0` + 3σ）。
# 附录

## A. 参考资料

1. 赵昊天, 何力新. *Hybrid-Gauge Real-Time Time-Dependent Density Functional Theory in the Numerical Atomic Orbital Basis.* J. Chem. Theory Comput. **2025**, 21, 3335–3341. DOI: 10.1021/acs.jctc.5c00111

2. ABACUS RT-TDDFT 官方教程：https://mcresearch.github.io/abacus-user-guide/abacus-tddft.html

3. ABACUS 输入参数文档（TDDFT 部分）：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#tddft-time-dependent-density-functional-theory

4. ABACUS 算例仓库（TDDFT 部分）：https://gitee.com/mcresearch/abacus-user-guide/tree/master/examples/tddft

5. Pemmaraju C D, et al. *Velocity-Gauge Real-Time TDDFT Within a Numerical Atomic Orbital Basis Set.* Comput. Phys. Commun. **2018**, 226, 30–38. DOI: 10.1016/j.cpc.2018.01.013

## B. 关键参数汇总

| 参数 | 用途 | 可选值 |
|------|------|--------|
| `esolver_type` | 求解器类型 | `tddft`（必须） |
| `basis_type` | 基组类型 | `lcao`（必须） |
| `gamma_only` | 是否 Gamma 点简化 | `0`（必须） |
| `calculation` | 计算类型 | `md`（必须） |
| `td_stype` | 规范选择 | `0`/`1`/`2` |
| `td_vext` | 外场开关 | `1` |
| `td_ttype` | 电场波包形状 | `0`=高斯 |
| `out_current` | 电流输出 | `1`（周期体系） |
| `out_dipole` | 偶极矩输出 | `1`（非周期体系） |

## C. 进阶学习方向

- **含 Ehrenfest 离子运动的 RT-TDDFT**：设置非零 `md_tfirst`，允许离子在电场下运动，模拟光激发驱动的结构动力学
- **多 k 点下的 HHG 计算**：增大 k 点网格以提升 HHG 谱的 k 空间分辨率
- **更高阶谐波**：增大场强或延长演化时间，观察 9 次以上谐波
- **非平衡态输运**：利用 RT-TDDFT 计算强场下的电荷转移和电流响应
