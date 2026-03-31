# ABACUS 态密度（DOS/PDOS）计算教程

## 前言

### 教程目标

本教程将帮助你掌握使用 ABACUS 计算和分析材料的态密度（DOS）和投影态密度（PDOS）。通过本教程，你将学会：

- 理解态密度的物理意义及其在材料研究中的应用
- 掌握 ABACUS 中 DOS/PDOS 计算的两步流程
- 熟悉关键 INPUT 参数的设置方法
- 使用 abacustest 工具进行结果后处理和可视化
- 分析和解读 DOS/PDOS 图谱

### 适用读者

本教程是 ABACUS 系列教程的第4篇，适合已完成前3篇教程学习的读者：

- 第1篇：ABACUS 基本介绍（INPUT/STRU/KPT 文件格式、结构优化）
- 第2篇：力学性质计算（弹性常数）
- 第3篇：能带结构计算

如果你已经熟悉 ABACUS 的基本操作和 SCF/NSCF 计算流程，可以直接学习本教程。

### 前置知识

学习本教程前，你需要了解：

- ABACUS 的基本文件结构（INPUT、STRU、KPT）
- 自洽场（SCF）和非自洽场（NSCF）计算的概念
- 基本的 Linux 命令行操作
- abacustest 工具的基本使用

### 本篇内容

本教程包含以下内容：

**第一章：态密度的物理基础**
介绍 DOS 和 PDOS 的定义、物理意义，以及它们与能带结构的关系。

**第二章：ABACUS 的 DOS/PDOS 计算方法**
详细讲解 ABACUS 中 DOS/PDOS 计算的两步流程和关键参数设置。

**第三章：MgO 态密度计算完整实战**
通过 MgO（氧化镁）案例，演示从输入文件准备到计算执行的完整流程。

**第四章：使用 abacustest 后处理**
介绍如何使用 abacustest 工具绘制和分析 DOS/PDOS 图谱。

**第五章：结果分析与常见问题**
讲解如何解读 DOS/PDOS 图，以及常见问题的排查方法。

### 案例说明

本教程使用 MgO（氧化镁）作为示例材料。MgO 是典型的离子晶体，具有岩盐结构，是宽带隙绝缘体。选择 MgO 的原因：

- 结构简单（2 原子），便于理解
- 电子结构特征明显：O 2p 轨道主导价带，Mg 3s 轨道贡献导带
- PDOS 可清晰展示不同元素和轨道的贡献
- 计算量适中，适合教学演示

通过 MgO 案例，你将学会如何计算和分析离子晶体的电子态密度。
# 第一章：态密度的物理基础

## 1.1 态密度（DOS）的定义与意义

### 什么是态密度

态密度（Density of States，简称 DOS）是指在能量为 E 附近单位能量间隔内可供电子占据的电子状态数目。数学上，DOS 定义为：

```
DOS(E) = dN/dE
```

其中 N 是能量小于 E 的电子态总数。

### 物理意义

态密度描述了电子态在能量空间的分布密度。通过 DOS，我们可以了解：

**材料的导电性**
- 金属：费米能级处 DOS 较高，电子容易激发，导电性好
- 半导体：费米能级位于带隙中，DOS 为零，导电性较差
- 绝缘体：费米能级位于较大带隙中，DOS 为零，几乎不导电

**电子态的分布**
- DOS 的峰值表示该能量附近电子态密集
- DOS 的谷值或零点表示该能量附近电子态稀疏或不存在
- 费米能级附近的 DOS 决定了材料的低能激发性质

**材料的光学和磁学性质**
- 光学跃迁强度与初态和终态的 DOS 相关
- 磁性材料的自旋极化 DOS 反映磁矩分布

### DOS 与材料性质的关系

不同类型材料的 DOS 特征：

| 材料类型 | 费米能级位置 | 费米能级处 DOS | 带隙 | 导电性 |
|---------|------------|--------------|------|--------|
| 金属 | 导带内 | 高 | 无 | 高 |
| 半导体 | 带隙中 | 零 | 小（<3 eV） | 中等 |
| 绝缘体 | 带隙中 | 零 | 大（>3 eV） | 低 |

## 1.2 DOS 与能带结构的关系

### 能带结构：E-k 关系

能带结构描述了材料中电子能量 E 与波矢 k 之间的关系，即 E(k)。能带结构揭示：

- 电子能级的分布
- 能带的宽度和色散关系
- 价带和导带之间的带隙

能带结构是在倒空间（k 空间）中描述电子态。

### DOS：对 k 空间积分

DOS 是通过对能带结构在整个布里渊区进行积分得到的：

```
DOS(E) = (1/V) Σ_n ∫_BZ δ(E - E_n(k)) dk
```

其中：
- V 是晶胞体积
- n 是能带指标
- BZ 表示布里渊区
- δ 是 Dirac delta 函数

这个积分的物理意义是：统计所有 k 点上能量为 E 的电子态数目。

### 从能带到 DOS 的转换

能带结构和 DOS 提供了互补的信息：

**能带结构的优势：**
- 显示电子态的 k 空间分布
- 揭示能带的色散关系
- 可以识别直接带隙和间接带隙

**DOS 的优势：**
- 统计所有 k 点的贡献，信息更全面
- 直接反映电子态的能量分布
- 便于计算热力学和输运性质

### 范霍夫奇点与 DOS 峰值

范霍夫奇点（Van Hove singularity）是能带结构中 ∇_k E(k) = 0 的点，即能带的极值点（最大值、最小值或鞍点）。

在范霍夫奇点附近：
- 能带较平坦，电子群速度接近零
- 大量 k 点对应相近的能量
- DOS 出现峰值或奇异性

DOS 的峰值通常对应能带结构中的范霍夫奇点，这些峰值在材料性质中起重要作用。

## 1.3 投影态密度（PDOS）

### PDOS 的定义

投影态密度（Partial DOS 或 Projected DOS，简称 PDOS）是将电子波函数投影到特定原子轨道上得到的态密度。PDOS 回答的问题是：某个能量处的电子态有多少贡献来自特定的原子或轨道。

数学上，PDOS 定义为：

```
PDOS_i(E) = Σ_n ∫_BZ |⟨ψ_n(k)|φ_i⟩|² δ(E - E_n(k)) dk
```

其中：
- ψ_n(k) 是第 n 条能带在 k 点的波函数
- φ_i 是第 i 个原子轨道
- |⟨ψ_n(k)|φ_i⟩|² 是投影系数

### 投影方式

PDOS 可以按不同方式投影：

**按元素投影**
- 将所有属于某元素的原子轨道贡献加和
- 例如：MgO 中的 Mg-PDOS 和 O-PDOS
- 用途：分析不同元素对电子态的贡献

**按轨道壳层投影**
- 按角量子数 l 分类：s (l=0)、p (l=1)、d (l=2)、f (l=3)
- 例如：O 的 2s-PDOS 和 2p-PDOS
- 用途：分析不同轨道类型的贡献

**按具体轨道投影**
- 按磁量子数 m 进一步细分
- 例如：p 轨道分为 p_x、p_y、p_z
- 用途：分析轨道取向和对称性

**按原子投影**
- 投影到特定原子的所有轨道
- 例如：第 1 个 O 原子的 PDOS
- 用途：分析局域电子态分布

### PDOS 的应用

**成键分析**
- 价带区域的 PDOS 重叠表示成键轨道
- 导带区域的 PDOS 重叠表示反键轨道
- PDOS 峰值的能量差反映轨道杂化强度

**轨道贡献分析**
- 识别价带顶和导带底的主要轨道成分
- 例如：过渡金属氧化物中 d 轨道对带隙的贡献
- 指导材料设计和性质调控

**化学键性质判断**
- 离子键：不同元素的 PDOS 分离明显
- 共价键：不同元素的 PDOS 重叠显著
- 金属键：费米能级处 PDOS 连续

### 总 DOS 与 PDOS 的关系

总 DOS 是所有 PDOS 的加和：

```
DOS(E) = Σ_i PDOS_i(E)
```

通过 PDOS，我们可以将总 DOS 分解为不同原子或轨道的贡献，从而深入理解电子结构。

## 1.4 费米能级的作用

### 费米能级的定义

费米能级（Fermi level）是绝对零度时电子占据的最高能级。在有限温度下，费米能级是化学势，定义为：

```
μ = (∂F/∂N)_T,V
```

其中 F 是自由能，N 是电子数。

在 DOS 图中，通常将费米能级设为能量零点，负能量表示费米能级以下（占据态），正能量表示费米能级以上（未占据态）。

### 费米能级的物理意义

**电子占据的分界线**
- 零温时：费米能级以下的态全部占据，以上的态全部空置
- 有限温度：费米能级附近的占据遵循费米-狄拉克分布

**材料类型的判据**
- 金属：费米能级穿过能带，DOS(E_F) > 0
- 半导体/绝缘体：费米能级位于带隙中，DOS(E_F) = 0

**电学性质的决定因素**
- 电导率与费米能级处的 DOS 成正比
- 费米能级附近的电子态参与低能激发过程

### 费米能级附近 DOS 的重要性

费米能级附近的 DOS 决定了材料的许多重要性质：

**电导率**
- 金属的电导率 σ ∝ DOS(E_F)
- 费米能级处 DOS 越高，电导率越高

**比热**
- 电子比热 C_e ∝ DOS(E_F) × T
- 低温下电子比热与费米能级处 DOS 成正比

**磁化率**
- 泡利顺磁磁化率 χ ∝ DOS(E_F)
- 费米能级处 DOS 高的材料磁化率大

**超导转变温度**
- BCS 理论：T_c ∝ exp(-1/[V×DOS(E_F)])
- 费米能级处 DOS 影响超导配对

### 不同材料的费米能级位置

**金属（如 Al、Cu）**
- 费米能级位于导带内
- DOS(E_F) 较高
- 价带和导带连续，无带隙

**半导体（如 Si、GaAs）**
- 费米能级位于带隙中
- 本征半导体：费米能级接近带隙中心
- n 型半导体：费米能级靠近导带底
- p 型半导体：费米能级靠近价带顶

**绝缘体（如 MgO、SiO₂）**
- 费米能级位于较大带隙中
- 带隙通常大于 3 eV
- 价带顶和导带底的 DOS 峰值远离费米能级

通过分析费米能级的位置和附近的 DOS 分布，我们可以判断材料的基本电学性质，这是 DOS 计算的重要应用之一。
# 第二章：ABACUS 的 DOS/PDOS 计算方法

## 2.1 计算流程：两步法

### 为什么需要两步计算

ABACUS 中计算 DOS/PDOS 采用两步法：

**第一步：SCF 自洽计算**
- 目的：求解基态电子密度
- 过程：迭代求解 Kohn-Sham 方程，直到电荷密度收敛
- 输出：收敛的电荷密度文件（SPIN1_CHG）

**第二步：NSCF 非自洽计算**
- 目的：在固定电荷密度下计算更密 k 点网格的本征值
- 过程：读入电荷密度，求解 Kohn-Sham 方程一次（不迭代）
- 输出：DOS 和 PDOS 数据文件

### 为什么不能一步完成

分两步的原因：

1. **k 点需求不同**
   - SCF：需要足够 k 点保证电荷密度收敛，但不需要太密
   - NSCF：需要更密的 k 点保证 DOS 曲线平滑

2. **计算效率**
   - SCF 每步都需要迭代，k 点太密会大幅增加计算量
   - NSCF 只计算一次，可以使用更密的 k 点

3. **灵活性**
   - 一次 SCF 可以用于多次 NSCF（不同 k 点、不同能量范围）
   - 便于参数调试和结果分析

### 计算流程图

```
准备输入文件
    ↓
SCF 自洽计算
    ├─ INPUT: calculation=scf, out_chg=1
    ├─ KPT: 常规 k 点网格
    └─ 输出: SPIN1_CHG
    ↓
修改 INPUT 参数
    ├─ calculation=nscf
    ├─ init_chg=file
    └─ out_dos=1 或 2
    ↓
NSCF 非自洽计算
    ├─ KPT: 更密的 k 点网格
    └─ 输出: DOS1, PDOS (如果 out_dos=2)
    ↓
后处理与可视化
```

## 2.2 关键 INPUT 参数详解

### out_dos 参数

`out_dos` 控制是否输出 DOS 以及输出类型。

| 取值 | 含义 | 输出文件 | 适用基组 |
|------|------|---------|---------|
| 0 | 不输出 DOS | 无 | 所有 |
| 1 | 输出总 DOS | DOS1, TDOS | 所有 |
| 2 | 输出总 DOS 和 PDOS | DOS1, TDOS, PDOS | 仅 LCAO |

**注意事项：**
- PDOS 只能在 LCAO 基组下输出（`basis_type = lcao`）
- 平面波基组（`basis_type = pw`）只能输出总 DOS（`out_dos = 1`）
- 原因：PDOS 需要投影到局域原子轨道，平面波基组没有局域轨道

**推荐设置：**
- LCAO 基组：`out_dos = 2`（同时获得 DOS 和 PDOS）
- 平面波基组：`out_dos = 1`（只能输出总 DOS）

### init_chg 参数

`init_chg` 控制初始电荷密度的来源。

| 取值 | 含义 | 适用计算类型 |
|------|------|------------|
| atomic | 从原子电荷密度叠加 | SCF |
| file | 从文件读取电荷密度 | NSCF |

**SCF 计算：**
```
calculation    scf
init_chg       atomic
out_chg        1
```
- `init_chg = atomic`：从原子电荷密度开始迭代
- `out_chg = 1`：输出收敛的电荷密度到 SPIN1_CHG

**NSCF 计算：**
```
calculation    nscf
init_chg       file
out_dos        2
```
- `init_chg = file`：读取 SCF 输出的 SPIN1_CHG
- 电荷密度固定，不再迭代

### calculation 参数

`calculation` 指定计算类型。

| 取值 | 含义 | 是否迭代 | 用途 |
|------|------|---------|------|
| scf | 自洽场计算 | 是 | 求解基态电荷密度 |
| nscf | 非自洽场计算 | 否 | 计算 DOS/能带 |

**参数组合：**
- SCF：`calculation = scf` + `init_chg = atomic`
- NSCF：`calculation = nscf` + `init_chg = file`

### symmetry 参数

`symmetry` 控制是否使用对称性简化 k 点。

| 取值 | 含义 | 推荐使用场景 |
|------|------|------------|
| -1 | 不进行对称性分析 | 调试 |
| 0 | 仅考虑时间反演对称性 | NSCF |
| 1 | 进行对称性分析 | SCF |

**SCF 计算：**
- `symmetry = 1`：利用对称性减少 k 点数量，提高效率

**NSCF 计算：**
- `symmetry = 0`：不使用对称性，保留完整 k 点网格
- 原因：DOS 需要完整的 k 点采样，对称性简化会导致 DOS 不准确

## 2.3 影响 DOS 质量的参数

### k 点密度

k 点密度直接影响 DOS 曲线的平滑度。

| k 点网格 | DOS 平滑度 | 计算量 | 推荐场景 |
|---------|-----------|--------|---------|
| 6×6×6 | 粗糙，有尖峰 | 低 | 快速测试 |
| 12×12×12 | 较平滑 | 中等 | 一般计算 |
| 18×18×18 | 平滑 | 较高 | 精确计算 |
| 24×24×24 | 非常平滑 | 高 | 高精度要求 |

**选择原则：**
- NSCF 的 k 点密度应高于 SCF
- 金属体系需要更密的 k 点（费米面附近态密度变化快）
- 绝缘体可以使用相对稀疏的 k 点

### smearing 相关参数

smearing 用于处理费米面附近的部分占据态，影响 SCF 收敛和 DOS 平滑度。

**smearing_method**

| 方法 | 适用体系 | 特点 |
|------|---------|------|
| gauss | 半导体、绝缘体 | 高斯展宽，平滑 |
| mp | 金属 | Methfessel-Paxton，减少展宽误差 |
| fd | 金属 | Fermi-Dirac，物理意义明确 |

**smearing_sigma**
- 单位：Rydberg
- 典型值：0.002-0.015 Ry（约 0.027-0.2 eV）
- 作用：控制展宽宽度
- 注意：仅影响 SCF 收敛，不影响 NSCF 的 DOS 输出

**dos_sigma**
- 单位：eV
- 作用：仅用于绘制 DOS 图时的高斯展宽
- 注意：不影响 SCF 或 NSCF 计算结果

**smearing_sigma vs dos_sigma 的区别：**

| 参数 | 作用阶段 | 影响对象 | 单位 |
|------|---------|---------|------|
| smearing_sigma | SCF 计算 | 电荷密度收敛 | Ry |
| dos_sigma | 后处理 | DOS 图平滑度 | eV |

## 2.4 LCAO 与平面波基组的差异

### 基组对 DOS 计算的影响

| 特性 | LCAO 基组 | 平面波基组 |
|------|----------|-----------|
| 总 DOS | ✓ 支持 | ✓ 支持 |
| PDOS | ✓ 支持 | ✗ 不支持 |
| 计算速度 | 快 | 慢 |
| 精度 | 依赖轨道质量 | 依赖 ecutwfc |

### 为什么平面波不能输出 PDOS

**LCAO 基组：**
- 波函数展开在局域原子轨道上：ψ = Σ c_i φ_i
- 投影系数 c_i 直接给出轨道贡献
- 可以自然地计算 PDOS

**平面波基组：**
- 波函数展开在平面波上：ψ = Σ c_k e^(ik·r)
- 平面波是非局域的，无法直接投影到原子轨道
- 需要额外的投影算符（如 PAW 方法），ABACUS 目前未实现

### 如何选择基组

**需要 PDOS 分析：**
- 必须使用 LCAO 基组
- 设置 `basis_type = lcao`
- 设置 `out_dos = 2`

**只需要总 DOS：**
- 可以使用平面波或 LCAO
- 平面波精度更高但速度慢
- LCAO 速度快但需要选择合适的轨道文件

## 2.5 PDOS 文件格式

### XML 文件结构

ABACUS 输出的 PDOS 文件是 XML 格式，位于 `OUT.suffix/PDOS`。

**文件头部：**
```xml
<pdos>
<nspin>1</nspin>
<norbitals>14</norbitals>
<energy_values units="eV">
0.00000
0.01000
0.02000
...
</energy_values>
```

- `<nspin>`：自旋通道数（1=非自旋极化，2=自旋极化）
- `<norbitals>`：轨道总数
- `<energy_values>`：能量网格点（横坐标）

### 轨道编号规则

每个轨道由以下量子数标识：

**l（角量子数）**
- l=0：s 轨道
- l=1：p 轨道
- l=2：d 轨道
- l=3：f 轨道

**m（磁量子数）**
- 取值范围：0 到 2l
- s 轨道：m=0（1 个）
- p 轨道：m=0,1,2（3 个，对应 p_z, p_x, p_y）
- d 轨道：m=0,1,2,3,4（5 个）

**z（径向轨道编号）**
- 表示同一角量子数的不同径向函数
- 例如：DZP 基组的 s 轨道有 z=1,2（双 zeta）

### 轨道数量示例

**Mg 的 2s1p 基组：**
- 2 个 s 轨道：l=0, m=0, z=1,2
- 3 个 p 轨道：l=1, m=0,1,2, z=1
- 共 5 个轨道

**O 的 2s2p1d 基组：**
- 2 个 s 轨道：l=0, m=0, z=1,2
- 6 个 p 轨道：l=1, m=0,1,2, z=1,2
- 5 个 d 轨道：l=2, m=0,1,2,3,4, z=1
- 共 13 个轨道

### 轨道数据格式

```xml
<orbital
index="1"
atom_index="1"
species="Mg"
l="0"
m="0"
z="1"
>
<data>
0.00000000 0.00000000
0.00123456 0.00123456
...
</data>
</orbital>
```

- `index`：轨道全局编号
- `atom_index`：原子编号（从 1 开始）
- `species`：元素符号
- `l, m, z`：量子数
- `<data>`：DOS 数据（两列对应自旋上/下）

### 数据提取

PDOS 文件可以手动解析，但推荐使用 abacustest 工具自动处理，将在第四章详细介绍。
# 第三章：MgO 态密度计算完整实战

## 3.1 案例背景与目标

### MgO 材料介绍

MgO（氧化镁）是典型的离子晶体，具有以下特点：

**晶体结构**
- 岩盐结构（NaCl 型）
- 面心立方晶格
- Mg²⁺ 和 O²⁻ 交替排列

**电子结构特征**
- 宽带隙绝缘体，实验带隙约 7.8 eV
- 强离子键：Mg 失去 2 个电子，O 获得 2 个电子
- 价带主要由 O 2p 轨道组成
- 导带主要由 Mg 3s 轨道组成

**为什么选择 MgO**
- 结构简单，只有 2 个原子
- 电子结构特征明显，便于理解
- PDOS 可清晰展示离子键特征
- 计算量适中，适合教学

### 计算目标

本案例将完成以下计算：

1. **SCF 自洽计算**：获得 MgO 的基态电荷密度
2. **NSCF 态密度计算**：计算总 DOS 和 PDOS
3. **结果分析**：
   - 确定带隙大小
   - 分析 Mg 和 O 的轨道贡献
   - 验证离子键特征

### 预期结果

**总 DOS 特征**
- 带隙约 4-5 eV（PBE 泛函低估带隙）
- 价带宽度约 5-6 eV
- 导带起始于费米能级以上约 4-5 eV

**PDOS 特征**
- 价带顶：O 2p 轨道主导
- 导带底：Mg 3s 轨道主导
- Mg 和 O 的 PDOS 分离明显（离子键特征）

## 3.2 输入文件详解

### 3.2.1 INPUT 文件

以下是完整的 INPUT 文件及逐参数说明：

```
INPUT_PARAMETERS

# System variables
calculation         scf          # 自洽场计算
symmetry            0            # 不使用对称性
kspacing            0.08         # k点间距，单位 1/Bohr
precision           double       # 双精度计算

# Plane wave related variables
ecutwfc             100          # 平面波截断能，单位 Ry

# Electronic structure
basis_type          lcao         # LCAO基组（用于输出PDOS）
ks_solver           genelpa      # Kohn-Sham方程求解器
smearing_method     gauss        # 高斯展宽
smearing_sigma      0.015        # 展宽参数，单位 Ry
mixing_type         broyden      # Broyden混合方法
mixing_beta         0.8          # 混合参数
scf_nmax            100          # 最大SCF迭代步数
scf_thr             1e-07        # SCF收敛阈值

# DOS output
init_chg            file         # 从文件读取电荷密度（NSCF用）
out_dos             2            # 输出总DOS和PDOS
```

**参数详细说明：**

| 参数 | 取值 | 说明 |
|------|------|------|
| calculation | scf | SCF自洽计算，求解基态电荷密度 |
| symmetry | 0 | 不使用对称性（NSCF需要完整k点） |
| kspacing | 0.08 | k点间距，自动生成k点网格 |
| precision | double | 双精度计算，提高数值精度 |
| ecutwfc | 100 | 平面波截断能，与轨道文件推荐值一致 |
| basis_type | lcao | LCAO基组，可输出PDOS |
| ks_solver | genelpa | 高效的本征值求解器 |
| smearing_method | gauss | 高斯展宽，适合绝缘体 |
| smearing_sigma | 0.015 | 展宽参数（约0.2 eV） |
| mixing_type | broyden | Broyden混合，加速SCF收敛 |
| mixing_beta | 0.8 | 混合系数 |
| scf_nmax | 100 | 最大迭代步数 |
| scf_thr | 1e-07 | 电荷密度收敛阈值 |
| init_chg | file | 从文件读取电荷密度（NSCF时使用） |
| out_dos | 2 | 输出总DOS和PDOS |

**注意：**
- SCF 计算时，`init_chg` 应设为 `atomic`，`out_dos` 可以不设置
- NSCF 计算时，需要将 `calculation` 改为 `nscf`，`init_chg` 改为 `file`

### 3.2.2 STRU 文件

以下是完整的 STRU 文件：

```
ATOMIC_SPECIES
Mg 24.305000 Mg.PD04.PBE.UPF
O 15.999400 O.upf

NUMERICAL_ORBITAL
Mg_gga_10au_100Ry_2s1p.orb
O_gga_6au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    2.97691954880     0.00000000000     0.00000000000
    1.48845977440     2.57808795428     0.00000000000
    1.48845977440     0.85936265143     2.43064463329

ATOMIC_POSITIONS
Cartesian

Mg
0.000000
1
    0.00000000000     0.00000000000     0.00000000000 1 1 1

O
0.000000
1
    2.97691954880     1.71872530285     1.21532231664 1 1 1
```

**结构说明：**

| 部分 | 内容 | 说明 |
|------|------|------|
| ATOMIC_SPECIES | Mg, O | 元素符号、原子量、赝势文件 |
| NUMERICAL_ORBITAL | 轨道文件 | Mg: 2s1p, O: 2s2p1d |
| LATTICE_CONSTANT | 1.889726 | Bohr转换因子（Bohr到Angstrom） |
| LATTICE_VECTORS | 3个矢量 | 晶格矢量（单位：LATTICE_CONSTANT） |
| ATOMIC_POSITIONS | Cartesian | 原子坐标（笛卡尔坐标系） |

**赝势文件：**
- `Mg.PD04.PBE.UPF`：Mg的ONCV赝势
- `O.upf`：O的赝势

**轨道文件：**
- `Mg_gga_10au_100Ry_2s1p.orb`：Mg的DZP基组（2个s轨道，1个p轨道）
- `O_gga_6au_100Ry_2s2p1d.orb`：O的TZDP基组（2个s轨道，2个p轨道，1个d轨道）

**晶格信息：**
- 晶格常数：约7.64 Bohr（4.04 Å）
- 结构类型：岩盐结构
- 原子数：2个（1个Mg，1个O）

### 3.2.3 KPT 文件

```
K_POINTS
0
Gamma
18 18 18 0 0 0
```

**参数说明：**
- 第一行：`K_POINTS` 关键字
- 第二行：`0` 表示自动生成k点
- 第三行：`Gamma` 表示Gamma中心网格
- 第四行：`18 18 18` 表示三个方向的k点数，`0 0 0` 表示不平移

**k点密度选择：**
- 18×18×18网格对应约5832个k点（考虑对称性后会减少）
- 对于MgO这样的绝缘体，该密度足以获得平滑的DOS曲线

## 3.3 第一步：SCF 自洽计算

### 准备SCF输入文件

将上述INPUT文件修改为SCF版本：

```
INPUT_PARAMETERS

calculation         scf          # SCF自洽计算
symmetry            1            # 使用对称性加速SCF
init_chg            atomic       # 从原子电荷密度开始
out_chg             1            # 输出电荷密度文件

# 其他参数与上面相同
basis_type          lcao
ecutwfc             100
ks_solver           genelpa
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-07
```

**关键修改：**
- `calculation = scf`
- `symmetry = 1`（利用对称性）
- `init_chg = atomic`
- `out_chg = 1`（输出电荷密度）

### 运行SCF计算

```bash
mpirun -np 8 abacus
```

或根据你的计算资源调整进程数。

### 检查SCF收敛性

查看 `OUT.ABACUS/running_scf.log` 文件，关注以下信息：

**收敛过程：**
```
ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)
GE1    -3.456789e+02  0.000000e+00   2.345e-01  1.234e+00
GE2    -3.457123e+02  -3.340000e-02  5.678e-02  1.123e+00
GE3    -3.457145e+02  -2.200000e-03  1.234e-03  1.098e+00
...
GE10   -3.457146e+02  -1.234567e-08  3.456e-08  1.045e+00
```

**收敛判据：**
- `EDIFF`：能量变化，应逐步减小
- `DRHO`：电荷密度变化，应小于 `scf_thr`（1e-07）
- 当 `DRHO < scf_thr` 时，SCF收敛

### 确认输出文件

检查 `OUT.ABACUS/` 目录下是否生成以下文件：

```bash
ls OUT.ABACUS/
```

**必需文件：**
- `SPIN1_CHG`：电荷密度文件（NSCF需要）
- `running_scf.log`：SCF运行日志
- `STRU_READIN_ADJUST.cif`：读入的结构文件

如果 `SPIN1_CHG` 文件存在，说明SCF计算成功，可以进行下一步。

## 3.4 第二步：NSCF 态密度计算

### 修改INPUT文件

将INPUT文件修改为NSCF版本，只需修改3个参数：

```
INPUT_PARAMETERS

calculation         nscf         # 改为NSCF
symmetry            0            # 改为0（不使用对称性）
init_chg            file         # 改为file（读取电荷密度）
out_dos             2            # 输出DOS和PDOS

# 其他参数保持不变
basis_type          lcao
ecutwfc             100
ks_solver           genelpa
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-07
```

**SCF vs NSCF 参数对比：**

| 参数 | SCF | NSCF | 说明 |
|------|-----|------|------|
| calculation | scf | nscf | 计算类型 |
| symmetry | 1 | 0 | NSCF需要完整k点 |
| init_chg | atomic | file | NSCF读取电荷密度 |
| out_chg | 1 | 不需要 | NSCF不输出电荷密度 |
| out_dos | 不需要 | 2 | NSCF输出DOS/PDOS |

### 运行NSCF计算

```bash
mpirun -np 8 abacus
```

NSCF计算通常比SCF快，因为只计算一次，不需要迭代。

### 确认DOS/PDOS文件生成

检查 `OUT.ABACUS/` 目录：

```bash
ls OUT.ABACUS/
```

**应生成的文件：**
- `DOS1`：原始DOS数据
- `DOS1_smearing.dat`：展宽后的DOS数据（常用于绘图）
- `TDOS`：总DOS数据
- `PDOS`：PDOS的XML文件
- `running_nscf.log`：NSCF运行日志

**文件说明：**
- `DOS1_smearing.dat`：三列数据（能量、DOS、积分DOS）
- `TDOS`：总DOS，可用于快速绘图
- `PDOS`：XML格式，包含所有轨道的PDOS信息

## 3.5 提取费米能级

### 使用grep命令提取

费米能级记录在运行日志中，可以用grep命令提取：

```bash
grep -i 'efermi' OUT.ABACUS/running_nscf.log
```

**输出示例：**
```
EFERMI = 8.27771540465 eV
```

或者从SCF日志中提取（结果相同）：

```bash
grep -i 'efermi' OUT.ABACUS/running_scf.log
```

### 费米能级的物理意义

对于MgO这样的绝缘体：

**费米能级位置**
- 位于带隙中间
- 价带顶在费米能级以下约2-3 eV
- 导带底在费米能级以上约2-3 eV

**在DOS图中的作用**
- 通常将费米能级设为能量零点
- 负能量：占据态（价带）
- 正能量：未占据态（导带）
- 费米能级处DOS为零（带隙）

### 费米能级在后处理中的应用

绘制DOS图时，需要将能量相对于费米能级平移：

```
E_plot = E_raw - E_fermi
```

abacustest工具会自动完成这一步，无需手动处理。

至此，MgO的DOS/PDOS计算已全部完成，下一章将介绍如何使用abacustest进行结果后处理和可视化。
# 第四章：使用 abacustest 后处理

## 4.1 abacustest dos-pdos 命令详解

### 基本语法

```bash
abacustest model dos-pdos -j <job_path> [options]
```

### 主要参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-j, --job` | ABACUS计算目录路径 | 必需 |
| `--range` | 能量范围（相对费米能级），单位eV | -10 10 |
| `--plot-type` | 绘图类型 | species |
| `--atom-index` | 原子编号（1-based） | 所有原子 |
| `--suffix` | 输出文件后缀 | 无 |
| `--no-save-data` | 不保存数据文件 | False |
| `--no-save-plot` | 不保存图片文件 | False |

### 绘图类型说明

| plot-type | 说明 | 示例 |
|-----------|------|------|
| species | 按元素投影 | Mg、O的PDOS |
| shell | 按轨道壳层投影 | Mg的s、p轨道 |
| orbital | 按具体轨道投影 | O的p_x、p_y、p_z |
| atom | 按原子编号投影 | 第1个原子的PDOS |

### 输出文件

| 文件名 | 内容 | 格式 |
|--------|------|------|
| DOS.dat | 总DOS数据 | 文本 |
| DOS.png | 总DOS图 | 图片 |
| PDOS.dat | PDOS数据 | 文本 |
| PDOS.png | PDOS图 | 图片 |

## 4.2 绘制总 DOS 图

### 基本用法

进入计算目录，执行：

```bash
cd <计算目录>
abacustest model dos-pdos -j .
```

这将自动：
1. 读取 `OUT.ABACUS/TDOS` 和 `OUT.ABACUS/PDOS` 文件
2. 从 `running_nscf.log` 提取费米能级
3. 生成 `DOS.png` 和 `PDOS.png`

### 自定义能量范围

聚焦费米能级附近：

```bash
abacustest model dos-pdos -j . --range -5 7
```

这将绘制费米能级以下5 eV到以上7 eV的范围。

### DOS 图的特征

**MgO 的总 DOS 图应显示：**
- 价带区域（负能量）：-6 eV 到 0 eV
- 带隙：约 4-5 eV（PBE泛函低估）
- 导带区域（正能量）：4 eV 以上
- 费米能级处 DOS 为零

## 4.3 绘制 PDOS 图

### 按元素投影

默认按元素投影：

```bash
abacustest model dos-pdos -j . --plot-type species
```

**MgO 的 PDOS 特征：**
- O 的 PDOS 主导价带（-6 eV 到 0 eV）
- Mg 的 PDOS 主导导带（4 eV 以上）
- 两者在价带区域几乎不重叠（离子键特征）

### 按轨道壳层投影

查看不同轨道的贡献：

```bash
abacustest model dos-pdos -j . --plot-type shell
```

**预期结果：**
- O 2p 轨道主导价带顶
- Mg 3s 轨道主导导带底
- O 2s 轨道在价带深处（约 -20 eV）

### 按具体轨道投影

查看 p 轨道的各个分量：

```bash
abacustest model dos-pdos -j . --plot-type orbital
```

这将显示 p_x、p_y、p_z 的独立贡献。

## 4.4 PDOS 图的物理意义

### 离子键的 PDOS 特征

MgO 的 PDOS 图清晰展示了离子键的特征：

**价带区域（占据态）**
- 主要由 O 2p 轨道组成
- Mg 的贡献极小
- 说明价电子主要定域在 O 原子上

**导带区域（未占据态）**
- 主要由 Mg 3s 轨道组成
- O 的贡献较小
- 说明激发态电子倾向于定域在 Mg 原子上

**PDOS 分离明显**
- Mg 和 O 的 PDOS 在价带区域几乎不重叠
- 表明电子转移完全，形成 Mg²⁺ 和 O²⁻
- 这是典型的离子键特征

### 与共价键的对比

如果是共价键材料（如 Si）：
- 价带区域不同原子的 PDOS 会显著重叠
- 表明电子在原子间共享
- PDOS 峰值能量接近

MgO 的 PDOS 分离明显，与共价键形成鲜明对比。

### 轨道贡献分析

通过 PDOS 可以确定：
- 价带顶的主要轨道成分：O 2p
- 导带底的主要轨道成分：Mg 3s
- 带隙跃迁类型：O 2p → Mg 3s

这些信息对理解材料的光学性质非常重要。
# 第五章：结果分析与常见问题

## 5.1 DOS/PDOS 图的解读

### 总 DOS 的分析

**判断材料类型**

通过费米能级处的 DOS 判断：
- DOS(E_F) > 0：金属
- DOS(E_F) = 0 且带隙 < 3 eV：半导体
- DOS(E_F) = 0 且带隙 > 3 eV：绝缘体

MgO 的带隙约 4-5 eV（PBE），属于绝缘体。

**带隙的确定**

从 DOS 图确定带隙：
1. 找到价带顶：费米能级以下 DOS 最后消失的位置
2. 找到导带底：费米能级以上 DOS 开始出现的位置
3. 带隙 = 导带底能量 - 价带顶能量

**注意：** PBE 泛函通常低估带隙，MgO 的实验带隙约 7.8 eV。

### PDOS 的分析

**轨道贡献识别**

通过 PDOS 确定：
- 价带主要由哪些轨道组成
- 导带主要由哪些轨道组成
- 不同元素的相对贡献

**成键性质判断**

| PDOS 特征 | 化学键类型 |
|----------|-----------|
| 不同元素 PDOS 分离 | 离子键 |
| 不同元素 PDOS 重叠 | 共价键 |
| 费米能级处 PDOS 连续 | 金属键 |

### MgO 案例的物理解释

**价带特征**
- O 2p 轨道主导，能量范围 -6 到 0 eV
- 对应 O²⁻ 的满壳层电子结构
- 6 个电子占据 3 个 p 轨道

**导带特征**
- Mg 3s 轨道主导，能量 > 4 eV
- 对应 Mg²⁺ 失去的 2 个电子的空轨道
- 激发态电子会回到 Mg 原子

**离子键形成**
- Mg 失去 2 个 3s 电子 → Mg²⁺
- O 获得 2 个电子填充 2p 轨道 → O²⁻
- PDOS 分离证实了完全的电子转移

## 5.2 常见问题排查

### Q1：PDOS 文件未生成

**症状：**
- `OUT.ABACUS/` 目录下没有 `PDOS` 文件
- 只有 `DOS1` 和 `TDOS` 文件

**可能原因：**
1. `basis_type` 不是 `lcao`
2. `out_dos` 设置为 1 而不是 2

**解决方法：**
```
# 检查 INPUT 文件
basis_type    lcao    # 必须是 lcao
out_dos       2       # 必须是 2
```

### Q2：DOS 曲线不平滑，有尖峰

**症状：**
- DOS 曲线呈锯齿状
- 有明显的尖峰

**可能原因：**
- k 点密度不够

**解决方法：**
1. 增加 k 点密度：
```
# KPT 文件
K_POINTS
0
Gamma
24 24 24 0 0 0    # 从 18×18×18 增加到 24×24×24
```

2. 或使用更小的 kspacing：
```
# INPUT 文件
kspacing    0.06    # 从 0.08 减小到 0.06
```

### Q3：费米能级提取错误

**症状：**
- grep 命令找不到 EFERMI
- 或者提取的值明显不合理

**可能原因：**
1. 查找的文件路径错误
2. 计算未正常完成

**解决方法：**
```bash
# 确认文件存在
ls OUT.ABACUS/running_nscf.log

# 使用正确的路径
grep -i 'efermi' OUT.ABACUS/running_nscf.log

# 或从 SCF 日志提取
grep -i 'efermi' OUT.ABACUS/running_scf.log
```

### Q4：NSCF 计算不收敛

**症状：**
- NSCF 计算报错或不收敛
- 提示找不到电荷密度文件

**可能原因：**
1. SCF 未完成或未输出电荷密度
2. `init_chg` 未设置为 `file`

**解决方法：**
1. 确认 SCF 已完成：
```bash
ls OUT.ABACUS/SPIN1_CHG
```

2. 检查 INPUT 参数：
```
calculation    nscf
init_chg       file    # 必须设置
```

## 5.3 进阶技巧

### 提高计算效率

**SCF 阶段：**
- 使用对称性：`symmetry = 1`
- 使用较稀疏的 k 点（12×12×12）
- 选择高效的求解器：`ks_solver = genelpa`

**NSCF 阶段：**
- 只计算需要的能量范围
- 使用并行计算

### 输出特定能量范围的 DOS

虽然 ABACUS 会输出完整能量范围的 DOS，但可以在后处理时限制范围：

```bash
abacustest model dos-pdos -j . --range -10 15
```

### 对比不同材料的 DOS

计算多个材料后，可以将 DOS 数据导出，使用绘图软件对比：

```bash
# 材料1
abacustest model dos-pdos -j material1/ --suffix mat1

# 材料2
abacustest model dos-pdos -j material2/ --suffix mat2
```

这将生成 `DOS_mat1.dat` 和 `DOS_mat2.dat`，便于对比。
# 附录

## 参数速查表

### DOS 计算关键参数

| 参数 | SCF | NSCF | 说明 |
|------|-----|------|------|
| calculation | scf | nscf | 计算类型 |
| symmetry | 1 | 0 | 对称性设置 |
| init_chg | atomic | file | 电荷密度来源 |
| out_chg | 1 | - | 输出电荷密度 |
| out_dos | - | 1或2 | 输出DOS（2包含PDOS） |
| basis_type | lcao | lcao | PDOS需要LCAO |

### smearing 参数推荐

| 体系类型 | smearing_method | smearing_sigma (Ry) |
|---------|----------------|-------------------|
| 金属 | mp 或 fd | 0.002-0.005 |
| 半导体 | gauss | 0.005-0.010 |
| 绝缘体 | gauss | 0.010-0.015 |

### k 点密度推荐

| 精度要求 | k 点网格 | 适用场景 |
|---------|---------|---------|
| 快速测试 | 6×6×6 | 初步测试 |
| 一般计算 | 12×12×12 | 常规计算 |
| 精确计算 | 18×18×18 | 发表论文 |
| 高精度 | 24×24×24 | 特殊需求 |

## 常用命令汇总

### ABACUS 运行命令

```bash
# 串行运行
abacus

# 并行运行（8进程）
mpirun -np 8 abacus

# 后台运行并保存日志
nohup mpirun -np 8 abacus > abacus.log 2>&1 &
```

### 费米能级提取

```bash
# 从 NSCF 日志提取
grep -i 'efermi' OUT.ABACUS/running_nscf.log

# 从 SCF 日志提取
grep -i 'efermi' OUT.ABACUS/running_scf.log
```

### abacustest 命令

```bash
# 基本用法
abacustest model dos-pdos -j .

# 指定能量范围
abacustest model dos-pdos -j . --range -5 7

# 按元素投影
abacustest model dos-pdos -j . --plot-type species

# 按轨道投影
abacustest model dos-pdos -j . --plot-type shell

# 查看帮助
abacustest model dos-pdos -h
```

## 常见问题汇总

| 问题 | 原因 | 解决方法 |
|------|------|---------|
| PDOS 文件未生成 | basis_type 不是 lcao 或 out_dos ≠ 2 | 设置 basis_type=lcao, out_dos=2 |
| DOS 曲线有尖峰 | k 点密度不够 | 增加 k 点网格密度 |
| 找不到费米能级 | 文件路径错误 | 检查 OUT.ABACUS/ 路径 |
| NSCF 不收敛 | 未读取电荷密度 | 设置 init_chg=file |
| 计算时间过长 | k 点过密或未用对称性 | SCF 用 symmetry=1 |

## 参考资料

### ABACUS 官方文档

- ABACUS 官网：http://abacus.ustc.edu.cn/
- 在线文档：https://abacus.deepmodeling.com/
- GitHub 仓库：https://github.com/deepmodeling/abacus-develop

### DOS/PDOS 相关文档

- DOS 计算教程：https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/dos.html
- PDOS 文件格式：https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/dos.html#pdos

### abacustest 工具

- GitHub：https://github.com/pxlxingliang/abacus-test
- Bohrium App：https://bohrium.dp.tech/apps/abacustest

## 进阶学习方向

完成本教程后，可以继续学习：

1. **杂化泛函计算**：使用 HSE06 获得更准确的带隙
2. **自旋极化 DOS**：计算磁性材料的自旋分辨 DOS
3. **投影能带**：结合能带和 PDOS 分析轨道贡献
4. **光学性质计算**：基于 DOS 计算介电函数和吸收谱
