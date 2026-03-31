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
