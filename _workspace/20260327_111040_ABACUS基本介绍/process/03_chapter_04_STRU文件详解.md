# 第四章：STRU 文件详解

STRU 文件包含晶体结构的全部几何信息，以及赝势和轨道文件的路径。类似于 VASP 的 POSCAR + POTCAR 合并体。

## 4.1 文件结构概览

STRU 文件由若干"块"组成，各块按固定顺序排列：

```
ATOMIC_SPECIES          # 元素种类与赝势
NUMERICAL_ORBITAL       # 轨道文件（LCAO 才需要）
LATTICE_CONSTANT        # 晶格常数（缩放因子）
LATTICE_VECTORS         # 晶格矢量
ATOMIC_POSITIONS        # 原子坐标
```

PW 基组不需要 `NUMERICAL_ORBITAL` 块，其余块相同。

## 4.2 ATOMIC_SPECIES 块

格式：

```
ATOMIC_SPECIES
元素名   原子质量   赝势文件名
```

示例（MgO，两种元素）：

```
ATOMIC_SPECIES
Mg   24.305   Mg_ONCV_PBE-1.0.upf
O    15.999   O_ONCV_PBE-1.0.upf
```

- **元素名**：与 `ATOMIC_POSITIONS` 块中一致
- **原子质量**：用于分子动力学，SCF/relax 计算中不影响结果
- **赝势文件名**：文件需放在 INPUT 中 `pseudo_dir` 指定的目录下

**赝势文件获取：** http://abacus.ustc.edu.cn/pseudo/list.htm
推荐使用 ONCV 模守恒赝势（PBE 泛函），命名规范为 `元素_ONCV_PBE-1.0.upf`。

## 4.3 NUMERICAL_ORBITAL 块（LCAO 专用）

使用 `basis_type lcao` 时，需要额外指定每个元素的轨道文件：

```
NUMERICAL_ORBITAL
Mg_gga_8au_100Ry_4s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb
```

- 每个元素对应一行，顺序与 `ATOMIC_SPECIES` 中一致
- 轨道文件名的含义：`元素_gga_截断半径au_截断能Ry_轨道配置.orb`
  - `8au`：轨道截断半径为 8 a.u.（原子单位）
  - `100Ry`：对应截断能 100 Ry
  - `4s2p1d`：轨道配置（s、p、d 轨道的数量）
- 轨道文件需与所用赝势匹配（同一套 PBE 泛函）
- 文件放在 INPUT 中 `orbital_dir` 指定的目录下

**PW 基组时删去此块即可，其余不变。**

## 4.4 LATTICE_CONSTANT 块

```
LATTICE_CONSTANT
1.8897261258369282
```

这是一个缩放因子。后续的 `LATTICE_VECTORS` 中的数值乘以此常数，得到实际的晶格矢量长度（单位：Bohr）。

常用的对应关系：

| LATTICE_CONSTANT 值 | LATTICE_VECTORS 单位 |
|---------------------|----------------------|
| `1.8897261258369282` | Å（埃）|
| `1.0` | Bohr（原子单位）|

最常见的做法是将 LATTICE_CONSTANT 设为 `1.8897261258369282`，然后 LATTICE_VECTORS 中直接填写以 **Å** 为单位的晶格矢量。

## 4.5 LATTICE_VECTORS 块

```
LATTICE_VECTORS
4.27957   0.00000   0.00000
0.00000   4.27957   0.00000
0.00000   0.00000   4.27957
```

三行，每行三个数，构成 3×3 矩阵，描述三个晶格基矢的方向和（相对）长度。

- 实际晶格矢量 = 矩阵数值 × LATTICE_CONSTANT（单位 Bohr）
- 支持非正交晶格（如六方、三斜）
- 若 LATTICE_CONSTANT 设为 `1.8897261258369282`，则直接填 Å 数值

六方晶格示例（MgB₂，a = 3.083 Å，c = 3.524 Å）：

```
LATTICE_VECTORS
3.083000   0.000000   0.000000
-1.541500  2.669500   0.000000
0.000000   0.000000   3.524000
```

## 4.6 ATOMIC_POSITIONS 块

```
ATOMIC_POSITIONS
Direct
元素名
初始磁矩
原子数
x  y  z  mx  my  mz
...
```

逐项说明：

**第1行：** 坐标系类型
- `Direct`：晶格坐标（分数坐标，0 到 1 之间）
- `Cartesian`：笛卡尔坐标，单位为 LATTICE_CONSTANT（Bohr 或 Å）

**元素名：** 与 ATOMIC_SPECIES 中一致。

**初始磁矩：** 该元素所有原子的初始磁矩值（Bohr-mag）。非磁性计算填 `0.0`。

**原子数：** 该元素的原子个数。

**坐标行：** 每个原子一行，格式为：
```
x  y  z  mx  my  mz
```
- `x y z`：坐标值
- `mx my mz`：各方向是否允许移动（`1` = 可移动，`0` = 固定），用于 relax 计算

## 4.7 完整示例

### 示例1：Si 晶体（PW 基组，Diamond 结构，8原子超胞）

```
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
5.43070   0.00000   0.00000
0.00000   5.43070   0.00000
0.00000   0.00000   5.43070

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

### 示例2：MgO 晶体（LCAO 基组，岩盐结构，8原子超胞）

```
ATOMIC_SPECIES
Mg   24.305   Mg_ONCV_PBE-1.0.upf
O    15.999   O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Mg_gga_8au_100Ry_4s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
4.27957   0.00000   0.00000
0.00000   4.27957   0.00000
0.00000   0.00000   4.27957

ATOMIC_POSITIONS
Direct
Mg
0.0
4
0.0   0.0   0.0   0  0  0
0.0   0.5   0.5   0  0  0
0.5   0.0   0.5   0  0  0
0.5   0.5   0.0   0  0  0
O
0.0
4
0.5   0.0   0.0   0  0  0
0.5   0.5   0.5   0  0  0
0.0   0.0   0.5   0  0  0
0.0   0.5   0.0   0  0  0
```

注意：MgO 示例中原子坐标后的 `0 0 0` 表示结构优化时固定所有原子（也可改为 `1 1 1` 允许移动）。

> **获取晶体结构：** 可从 Materials Project（https://materialsproject.org）、ICSD 等数据库下载 CIF 文件，再通过工具（如 ASE、abacustest）转换为 STRU 格式。下一章的案例中会展示 abacustest 的转换用法。
