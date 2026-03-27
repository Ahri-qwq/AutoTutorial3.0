## 第一章 STRU：晶体结构文件

`STRU` 文件存储晶体的结构信息，包括元素种类、赝势/轨道文件名、晶格参数和原子坐标。ABACUS 中必须将这个文件命名为 `STRU`，不可更改。

### 1.1 文件整体结构

STRU 文件由若干以关键字开头的块组成，按顺序排列。对于平面波（PW）计算，必选块有四个；LCAO 计算需额外添加 `NUMERICAL_ORBITAL` 块：

```
# 注释行（可以用 # 开头）

ATOMIC_SPECIES          # 元素种类和赝势文件（必选）

NUMERICAL_ORBITAL       # 轨道文件（LCAO 专用，PW 不需要）

LATTICE_CONSTANT        # 晶格常数（必选）

LATTICE_VECTORS         # 晶格向量（必选）

ATOMIC_POSITIONS        # 原子坐标（必选）
```

各块的顺序不能颠倒。

### 1.2 各块详解

#### ATOMIC_SPECIES

```
ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf
```

每行格式：`元素名  原子量  赝势文件名`

- 元素名：需与 ATOMIC_POSITIONS 中的标签一致
- 原子量：用于分子动力学中的质量，SCF 计算结果与此值无关
- 赝势文件名：ABACUS 使用模守恒（ONCV）赝势，UPF 格式。文件放在 INPUT 中 `pseudo_dir` 指定的目录下。推荐使用 [SG15 ONCV 赝势库](http://abacus.ustc.edu.cn/pseudo/list.htm)

多元素体系每行写一个：

```
ATOMIC_SPECIES
Mg  24.305  Mg_ONCV_PBE-1.0.upf
O   15.999  O_ONCV_PBE-1.0.upf
```

#### NUMERICAL_ORBITAL（LCAO 专用）

```
NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb
```

每行一个轨道文件，顺序与 ATOMIC_SPECIES 中的元素对应。文件放在 INPUT 中 `orbital_dir` 指定目录下。PW 计算删去这整个块。

轨道文件命名格式：`元素_泛函_截断半径_能量截断_基组规模.orb`
- 例：`Al_gga_7au_100Ry_4s4p1d.orb` 表示 GGA 泛函、截断半径 7 Bohr、100 Ry 能量截断、4s4p1d 基函数

#### LATTICE_CONSTANT

```
LATTICE_CONSTANT
1.8897259886
```

晶格常数的单位因子（单位：Bohr）。LATTICE_VECTORS 中的数值乘以这个因子才是实际长度。

常用换算：
- `1.8897259886` Bohr = 1.0 Å（最常用写法：把长度写成 Å，LATTICE_CONSTANT 写 1.8897259886）
- 直接写 `1.0`：则 LATTICE_VECTORS 中的单位为 Bohr

#### LATTICE_VECTORS

```
LATTICE_VECTORS
4.03460  0.00000  0.00000
0.00000  4.03460  0.00000
0.00000  0.00000  4.03460
```

三行分别是晶格向量 **a₁、a₂、a₃**，单位为 LATTICE_CONSTANT。对于立方晶系，三个向量互相垂直且等长。

实际晶格常数 = LATTICE_VECTORS 中的值 × LATTICE_CONSTANT。上例中 Al FCC 晶格常数 = 4.03460 × 1.8897259886 Bohr = 4.03460 Å。

#### ATOMIC_POSITIONS

```
ATOMIC_POSITIONS
Direct            # 坐标类型

Al                # 元素名（与 ATOMIC_SPECIES 对应）
0                 # 初始磁矩（Bohr mag，非磁体系写 0）
4                 # 原子数目
0.0  0.0  0.0  0 0 0   # x y z  move_x move_y move_z
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

**坐标类型**（第一行）：
- `Direct`：分数坐标（推荐）。坐标值在 0~1 之间，以晶格向量为基
- `Cartesian`：直角坐标，单位为 LATTICE_CONSTANT
- `Cartesian_angstrom`：直角坐标，单位为 Å

**move_x move_y move_z**：控制原子在结构优化时是否可以移动（1=可移动，0=固定）。SCF 计算中这三个数没有实际作用，习惯写 `0 0 0` 或 `1 1 1` 均可。

多元素体系，各元素块依次排列（元素名→磁矩→原子数→坐标行）。

### 1.3 Al FCC 完整示例：PW 写法

铝（Al）面心立方（FCC）结构，原胞含 4 个原子，实验晶格常数 4.050 Å，计算优化值约 4.035 Å。

```
# Al FCC，PW 计算，4 原子原胞

ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897259886

LATTICE_VECTORS
4.03460  0.00000  0.00000
0.00000  4.03460  0.00000
0.00000  0.00000  4.03460

ATOMIC_POSITIONS
Direct
Al
0
4
0.0  0.0  0.0  0 0 0
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

注意：PW 写法不包含 `NUMERICAL_ORBITAL` 块。

### 1.4 Al FCC 完整示例：LCAO 写法

```
# Al FCC，LCAO 计算，4 原子原胞

ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.8897259886

LATTICE_VECTORS
4.03460  0.00000  0.00000
0.00000  4.03460  0.00000
0.00000  0.00000  4.03460

ATOMIC_POSITIONS
Direct
Al
0
4
0.0  0.0  0.0  0 0 0
0.5  0.5  0.0  0 0 0
0.5  0.0  0.5  0 0 0
0.0  0.5  0.5  0 0 0
```

与 PW 写法唯一的区别是添加了 `NUMERICAL_ORBITAL` 块。同时，INPUT 中需要将 `basis_type` 改为 `lcao`，并填写 `orbital_dir` 路径。

### 1.5 PW 与 LCAO 写法对比

| 项目 | PW | LCAO |
|------|----|----|
| NUMERICAL_ORBITAL 块 | 不需要 | 必须 |
| INPUT 中 basis_type | pw | lcao |
| INPUT 中 orbital_dir | 不需要 | 必须 |
| ecutwfc 意义 | 控制基组大小 | 控制辅助平面波格点密度 |
| 适用场景 | 收敛性测试、高精度 | 大体系、O(N) 线性标度 |
