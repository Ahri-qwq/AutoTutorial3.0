# 三、案例一：水分子 ELF（自旋非极化）

本案例对真空中的水分子（H₂O）计算 ELF，分别使用平面波（PW）和原子轨道（LCAO）两种基组。

算例文件：
- PW：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/H20-pw
- LCAO：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/H20-lcao

## 3.1 平面波计算

### 输入文件

**INPUT：**

```
INPUT_PARAMETERS
suffix              autotest
calculation         scf
ntype               2
basis_type          pw
pseudo_dir          ./

ecutwfc             100          # 截断能，单位 Ry
scf_thr             1e-6
scf_nmax            100

nspin               1            # 自旋非极化
symmetry            0

out_elf             1 3          # 输出 ELF，保留 3 位有效数字
```

**STRU：**

水分子置于约 14.82 Å 的真空立方盒内，分子坐标来自 ELF.cube 头部数据（单位 Bohr）：

```
ATOMIC_SPECIES
H   1.008    H_ONCV_PBE-1.0.upf
O  15.999    O_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.889726         # Bohr/Å

LATTICE_VECTORS
14.82   0.00    0.00
 0.00  14.82    0.00
 0.00   0.00   14.82

ATOMIC_POSITIONS
Direct

H
0.0
2
0.4274  0.6688  0.2972  0  0  0    # H1
0.5087  0.7379  0.2679  0  0  0    # H2

O
0.0
1
0.4841  0.7029  0.3264  0  0  0    # O
```

**KPT：**

水分子在真空盒内，只需 Γ 点采样：

```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

### 运行计算

```bash
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

### 输出文件

计算结束后，在 `OUT.autotest/` 目录下生成 `ELF.cube`。文件头部如下：

```
STEP: 0  Cubefile created from ABACUS. Inner loop is z, followed by y and x
1 (nspin)
3 0.0 0.0 0.0
180 0.155556 0.000000 0.000000
180 0.000000 0.155556 0.000000
180 0.000000 0.000000 0.155556
 1 1.000000 11.967714 18.726812 8.322864
 1 1.000000 14.244302 20.657387 7.503683
 8 6.000000 13.550429 19.681072 9.139011
 0.000e+00 0.000e+00 0.000e+00 0.000e+00 ...
```

- 180×180×180 FFT 格点，步长 0.155556 Bohr（约 0.082 Å）
- 对应约 28 Bohr（~14.8 Å）的真空盒子
- 三行原子信息：H（序数 1）、H（序数 1）、O（序数 8），坐标单位 Bohr
- 远离分子的真空区域 ELF≈0

### VESTA 可视化

**等高面图（3D）：**

1. 用 VESTA 打开 `ELF.cube`（File → Open）
2. 菜单 Properties → Isosurfaces
3. 添加等值面，推荐 ELF = 0.75（捕捉成键区域）
4. 调整颜色和透明度

**截面图（2D）：**

1. Properties → Sections/Isosurfaces
2. 选择 Section
3. 勾选 "Specify from three atoms"，选取 H1、O、H2 三个原子，得到 H-O-H 平面截面

**可视化结果解读：**

| 区域 | ELF 值 | 物理含义 |
|------|--------|---------|
| O-H 键区（两核之间） | ~0.7–0.8 | 共价键电子，中等局域化 |
| O 原子两侧（孤对电子） | 接近 1 | 高度局域的孤对电子 |
| 远离分子的真空区域 | ≈ 0 | 无电子分布 |

水分子是典型的极性共价分子，O 原子有两个孤对电子，ELF 图像可以直观地显示出这两个孤对电子的"兔耳"形分布。

## 3.2 原子轨道计算（LCAO）

算例地址：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/H20-lcao

### 与 PW 的差异

LCAO 计算只需在 INPUT 和 STRU 中做少量修改，STRU 中的坐标与 PW 相同：

**INPUT（LCAO 版本）：**

```
INPUT_PARAMETERS
suffix              autotest
calculation         scf
ntype               2
basis_type          lcao          # 改为 LCAO
pseudo_dir          ./
orbital_dir         ./            # 新增：指定轨道文件目录

ecutwfc             100
scf_thr             1e-6
scf_nmax            100

nspin               1
symmetry            0

out_elf             1 3
```

**STRU（LCAO 版本，新增 NUMERICAL_ORBITAL 段）：**

```
ATOMIC_SPECIES
H   1.008    H_ONCV_PBE-1.0.upf
O  15.999    O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_6au_100Ry_2s1p.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
14.82   0.00    0.00
 0.00  14.82    0.00
 0.00   0.00   14.82

ATOMIC_POSITIONS
Direct

H
0.0
2
0.4274  0.6688  0.2972  0  0  0
0.5087  0.7379  0.2679  0  0  0

O
0.0
1
0.4841  0.7029  0.3264  0  0  0
```

轨道文件命名含义：以 `O_gga_7au_100Ry_2s2p1d.orb` 为例：
- `gga`：GGA 泛函
- `7au`：截断半径 7 Bohr
- `100Ry`：与 ecutwfc = 100 Ry 匹配
- `2s2p1d`：2 个 s 轨道，2 个 p 轨道，1 个 d 轨道

### LCAO 稳定性修正

LCAO 基组存在不完备性，在远离原子核的区域，动能密度的数值行为可能出现异常，导致真空区出现非物理的非零 ELF。ABACUS 采用文献 [2] 的建议，在动能密度的分子部分加入 10⁻⁵ 的小量修正。该修正确保了远场行为的正确性，同时不影响近核区域的结果。PW 基组没有此问题。

### PW vs LCAO 结果对比

LCAO 计算得到的 ELF 与 PW 结果大体一致，O-H 键区和孤对电子均有清晰的高 ELF 峰，但由于基组不同，细节上存在差异：

- 化学键区（近核区域）：两种基组结果接近
- 真空区域：PW 结果干净（ELF≈0）；LCAO 经修正后也正常
- 孤对电子形状：LCAO 略有差异，依赖基组大小（DZP vs TZP）

对于精确的 ELF 分析，PW 基组更可靠；LCAO 计算速度更快，适合大体系的定性分析。
