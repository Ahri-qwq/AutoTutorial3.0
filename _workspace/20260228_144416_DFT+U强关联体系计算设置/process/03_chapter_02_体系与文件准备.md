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
Ni_gga_8au_100Ry_4s2p2d1f.orb     # Ni1的轨道文件
Ni_gga_8au_100Ry_4s2p2d1f.orb     # Ni2共用同一轨道文件
O_gga_6au_100Ry_2s2p1d.orb
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
smearing_method         gauss
smearing_sigma          0.01

#Parameters (Mixing)
mixing_type             broyden
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
| `Ni_gga_8au_100Ry_4s2p2d1f.orb` | Ni 的数值原子轨道（截断半径8 au，能量截断100 Ry） |
| `O_gga_6au_100Ry_2s2p1d.orb` | O 的数值原子轨道 |

轨道文件中，Ni 只有一个 d 基组（zeta=0），这在后续分析 occupation matrix 时会用到。
