## 三、输入文件准备

### 3.1 STRU 文件：反铁磁结构的关键设置

本案例计算 II 型反铁磁 NiO，原胞中两个 Ni 原子的磁矩大小相等、方向相反。

在 ABACUS 的 STRU 文件中，处理反铁磁的常用做法是**将不同初始磁矩的同种元素定义为不同的 atomic species**。虽然它们使用相同的赝势和轨道文件，但通过分别定义 Ni1 和 Ni2，可以对每种 species 统一设置磁矩方向：

```
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf
Ni2 58.693 Ni_ONCV_PBE-1.0.upf
O   15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Ni_gga_9au_100Ry_4s2p2d1f.orb
Ni_gga_9au_100Ry_4s2p2d1f.orb
O_gga_7au_100Ry_2s2p1d.orb
```

在 `ATOMIC_POSITIONS` 部分，分别为 Ni1 和 Ni2 指定初始磁矩：

```
ATOMIC_POSITIONS
Cartesian_angstrom

Ni1
2.0     // 初始磁矩 +2.0 μB
...     // 原子坐标

Ni2
-2.0    // 初始磁矩 -2.0 μB
...     // 原子坐标

O
0.0
...
```

这里 `Ni1` 和 `Ni2` 共用同一套赝势（`Ni_ONCV_PBE-1.0.upf`）和轨道文件（`Ni_gga_9au_100Ry_4s2p2d1f.orb`），差别仅在于初始磁矩符号相反。

> **说明：** 使用双 species 方案的好处是，DFT+U 的参数（`orbital_corr`、`hubbard_u`）按 species 顺序指定，可以独立控制每种 Ni 的 U 值。两个 Ni species 使用相同的 U 值是物理上正确的——它们本质上是同一种元素。

---

### 3.2 INPUT 文件：DFT+U 参数逐一解析

完整的 INPUT 文件如下（仅列出与 DFT+U 相关的关键部分）：

```
INPUT_PARAMETERS
suffix              NiO
calculation         scf
basis_type          lcao
nspin               2

#Parameter DFT+U
dft_plus_u          1
orbital_corr        2 2 -1
hubbard_u           5.0 5.0 0.0

#输出控制
out_bandgap         1
out_mul             1
out_chg             1
```

**核心三参数：**

| 参数 | 含义 | 本例设置 |
|------|------|----------|
| `dft_plus_u` | DFT+U 总开关，1 为开启，0 为关闭（即使设置了其他+U参数也不会生效） | 1 |
| `orbital_corr` | 数组，长度为 ntype（atomic species 数），指定每种 species 施加 +U 的 l 量子数；-1 表示不施加 | 2 2 -1 |
| `hubbard_u` | 数组，长度与 orbital_corr 相同，单位 eV，指定每种 species 的 U 值 | 5.0 5.0 0.0 |

本例中 `ntype = 3`（Ni1、Ni2、O），`orbital_corr = 2 2 -1` 表示：
- Ni1：在 l=2（d 轨道）施加 +U，U = 5.0 eV
- Ni2：在 l=2（d 轨道）施加 +U，U = 5.0 eV
- O：不施加 +U（-1）

**输出控制参数：**

| 参数 | 作用 |
|------|------|
| `out_bandgap 1` | 在日志中输出每个自旋通道的能隙 |
| `out_mul 1` | 输出 Mulliken 布居分析（mulliken.txt），含原子磁矩 |
| `out_chg 1` | 输出电荷密度及末步的 onsite.dm 文件 |

---

### 3.3 KPT 文件

```
K_POINTS
0
Gamma
4 4 4 0 0 0
```

4×4×4 的 Gamma 中心 k 点网格，适用于 NiO 的 SCF 计算。

---

**文件树总览：**

```
ABACUS_DFT+U/
├── INPUT
├── STRU
├── KPT
├── Ni_ONCV_PBE-1.0.upf
├── Ni_gga_9au_100Ry_4s2p2d1f.orb
├── O_ONCV_PBE-1.0.upf
└── O_gga_7au_100Ry_2s2p1d.orb
```
