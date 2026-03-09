# 四、案例二：BCC 铁 ELF（自旋极化）

算例地址：https://github.com/MCresearch/abacus-user-guide/tree/master/examples/elf/bcc-Fe-pw

体心立方铁（BCC Fe）是铁磁性金属，使用 PW 基组计算自旋极化的 ELF，可以得到自旋向上和自旋向下电子各自的局域化图像。

## 4.1 输入文件

**INPUT：**

```
INPUT_PARAMETERS
suffix              autotest
calculation         scf
ntype               1
basis_type          pw
pseudo_dir          ./

ecutwfc             100          # 截断能，单位 Ry
scf_thr             1e-6
scf_nmax            200

nspin               2            # 自旋极化（共线磁矩）
smearing_method     gauss
smearing_sigma      0.01

mixing_type         broyden
mixing_beta         0.7

out_elf             1 3          # 输出 ELF（总/自旋↑/自旋↓）
```

**STRU：**

BCC Fe 传统晶胞含 2 个原子（格点参数 a = 2.866 Å），初始磁矩按元素设置为 2.3 μB：

```
ATOMIC_SPECIES
Fe  55.845  Fe_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.889726         # Bohr/Å

LATTICE_VECTORS
2.866  0.000  0.000
0.000  2.866  0.000
0.000  0.000  2.866

ATOMIC_POSITIONS
Direct

Fe
2.3              # 初始磁矩，单位 μB
2
0.000  0.000  0.000  m  1  1  1
0.500  0.500  0.500  m  1  1  1
```

初始磁矩设置为 2.3 μB（接近铁的实验磁矩 ~2.2 μB）。设置合理的初始磁矩对收敛到铁磁基态至关重要——过小的初始磁矩可能导致 SCF 收敛到非磁态。

**KPT：**

```
K_POINTS
0
Gamma
8 8 8 0 0 0
```

## 4.2 运行与输出文件

```bash
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

`nspin = 2` 时，ABACUS 输出三个 cube 文件：

| 文件名 | 内容 |
|--------|------|
| `ELF.cube` | 总 ELF（自旋向上+向下的综合） |
| `ELF_SPIN1.cube` | 自旋向上（↑）电子的 ELF |
| `ELF_SPIN2.cube` | 自旋向下（↓）电子的 ELF |

## 4.3 VESTA 可视化

对每个 cube 文件分别在 VESTA 中打开，查看 (100) 晶面截面图：

1. File → Open → 选择 `ELF.cube`（或 `ELF_SPIN1.cube` / `ELF_SPIN2.cube`）
2. Properties → Sections
3. 设置截面方向：选择 (100) 平面（设置法向量为 [1, 0, 0]）
4. 调整颜色映射（推荐蓝-白-红方案，蓝=0，红=1）

## 4.4 结果解读

对三种 ELF 的 (100) 面截面图：

**总 ELF（ELF.cube）：**
- Fe 原子核附近 ELF 接近 1（强烈局域的芯电子/内层电子）
- 原子间区域 ELF 较低（金属 d 电子相对离域）
- 整体图像具有 BCC 晶体的四重对称性

**自旋向上 ELF（ELF_SPIN1.cube）：**
- 自旋向上电子在 Fe 原子附近的局域化分布
- 反映了 d 轨道中自旋向上电子的空间分布特征

**自旋向下 ELF（ELF_SPIN2.cube）：**
- 与自旋向上 ELF 形状相似，但数值略有不同
- 自旋向上和向下的 ELF 差异体现了铁磁性——自旋向上电子更多（磁矩约 2.2 μB），其 ELF 峰值略高

BCC Fe 与水分子的对比：金属体系的 ELF 在原子间区域通常在 0.3–0.5 之间（类均匀电子气），而分子的共价键区域 ELF 可达 0.7 以上。这种对比直观地说明了金属键与共价键在电子局域化程度上的本质差异。
