# 第四章：Step 1 — 计算二阶力常数

进入 `2nd` 目录。本步骤使用 ABACUS + Phonopy 计算二阶力常数，最终得到 `FORCE_CONSTANTS_2ND`。

## 4.1 结构优化

计算力常数前，先对体系进行结构优化（`calculation = relax`），确保原子处于受力平衡的平衡构型。
本例使用 2×2×2 k 点采样，ecut = 100 Ry（注意：实际研究中应测试 k 点收敛性）。

结构优化完成后，得到以下 STRU 文件：

```
ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
0 2.81594778072 2.81594778072 #latvec1
2.81594778072 0 2.81594778072 #latvec2
2.81594778072 2.81594778072 0 #latvec3

ATOMIC_POSITIONS
Direct # direct coordinate

Si #label
0 #magnetism
2 #number of atoms
0.875  0.875  0.875  m  0  0  0
0.125  0.125  0.125  m  0  0  0
```

说明：
- `LATTICE_CONSTANT = 1.88972612546` 单位为 Bohr（等于 1 Å），是格矢量的缩放因子
- 格矢量乘以此缩放因子，得到 FCC Si 的实际格矢（约 5.31 Bohr ≈ 2.82 Å）
- `m 0 0 0` 表示原子坐标固定，为 relax 优化后的最终平衡位置
- Si 的质量（28.0855）在力常数计算中不起作用

## 4.2 Phonopy 产生超胞微扰构型

在 `2nd` 目录下执行：

```bash
phonopy setting.conf --abacus -d
```

其中 `setting.conf` 的内容为：

```
DIM = 2 2 2
ATOM_NAME = Si
```

`DIM = 2 2 2` 表示将 2 原子原胞扩展为 2×2×2 超胞（共 16 个原子）。
Phonopy 会生成需要计算受力的微扰超胞构型文件。

Si 金刚石结构对称性高，只需产生 **1 个**微扰构型 `STRU-001`。

## 4.3 ABACUS SCF 计算受力

对 `STRU-001` 进行单点 SCF 计算，获取原子受力。INPUT 中的关键设置：

```
calculation   scf
cal_force     1        # 必须开启受力输出
stru_file     STRU-001 # 指向 Phonopy 生成的超胞构型
```

**小技巧**：通过 `stru_file` 参数直接指定 Phonopy 生成的构型文件，无需手动重命名。

计算完成后，用以下命令从 ABACUS 输出中提取受力，生成 `FORCE_SETS` 文件：

```bash
phonopy -f OUT.DIA-50/running_scf.log
```

其中 `OUT.DIA-50` 是 ABACUS 的输出目录（目录名由 INPUT 中的 `suffix` 参数决定，
若 `suffix = DIA-50` 则输出目录为 `OUT.DIA-50`）。

## 4.4 计算声子谱和二阶力常数

执行：

```bash
phonopy -p band.conf --abacus
```

`band.conf` 的完整内容为：

```
ATOM_NAME = Si
DIM = 2 2 2
MESH = 8 8 8
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.5  0.625 0.25 0.625, 0.375 0.375 0.75  0.0 0.0 0.0  0.5 0.5 0.5
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
FORCE_CONSTANTS = WRITE
FULL_FORCE_CONSTANTS = .TRUE.
```

参数说明：
- `MESH = 8 8 8`：声子谱积分的 q 点网格
- `PRIMITIVE_AXES`：原始布里渊区的基矢，这里 Si FCC 取标准设置
- `BAND`：高对称路径点，覆盖 Γ-X-K-L-Γ-L 路径
- `FORCE_CONSTANTS = WRITE`：将力常数写入文件
- **`FULL_FORCE_CONSTANTS = .TRUE.`**：输出完整力常数矩阵（必须设置，见注意事项）

命令执行后，Phonopy 输出：
- `band.yaml`：声子谱数据（可用于绘图）
- `FORCE_CONSTANTS`：二阶力常数文件

如需导出 gnuplot 格式的声子谱：

```bash
phonopy-bandplot --gnuplot > pho.dat
```

## 4.5 单位转换：生成 FORCE_CONSTANTS_2ND

ShengBTE 要求 `FORCE_CONSTANTS_2ND` 的单位为 **eV/Å²**，
但 ABACUS 结合 Phonopy 输出的 `FORCE_CONSTANTS` 单位为 **eV/(Å·au)**，其中 1 au = 0.52918 Å。

运行单位转换脚本：

```bash
python au2si.py
```

脚本将 `FORCE_CONSTANTS` 乘以相应的单位换算系数，输出 `FORCE_CONSTANTS_2ND`，即 ShengBTE 所需的格式。

至此，第一步完成。`shengbte/` 目录中已提供参考用的 `FORCE_CONSTANTS_2ND` 文件，可用于对照验证。
