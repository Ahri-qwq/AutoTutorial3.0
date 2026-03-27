# 第二章：案例一——硅的能带结构

## 二、案例一：硅（Si）的能带结构

### 2.1 体系说明

硅是金刚石立方结构（空间群 Fd$\bar{3}$m），晶格常数 5.4307 Å，每个原胞含 2 个 Si 原子。硅是典型的间接带隙半导体：价带顶位于布里渊区中心（Γ 点），导带底位于 Γ→X 方向上约 0.85 的位置（Δ 点）。

PBE 泛函计算得到的带隙约为 **0.57 eV**，低于实验值 1.17 eV——这是 PBE 的系统性低估，属于正常现象，不是计算错误。

本案例使用 LCAO 基组，结构为 2 原子原胞（primitive cell）。

---

### 2.2 输入文件准备

#### STRU 文件

以下是 Si 原胞的结构文件（2 个原子，FCC 原胞向量）：

```
# STRU
ATOMIC_SPECIES
Si  28.085  Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
0.000000000  2.715350000  2.715350000  #latvec1
2.715350000  0.000000000  2.715350000  #latvec2
2.715350000  2.715350000  0.000000000  #latvec3

ATOMIC_POSITIONS
Direct

Si  #label
0   #magnetism
2   #number of atoms
0.00  0.00  0.00  m  1  1  1
0.25  0.25  0.25  m  1  1  1
```

赝势和轨道文件说明：
- `Si_ONCV_PBE-1.0.upf`：模守恒赝势（ONCV），PBE 泛函
- `Si_gga_8au_60Ry_2s2p1d.orb`：GGA 泛函，截断半径 8 a.u.，DZP 基组

> 赝势和轨道文件可从 ABACUS 的 [PP_ORB 目录](https://github.com/deepmodeling/abacus-develop/tree/develop/tests/PP_ORB)下载，或直接使用计算环境中 `/abacus-develop/tests/PP_ORB/` 路径。

#### INPUT 文件（SCF）

```
# INPUT（SCF，自洽计算）
INPUT_PARAMETERS
suffix          Si
pseudo_dir      ./
orbital_dir     ./

calculation     scf
ecutwfc         50
scf_thr         1e-7
scf_nmax        300

basis_type      lcao

smearing_method gauss
smearing_sigma  0.01

mixing_type     broyden
mixing_beta     0.4

symmetry        1
out_chg         1        # 输出电荷密度，供 NSCF 步骤读取
```

关键参数说明：

- `out_chg 1`：在 `OUT.Si/` 目录下输出 `SPIN1_CHG.cube`，NSCF 步骤必须用到
- `ecutwfc 50`：平面波截断能（Ry），LCAO 计算中该参数控制电子密度网格精度
- `smearing_sigma 0.01`：Si 是半导体，展宽可取较小值（0.01 Ry ≈ 0.136 eV）

#### KPT 文件（SCF）

SCF 阶段使用均匀 k 网格：

```
# KPT（SCF）
K_POINTS
0
Gamma
9 9 9 0 0 0
```

`9×9×9` 的 Gamma 中心网格，对 Si 原胞（2原子）能够给出收敛的电子密度。

#### INPUT 文件（NSCF 能带）

```
# INPUT（NSCF，非自洽能带计算）
INPUT_PARAMETERS
suffix          Si
pseudo_dir      ./
orbital_dir     ./

init_chg        file     # 读取 SCF 输出的电荷密度
calculation     nscf
ecutwfc         50
scf_thr         1e-7
scf_nmax        300

basis_type      lcao

symmetry        0        # NSCF 必须关闭对称性
out_band        1        # 输出能带文件 BANDS_1.dat
out_bandgap     1        # 在 running_nscf.log 中输出带隙信息
```

与 SCF INPUT 的差异：
- `calculation` 改为 `nscf`
- `init_chg file`：从 `OUT.Si/` 读取 `SPIN1_CHG.cube`
- `symmetry 0`：关闭对称性（原因见第一章）
- 添加 `out_band 1` 和 `out_bandgap 1`
- 去掉 `out_chg 1`（NSCF 不需要再输出电荷密度）

#### KPT 文件（NSCF 能带）

Si 原胞属于 FCC 布里渊区，高对称路径选取 Γ-X-U|K-Γ-L-W-X：

```
# KPT（NSCF 能带，FCC 布里渊区高对称路径）
K_POINTS
8
Line
0.000  0.000  0.000  30  # Gamma
0.500  0.000  0.500  20  # X
0.625  0.250  0.625  20  # U
0.375  0.375  0.750  50  # K
0.000  0.000  0.000  50  # Gamma
0.500  0.500  0.500  20  # L
0.500  0.250  0.750  20  # W
0.500  0.000  0.500   1  # X
```

格式说明：
- 第一行 `K_POINTS` 为关键字
- 第二行 `8`：高对称点数目
- 第三行 `Line`：能带计算专用模式
- 从第四行起：每行为一个高对称点的分数坐标 + 到下一点的插值 k 点数 + 注释
- **最后一行** 的插值点数填 `1`，表示路径在此结束

> `U` 和 `K` 是同一点在不同路径段的端点（在 FCC 布里渊区中两者等价，坐标略不同），中间无插值，因此 U 行的点数填 `20`（到下一段 K 的距离），而不是 0。

---

### 2.3 运行计算

**目录结构建议：**

```
Si_band/
├── STRU
├── INPUT           # 运行时根据步骤切换
├── KPT             # 运行时根据步骤切换
├── INPUT_scf
├── INPUT_nscf
├── KPT_scf
├── KPT_nscf
├── Si_ONCV_PBE-1.0.upf
└── Si_gga_8au_60Ry_2s2p1d.orb
```

**Step 1：SCF 计算**

```bash
cp INPUT_scf INPUT
cp KPT_scf KPT
mpirun -np 4 abacus | tee scf.log
```

计算完成后，检查收敛：

```bash
grep "convergence" OUT.Si/running_scf.log
# 应出现：convergence is achieved
```

**Step 2：NSCF 能带计算**

```bash
cp INPUT_nscf INPUT
cp KPT_nscf KPT
mpirun -np 4 abacus | tee nscf.log
```

NSCF 通常比 SCF 快得多（不迭代电荷密度），运行完成后在 `OUT.Si/` 目录下会生成：

- `BANDS_1.dat`：能带数据
- `running_nscf.log`：含带隙信息

---

### 2.4 费米能级提取

绘制能带图需要把能量零点移到费米能级。从日志文件提取：

```bash
grep "EFERMI" OUT.Si/running_nscf.log
# 输出示例：EFERMI = 6.389520000 eV
```

也可以从 SCF 日志中提取（两者应一致）：

```bash
grep "EFERMI" OUT.Si/running_scf.log
```

记录该值，后续配置 `config.json` 时使用。

---

### 2.5 Atomkit 后处理——绘制能带图

进入输出目录，将 KPT 文件复制进去（`abacus-plot` 需要 KPT 来确定高对称点位置）：

```bash
cp KPT_nscf OUT.Si/KPT
cd OUT.Si/
```

创建 `config.json`：

```json
{
    "bandfile": "BANDS_1.dat",
    "efermi": 6.389520000,
    "energy_range": [-6, 6],
    "bandfig": "band.png",
    "kptfile": "KPT",
    "dpi": 300
}
```

参数说明：
- `efermi`：从上一步提取的费米能（eV），能带图以此为零点
- `energy_range`：图的纵坐标范围（eV），`[-6, 6]` 通常足够展示主要特征
- `kptfile`：高对称点标注所需的 KPT 文件路径

运行绘图命令：

```bash
abacus-plot -b
```

生成 `band.png`。

---

### 2.6 结果解读

从能带图和日志文件可以提取关键信息：

**带隙查询：**

```bash
grep -i "bandgap" OUT.Si/running_nscf.log
# 示例输出：E_bandgap  0.5700  eV
```

PBE 计算的 Si 带隙约为 **0.57 eV**，间接带隙特征表现为：
- 价带顶（VBM）：位于 Γ 点
- 导带底（CBM）：位于 Γ→X 路径上（Δ 点，约 0.85 × Γ→X 处）

这与实验值 1.17 eV 的差距是 PBE 的系统性低估，并非计算错误。若需要更准确的带隙，可使用 HSE06 杂化泛函（见第四章进阶提示）。

**能带特征：**
- Γ 点附近：价带的最高三条带简并（重空穴带、轻空穴带、自旋轨道劈裂带），PBE 下无 SOC 时三重简并
- L 点附近：导带存在明显极小值
- 整体带宽约 12 eV（sp 价带区）
