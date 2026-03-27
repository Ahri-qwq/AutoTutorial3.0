
# （修改后完整稿件）

---

ABACUS 的电子结构计算包含两类核心任务：电子能带结构（band structure）和态密度（density of states，DOS）。两类计算共享同一套流程骨架：先做自洽计算（SCF）获得收敛的电子密度，再做非自洽计算（NSCF）读入该密度，在任意 k 点上求解 Kohn-Sham 方程。本教程通过两个典型案例演示完整流程：

- **案例一：硅（Si）** — 计算能带结构，验证 PBE 间接带隙约 0.57 eV
- **案例二：铝（Al）** — 计算 TDOS 和 PDOS，分析金属特征电子结构

后处理统一使用 Atomkit 的 `abacus-plot` 命令完成绘图。

**适用读者：** 已安装 ABACUS 并能运行基础 SCF 计算的用户

**前置知识：** ABACUS 的 INPUT/STRU/KPT 文件格式，DFT 基本概念

**计算环境：** 案例一使用 LCAO 基组，案例二使用平面波（PW）基组，均基于 ABACUS v3.10.1

---

## 一、DOS 与能带结构

电子能带结构描述材料中电子的色散关系——能量 $E$ 随波矢 **k** 的变化。通过能带图可以判断：

- **带隙**：价带顶与导带底之间的能量差，决定材料是金属（无带隙）、半导体还是绝缘体
- **间接/直接带隙**：价带顶和导带底是否位于同一 k 点，直接影响光学跃迁效率
- **有效质量**：能带曲率，影响载流子迁移率

态密度 $g(E)$ 是能量 $E$ 附近单位能量区间内可供电子占据的状态数，由能带结构对整个布里渊区积分得到。投影态密度（PDOS）进一步将 $g(E)$ 分解到各原子的各轨道角动量通道（s、p、d……），用于分析成键和磁性来源。

**两步计算框架：**

```
Step 1: SCF 自洽计算
   目的：获得收敛的基态电子密度
   关键设置：out_chg 1（输出电荷密度文件）
   输出：OUT.suffix/SPIN1_CHG.cube

Step 2: NSCF 非自洽计算
   目的：固定电子密度，在指定 k 点求解 KS 方程
   关键设置：
     init_chg  file    （读入 Step 1 的电荷密度）
     symmetry  0       （关闭对称性，保留所有指定 k 点）
     calculation nscf
   输出（能带）：OUT.suffix/BANDS_1.dat
   输出（DOS） ：OUT.suffix/TDOS, PDOS
```

两步计算的关键区别是 **KPT 文件不同**：

| 用途 | KPT 模式 | 说明 |
|------|---------|------|
| SCF | MP 网格（Gamma 或 MP） | 保证电荷密度收敛 |
| NSCF 能带 | Line 模式 | 沿高对称路径采 k 点 |
| NSCF DOS | 密 MP 网格 | k 点越密，DOS 越光滑 |

> **为什么 NSCF 要关闭对称性？** SCF 阶段利用晶体对称性减少计算量是合理的。NSCF 阶段用户手动指定了高对称路径上的 k 点，如果开启对称性，程序会将部分 k 点折叠，导致能带路径不完整。

---

## 二、案例一：硅（Si）的能带结构

### 2.1 体系说明

硅是金刚石立方结构（空间群 Fd$\bar{3}$m），晶格常数 5.4307 Å，每个原胞含 2 个 Si 原子。硅是典型的间接带隙半导体：价带顶位于布里渊区中心（Γ 点），导带底位于 Γ→X 方向上约 0.85 的位置（Δ 点）。

PBE 泛函计算得到的带隙约为 **0.57 eV**，低于实验值 1.17 eV——这是 PBE 的系统性低估，属于正常现象，不是计算错误。

本案例使用 LCAO 基组，结构为 2 原子原胞（primitive cell）。

---

### 2.2 输入文件准备

#### STRU 文件

Si 原胞的结构文件（2 个原子，FCC 原胞向量）：

```
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

文件说明：
- `Si_ONCV_PBE-1.0.upf`：模守恒赝势（ONCV），PBE 泛函
- `Si_gga_8au_60Ry_2s2p1d.orb`：GGA 泛函，截断半径 8 a.u.，DZP 基组
- 晶格向量为 FCC 原胞向量，即 $\frac{a}{2}[0,1,1]$，$\frac{a}{2}[1,0,1]$，$\frac{a}{2}[1,1,0]$（a=5.4307 Å）

> 赝势和轨道文件可从计算环境的 `/abacus-develop/tests/PP_ORB/` 路径获取。

#### INPUT 文件（SCF）

```
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
smearing_sigma  0.01     # 单位 Ry，半导体可取较小值

mixing_type     broyden
mixing_beta     0.4

symmetry        1
out_chg         1        # 输出电荷密度，供 NSCF 步骤读取
```

关键参数说明：

- `out_chg 1`：在 `OUT.Si/` 目录下输出 `SPIN1_CHG.cube`，NSCF 步骤必须用到
- `ecutwfc 50`：截断能（Ry），LCAO 计算中控制电子密度网格精度
- `smearing_sigma 0.01 Ry`（约 0.136 eV）：Si 是半导体，展宽取较小值即可

#### KPT 文件（SCF）

```
K_POINTS
0
Gamma
9 9 9 0 0 0
```

`9×9×9` Gamma 中心网格，对 Si 原胞（2原子）足以给出收敛的电子密度。

#### INPUT 文件（NSCF 能带）

```
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

与 SCF INPUT 相比，主要改动：`calculation nscf`、`init_chg file`、`symmetry 0`，添加 `out_band 1` 和 `out_bandgap 1`，去掉 `out_chg 1`。

#### KPT 文件（NSCF 能带）

Si 原胞属于 FCC 布里渊区，高对称路径选取 Γ-X-U|K-Γ-L-W-X：

```
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
- 第二行 `8`：高对称点数目
- 第三行 `Line`：能带计算专用模式，每行最后一列为到下一高对称点的插值 k 点数
- 最后一行点数填 `1`，表示路径在此结束

> `U`（0.625,0.250,0.625）和 `K`（0.375,0.375,0.750）是 FCC 布里渊区中同一条棱的两端点，坐标略有不同但物理等价，中间段只需 20 个插值点。

---

### 2.3 运行计算

目录结构建议：

```
Si_band/
├── STRU
├── INPUT_scf          INPUT_nscf
├── KPT_scf            KPT_nscf
├── Si_ONCV_PBE-1.0.upf
└── Si_gga_8au_60Ry_2s2p1d.orb
```

**Step 1：SCF 计算**

```bash
cp INPUT_scf INPUT && cp KPT_scf KPT
mpirun -np 4 abacus | tee scf.log
grep "convergence" OUT.Si/running_scf.log   # 应出现：convergence is achieved
```

**Step 2：NSCF 能带计算**

```bash
cp INPUT_nscf INPUT && cp KPT_nscf KPT
mpirun -np 4 abacus | tee nscf.log
```

NSCF 不迭代电荷密度，通常比 SCF 快得多。完成后 `OUT.Si/` 目录下会生成 `BANDS_1.dat` 和 `running_nscf.log`。

---

### 2.4 费米能级提取

```bash
grep "EFERMI" OUT.Si/running_nscf.log
# 输出示例：EFERMI = 6.389520000 eV
```

记录该值，后续配置 `config.json` 时使用。SCF 日志中也可提取（两者应一致）：

```bash
grep "EFERMI" OUT.Si/running_scf.log
```

---

### 2.5 Atomkit 后处理——绘制能带图

将 KPT 文件复制到输出目录（`abacus-plot` 需要它来确定高对称点位置）：

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

- `efermi`：上一步提取的费米能（eV），能带图以此为零点
- `energy_range`：图的纵坐标范围（相对费米能），`[-6, 6]` 通常足够展示主要特征
- `kptfile`：用于标注高对称点名称

```bash
abacus-plot -b
```

生成 `band.png`。

---

### 2.6 结果解读

```bash
grep -i "bandgap" OUT.Si/running_nscf.log
# 示例输出：E_bandgap  0.5700  eV
```

PBE 计算的 Si 带隙约为 **0.57 eV**，间接带隙特征表现为：

- 价带顶（VBM）：Γ 点
- 导带底（CBM）：Γ→X 路径上（Δ 点，约 0.85 × Γ→X 处）

这与实验值 1.17 eV 的差距是 PBE 的系统性低估，并非计算错误。若需更准确的带隙，可使用 HSE06 杂化泛函（见第四章）。

能带图其他特征：
- Γ 点附近：价带最高三条带简并（无 SOC 时）
- 整体带宽约 12 eV（sp 价带区）

---

## 三、案例二：铝（Al）的态密度

### 3.1 体系说明

铝是面心立方结构（FCC），晶格常数 4.0451 Å，空间群 Fm$\bar{3}$m。铝是典型的简单金属：费米能级附近无带隙，DOS 在费米面处连续，导电性好。

本案例目标：
- TDOS：验证金属特征（费米面处 DOS 不为零）
- PDOS：分解到 Al 3s 和 3p 轨道，分析轨道贡献
- 费米能级提取（EFERMI，以计算结果为准）

本案例使用平面波（PW）基组，nspin=1（无磁性）。

---

### 3.2 输入文件准备

#### STRU 文件

Al FCC 4 原子单胞。DOS 计算推荐先用 Atomkit 转换为 1 原子原胞：

```bash
echo -e "2\n202\n101 STRU\n101" | atomkit
# 生成 PRIMCELL.STRU（1 原子 FCC 原胞），在 INPUT 中用 stru_file PRIMCELL.STRU 指定
```

4 原子单胞 STRU：

```
ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf  upf201

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
4.0450551637  0.0000000000  0.0000000000  #latvec1
0.0000000000  4.0450551637  0.0000000000  #latvec2
0.0000000000  0.0000000000  4.0450551637  #latvec3

ATOMIC_POSITIONS
Direct
Al  #label
0   #magnetism
4   #number of atoms
0.0000000000  0.0000000000  0.0000000000  m  0  0  0
0.5000000000  0.5000000000  0.0000000000  m  0  0  0
0.5000000000  0.0000000000  0.5000000000  m  0  0  0
0.0000000000  0.5000000000  0.5000000000  m  0  0  0
```

#### INPUT 文件（SCF）

金属计算必须使用 smearing，否则费米面处占据数不连续，SCF 难以收敛：

```
INPUT_PARAMETERS
suffix          Al
pseudo_dir      ./

calculation     scf
ecutwfc         60
scf_thr         1e-9
scf_nmax        100

basis_type      pw

nspin           1

smearing_method gauss    # 高斯展宽，适用于金属
smearing_sigma  0.02     # 单位 Ry，约 0.272 eV，金属常用值

mixing_type     broyden
mixing_beta     0.4

symmetry        1
out_chg         1
```

关键参数说明：

- `smearing_sigma 0.02 Ry`（约 0.27 eV）：金属取 0.01–0.05 Ry。值越小，基态越准确，但收敛越难；Al 取 0.02 是常用值
- `ecutwfc 60 Ry`：Al 的 ONCV 赝势通常需要 40–80 Ry

#### KPT 文件（SCF）

```
K_POINTS
0
Gamma
12 12 12 0 0 0
```

#### INPUT 文件（NSCF，DOS）

```
INPUT_PARAMETERS
suffix          Al
pseudo_dir      ./

init_chg        file
calculation     nscf
ecutwfc         60
scf_thr         1e-9
scf_nmax        100

basis_type      pw

nspin           1

smearing_method gauss
smearing_sigma  0.02

symmetry        0        # NSCF 必须关闭对称性

out_dos         2        # 1=只输出 TDOS；2=同时输出 TDOS 和 PDOS
dos_sigma       0.07     # DOS 展宽（eV），默认 0.07
```

`out_dos` 取值说明：
- `1`：输出 TDOS，文件为 `OUT.Al/DOS1` 和 `DOS1_smearing.dat`
- `2`：同时输出 TDOS 和 PDOS，额外生成 `OUT.Al/PDOS`（xml 格式）

`dos_sigma`（eV）：控制输出 DOS 曲线的平滑程度，默认 0.07 eV 对多数体系适用；k 点较少时可适当增大（0.1–0.2）。

#### KPT 文件（NSCF，DOS）

```
K_POINTS
0
Gamma
20 20 20 0 0 0
```

DOS 计算需要比 SCF 更密的 k 网格（建议至少 1.5 倍），不使用 Line 模式。

---

### 3.3 运行计算

```
Al_dos/
├── STRU
├── INPUT_scf          INPUT_nscf
├── KPT_scf            KPT_nscf
└── Al_ONCV_PBE-1.0.upf
```

**Step 1：SCF 计算**

```bash
cp INPUT_scf INPUT && cp KPT_scf KPT
mpirun -np 4 abacus | tee scf.log
grep "convergence" OUT.Al/running_scf.log
```

**Step 2：NSCF DOS 计算**

```bash
cp INPUT_nscf INPUT && cp KPT_nscf KPT
mpirun -np 4 abacus | tee nscf.log
```

完成后，`OUT.Al/` 目录下应有 `TDOS`、`DOS1_smearing.dat`（总态密度）和 `PDOS`（投影态密度）。

---

### 3.4 费米能级提取与 TDOS 绘图

**提取费米能级：**

```bash
grep "EFERMI" OUT.Al/running_nscf.log
# 输出示例：EFERMI = 10.963171515 eV
```

这里的 EFERMI 是 ABACUS 内部能量参考系下的绝对值（含赝势贡献），用于绘图时将能量零点移到费米面。

**绘制 TDOS 图：**

```bash
cd OUT.Al/
```

创建 `config.json`：

```json
{
    "tdosfile": "TDOS",
    "efermi": 10.963171515,
    "energy_range": [-10, 10],
    "dos_range": [0, 5],
    "tdosfig": "tdos.png",
    "dpi": 300
}
```

- `energy_range`：横坐标范围（相对费米能，eV），`[-10, 10]` 可展示 Al 的主要 sp 带
- `dos_range`：纵坐标范围，需根据实际 DOS 值调整

```bash
abacus-plot -d
```

生成 `tdos.png`。Al 金属 TDOS 的典型特征：费米能级处 DOS 连续不为零，曲线平滑，无带隙。

---

### 3.5 PDOS 绘图

PDOS 将态密度分解到各元素的各角动量通道，Al 的主要贡献来自 3s（低能区）和 3p（费米面附近）轨道。

修改 `config.json`：

```json
{
    "pdosfile": "PDOS",
    "efermi": 10.963171515,
    "energy_range": [-10, 10],
    "dos_range": [0, 5],
    "figsize": [14, 10],
    "species": {"Al": [0, 1]},
    "pdosfig": "pdos.png"
}
```

`species` 字段：
- `"Al": [0, 1]`：输出 Al 的 l=0（s）和 l=1（p）通道
- l 的对应：`0=s, 1=p, 2=d, 3=f`

```bash
abacus-plot -d -p -o
```

命令参数：`-d` 处理 DOS/PDOS，`-p` 分轨道模式，`-o` 输出数据文件到 `PDOS_FILE/` 目录。生成 `pdos.png`。

---

### 3.6 结果解读

**TDOS 特征：**
- 费米能级附近 DOS 连续，确认 Al 是金属
- 约 -10 eV 处有一个分离的窄峰，对应 Al 3s 态（近自由电子特征）
- -7 eV 至费米面之间是 Al 3p 主带，DOS 较平滑

**PDOS 特征：**
- Al 3s 轨道：贡献主要在 -10 eV 附近的低能区
- Al 3p 轨道：覆盖从约 -7 eV 到费米面，是费米面处 DOS 的主要来源
- Al 是近自由电子金属，s-p 杂化程度低于过渡金属

**费米能量说明：**
`EFERMI`（如 10.963 eV）是 ABACUS 内部参考系下的绝对值，含赝势贡献。绘图时已将其设为零点，因此 DOS 图横坐标是相对费米能：负值为已占据态，正值为空态。

---

## 四、拓展提示

### 4.1 磁性体系（Fe，nspin=2）

对铁磁/反铁磁体系（如 bcc Fe），在 INPUT 中增加 `nspin 2`，并在 STRU 原子行设置初始磁矩：

```
Fe  #label
4.0 #magnetism（初始磁矩，μB）
1   #number of atoms
0.0  0.0  0.0  m  0  0  0
```

`nspin 2` 时，NSCF 输出两套结果：
- 能带：`BANDS_1.dat`（自旋向上）、`BANDS_2.dat`（自旋向下）
- DOS：`DOS1`（自旋向上）、`DOS2`（自旋向下）

两套 DOS 共享同一费米能级，绘图时通常将自旋向下取负值叠加显示。

bcc Fe 参考参数：晶格常数 2.866 Å，收敛磁矩约 2.2 μB，smearing_sigma 建议 0.002 Ry（较小展宽以准确描述自旋极化）。

完整 Al/Fe 算例：
```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
cd abacus-user-guide/examples/dos_band   # 含 Al/ 和 Fe/ 两个目录
```

### 4.2 带隙修正：HSE06

PBE 系统性低估带隙（Si: 0.57 eV vs 实验 1.17 eV）。使用 HSE06 可以接近实验值，在 INPUT 中添加：

```
dft_functional  hse
hse_omega       0.11    # 屏蔽参数（Bohr^-1），HSE06 默认值
```

HSE06 需要 LCAO 基组，计算量约为 PBE 的 5–10 倍，且依赖 LibRI 库编译支持。

### 4.3 高对称路径自动生成

Atomkit 可以从结构文件自动生成标准 KPT，避免手动查找高对称点出错：

```bash
echo -e "3\n301\n3\n101 PRIMCELL.STRU\n0.06" | atomkit
# 生成 KLINES 文件（即 NSCF 能带的 KPT）
# 0.06 为 kspacing（Å⁻¹），值越小能带越平滑
```

### 4.4 常见问题

**NSCF 报"找不到 SPIN1_CHG.cube"：**
1. 检查 NSCF 的 `suffix` 与 SCF 是否完全一致（ABACUS 在 `OUT.suffix/` 下找文件）
2. 若 SCF 和 NSCF 在不同目录，在 INPUT 中设置 `read_file_dir /path/to/OUT.suffix/`
3. 先确认 SCF 收敛：`grep "convergence" OUT.*/running_scf.log`

---

## 附录

### A. 关键参数速查表

#### SCF/NSCF 通用参数

| 参数 | SCF 典型值 | NSCF 典型值 | 说明 |
|------|-----------|------------|------|
| `calculation` | scf | nscf | 计算类型 |
| `symmetry` | 1 | 0 | NSCF 必须为 0 |
| `out_chg` | 1 | 不设 | SCF 必须输出电荷密度 |
| `init_chg` | atomic | file | NSCF 必须读取文件 |
| `out_band` | 不设 | 1 | 输出能带（NSCF） |
| `out_dos` | 不设 | 1 或 2 | 1=TDOS；2=TDOS+PDOS |
| `out_bandgap` | 不设 | 1 | 在 log 中输出带隙 |
| `dos_sigma` | — | 0.07（默认，eV） | DOS 展宽 |

#### 金属 vs 半导体参数对比

| 参数 | 金属（Al） | 半导体（Si） | 说明 |
|------|-----------|------------|------|
| `smearing_method` | gauss | gauss | 金属必须用 smearing |
| `smearing_sigma` | 0.02 Ry | 0.01 Ry | 金属可取较大值 |
| `nspin` | 1（非磁） | 1 | 磁性体系设 2 |
| SCF k 网格 | 12×12×12 | 9×9×9 | 金属需更密 |
| DOS k 网格 | 20×20×20 | 15×15×15 | DOS 比 SCF 更密 |

#### Atomkit 命令速查

| 命令 | 用途 | config.json 必填字段 |
|------|------|---------------------|
| `abacus-plot -b` | 能带图 | bandfile, efermi, kptfile |
| `abacus-plot -d` | TDOS 图 | tdosfile, efermi |
| `abacus-plot -d -p -o` | PDOS 图 | pdosfile, efermi, species |

### B. 参考资料

1. ABACUS 官方文档 - DOS 和能带：https://mcresearch.github.io/abacus-user-guide/abacus-dos.html
2. ABACUS 官方文档 - PDOS：https://mcresearch.github.io/abacus-user-guide/abacus-pdos.html
3. ABACUS 官方文档 - INPUT 参数：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html
4. 布里渊区高对称点参考：Setyawan & Curtarolo, Comp. Mater. Sci. 49, 299–312 (2010)
5. Atomkit 文档：https://vaspkit.com/atomkit.html
