# 第三章：案例二——铝的态密度

## 三、案例二：铝（Al）的态密度

### 3.1 体系说明

铝是面心立方结构（FCC），晶格常数 4.0451 Å，空间群 Fm$\bar{3}$m。铝是典型的简单金属：费米能级附近无带隙，DOS 在费米面处连续且较高，导电性好。

本案例计算目标：
- TDOS（Total DOS）：Al 的全态密度，验证金属特征（费米面处 DOS 不为零）
- PDOS（Partial DOS）：分解到 Al 3s 和 3p 轨道，分析轨道贡献比例
- 费米能级提取：EFERMI ≈ 10.963 eV（绝对值，以 ABACUS 内部能量零点为参考）

本案例使用平面波（PW）基组，nspin=1（无磁性）。

---

### 3.2 输入文件准备

#### STRU 文件

Al FCC 结构，使用 4 原子常规单胞（conventional cell）。实际 DOS 计算推荐使用 1 原子原胞（primitive cell）以减少计算量，可用 Atomkit 转换：

```bash
# 先准备好 STRU 文件（4原子单胞），再用 atomkit 提取原胞
echo -e "2\n202\n101 STRU\n101" | atomkit
# 生成 PRIMCELL.STRU，即 1 原子 FCC 原胞
```

4 原子单胞的 STRU 文件如下（如直接用 1 原子原胞则跳过转换步骤）：

```
# STRU（Al FCC，4 原子单胞）
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

结构参数说明：
- 晶格常数 `1.88972612546`（Bohr）× `4.0450551637`（Å）= 4.0451 Å，与实验值一致
- 磁矩设为 `0`（Al 无磁性，nspin=1）
- 原子坐标固定（`m 0 0 0`），不进行结构优化

#### INPUT 文件（SCF）

金属 Al 的 SCF 计算需要使用合适的 smearing 方法：

```
# INPUT（Al SCF）
INPUT_PARAMETERS
suffix          Al
pseudo_dir      ./

calculation     scf
ecutwfc         60
scf_thr         1e-9
scf_nmax        100

basis_type      pw

nspin           1        # 非磁性，不区分自旋

smearing_method gauss    # 高斯展宽，适用于金属
smearing_sigma  0.02     # 展宽参数（Ry），约 0.272 eV

mixing_type     broyden
mixing_beta     0.4

symmetry        1
out_chg         1        # 输出电荷密度
```

关键参数说明：

- `smearing_method gauss`：金属计算必须使用 smearing，否则费米面处的占据数不连续，SCF 难以收敛。高斯展宽（gauss/gaussian）是最通用的选择
- `smearing_sigma 0.02`：展宽参数（单位 Ry）。金属一般取 0.01–0.05 Ry；值越小，基态越准确，但收敛越难；值越大，收敛越快，但总能量偏差增大。Al 取 0.02 Ry（约 0.27 eV）是常用值
- `ecutwfc 60`：PW 截断能（Ry），Al 的赝势（ONCV）通常需要 40–80 Ry
- `nspin 1`：无磁性体系，默认值为 1，可省略，此处显式写出以说明

#### KPT 文件（SCF）

```
# KPT（Al SCF，均匀 k 网格）
K_POINTS
0
Gamma
12 12 12 0 0 0
```

Al 金属的费米面精度对 k 点密度敏感，SCF 阶段使用 `12×12×12`（单胞）或 `16×16×16`（原胞）。

#### INPUT 文件（NSCF，DOS）

```
# INPUT（Al NSCF，DOS 计算）
INPUT_PARAMETERS
suffix          Al
pseudo_dir      ./

init_chg        file     # 读取 SCF 电荷密度
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
dos_sigma       0.07     # DOS 展宽（eV），控制输出 DOS 的平滑程度，默认 0.07
```

`out_dos` 的取值：
- `1`：输出总态密度（TDOS），文件为 `OUT.Al/DOS1` 和 `DOS1_smearing.dat`
- `2`：同时输出 TDOS 和 PDOS，额外生成 `OUT.Al/PDOS`（xml 格式）

`dos_sigma` 是高斯展宽参数（eV），默认值 0.07 eV 对多数金属和半导体均适用；k 点较少时可适当增大（0.1–0.2）以使 DOS 曲线更平滑。

#### KPT 文件（NSCF，DOS）

DOS 计算需要比 SCF 更密的 k 网格，以确保对布里渊区的积分精度：

```
# KPT（Al NSCF，DOS 密网格）
K_POINTS
0
Gamma
20 20 20 0 0 0
```

> 对于 DOS 计算，不需要 Line 模式，仍然使用均匀网格，但点数应显著大于 SCF 阶段（建议至少 1.5 倍）。

---

### 3.3 运行计算

**目录结构：**

```
Al_dos/
├── STRU              # 或 PRIMCELL.STRU（原胞）
├── INPUT
├── KPT
├── INPUT_scf
├── INPUT_nscf
├── KPT_scf
├── KPT_nscf
└── Al_ONCV_PBE-1.0.upf
```

**Step 1：SCF 计算**

```bash
cp INPUT_scf INPUT
cp KPT_scf KPT
mpirun -np 4 abacus | tee scf.log
```

检查收敛：

```bash
grep "convergence" OUT.Al/running_scf.log
# 应出现：convergence is achieved
```

**Step 2：NSCF DOS 计算**

```bash
cp INPUT_nscf INPUT
cp KPT_nscf KPT
mpirun -np 4 abacus | tee nscf.log
```

计算完成后，`OUT.Al/` 目录下应有：
- `TDOS`、`DOS1_smearing.dat`：总态密度
- `PDOS`：投影态密度（xml 格式，`out_dos 2` 时生成）

---

### 3.4 费米能级提取与 TDOS 绘图

**提取费米能级：**

```bash
grep "EFERMI" OUT.Al/running_nscf.log
# 输出示例：EFERMI = 10.963171515 eV
```

注意这里的 EFERMI 是 ABACUS 内部能量参考零点（包含赝势贡献）下的绝对值，用于绘图时将能量归零到费米面。

**绘制 TDOS 图：**

进入输出目录，创建 `config.json`：

```bash
cd OUT.Al/
```

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

参数说明：
- `tdosfile`：TDOS 文件路径
- `energy_range`：横坐标范围（相对费米能的 eV），`[-10, 10]` 可展示 Al 的主要 sp 带
- `dos_range`：纵坐标（态密度）的范围，需根据实际 DOS 值调整

运行绘图：

```bash
abacus-plot -d
```

生成 `tdos.png`。Al 金属 TDOS 的典型特征是费米能级处 DOS 值不为零（约 0.3–0.4 states/eV/cell），曲线平滑连续，无带隙。

---

### 3.5 PDOS 绘图

PDOS 将态密度分解到各元素的各角动量通道，Al 的主要贡献来自 3s 轨道（低能区，约 -10 eV）和 3p 轨道（费米面附近）。

在 `OUT.Al/` 目录下，修改 `config.json`：

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

`species` 字段说明：
- `"Al": [0, 1]`：对 Al 元素，输出角动量量子数 `l=0`（s 轨道）和 `l=1`（p 轨道）的投影态密度
- `l` 的对应关系：`0=s, 1=p, 2=d, 3=f`
- 若要输出所有轨道：`"Al": [0, 1, 2]`（Al 无 d 电子，d 通道贡献接近 0）

运行绘图：

```bash
abacus-plot -d -p -o
```

命令参数说明：
- `-d`：处理 DOS/PDOS 相关计算
- `-p`：分轨道输出（partial，即 PDOS 模式）
- `-o`：输出 PDOS 数据文件（在 `PDOS_FILE/` 目录下）

生成 `pdos.png`。

---

### 3.6 结果解读

**TDOS 特征：**
- 费米能级附近 DOS 连续，确认 Al 是金属
- 约 -10 eV 处有一个分离的窄峰，对应 Al 的 3s 态（近自由电子特征）
- -7 eV 至费米面之间是 Al 3p 主带，DOS 较平滑（自由电子特征）

**PDOS 特征：**
- Al 3s 轨道：贡献主要在 -10 eV 附近的低能区
- Al 3p 轨道：贡献覆盖从约 -7 eV 到费米面，是费米面处 DOS 的主要来源
- 两者在某些能量区间有混合（s-p 杂化），但 Al 是近自由电子金属，杂化程度低于过渡金属

**费米能量的物理意义：**
`EFERMI = 10.963 eV` 是 ABACUS 内部参考系下的费米能绝对值，含赝势贡献。绘图时已将其设为能量零点，因此 DOS 图的横坐标是相对费米能，负值为费米面以下的已占据态，正值为费米面以上的空态。

> **磁性金属（Fe）的 DOS 计算：** 若体系有磁性（如 bcc Fe），设置 `nspin 2`，并在 STRU 的原子行设置初始磁矩（Fe 通常设 4.0）。此时 NSCF 计算会输出 `DOS1`（自旋向上）和 `DOS2`（自旋向下）两套文件，分别绘图后叠加对称即可得到自旋极化 DOS。详见第四章进阶提示。
