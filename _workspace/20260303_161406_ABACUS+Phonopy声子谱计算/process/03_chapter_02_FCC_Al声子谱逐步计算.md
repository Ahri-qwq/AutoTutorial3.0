# 二、FCC Al 声子谱：逐步计算

## 2.1 准备工作

### 安装 Phonopy

```bash
git clone https://github.com/phonopy/phonopy.git
cd phonopy
python3 setup.py install
```

也可以用 pip 安装：

```bash
pip install phonopy
```

安装完成后执行 `phonopy --version` 确认安装成功。

### 获取案例文件

案例文件可以从 Gitee 下载：

```bash
git clone https://gitee.com/mcresearch/abacus-user-guide.git
```

进入 `abacus-user-guide/examples/interface_Phonopy/1_Al` 目录，即可看到本教程所用的所有文件。

### 目录结构

开始计算前，工作目录的文件结构如下：

```
1_Al/
├── STRU                             # FCC Al 初始结构
├── INPUT                            # ABACUS 计算参数（SCF + 受力计算）
├── KPT                              # k 点设置
├── band.conf                        # Phonopy 声子谱配置
├── plot_pho.gp                      # gnuplot 绘图脚本
└── psp/
    ├── Al_ONCV_PBE-1.0.upf          # Al 赝势
    └── Al_gga_7au_100Ry_4s4p1d.orb  # Al 数值原子轨道
```

`band.conf` 和 `plot_pho.gp` 文件在后续步骤中使用，内容见第五步和第六步。

---

## 2.2 第一步：准备初始结构（STRU）

FCC Al 已经过结构优化，晶格常数为 4.035 Å。STRU 文件内容如下：

```
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf upf201

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
4.03459549706 0 0 #latvec1
0 4.03459549706 0 #latvec2
0 0 4.03459549706 #latvec3

ATOMIC_POSITIONS
Direct

Al #label
0 #magnetism
4 #number of atoms
0  0  0  m  0  0  0
0.5  0.5  0  m  0  0  0
0.5  0  0.5  m  0  0  0
0  0.5  0.5  m  0  0  0
```

几点说明：
- `LATTICE_CONSTANT` 值为 1.88972612546（即 1 Å 对应的 Bohr 数），这是 ABACUS 的单位转换因子。实际晶格常数为 `4.03459549706 × 1.88972612546 Å ≈ 7.625 Bohr ≈ 4.035 Å`，与 FCC Al 的实验值吻合。
- 4 个 Al 原子占据 FCC 常规胞的 4 个顶角/面心位置，坐标以分数坐标表示。
- `m 0 0 0` 中，`m` 是固定标志关键字，后面三个数字分别对应 x、y、z 方向：`0` 表示不约束（原子可自由移动），`1` 表示固定。`0 0 0` 即三个方向均不约束。

**重要提醒：** 在做声子计算前，结构必须充分优化（力收敛到 < 1×10⁻³ eV/Å 量级）。此案例已提供优化好的结构，直接使用即可。

---

## 2.3 第二步：Phonopy 产生微扰构型

在工作目录下执行：

```bash
phonopy -d --dim="2 2 2" --abacus
```

**参数说明：**
- `-d`：生成位移（displacement）构型
- `--dim="2 2 2"`：对三个方向各扩 2 倍，生成 2×2×2 超胞。常规 FCC 胞有 4 个原子，2×2×2 超胞有 32 个原子。
- `--abacus`：使用 ABACUS 的结构文件格式（STRU）

执行后 Phonopy 会在当前目录产生：
- `phonopy_disp.yaml`：位移信息记录文件
- `STRU-001`：微扰构型（对某个 Al 原子施加微小位移后的超胞结构）

由于 FCC 的晶格对称性极高，这里**只产生 1 个微扰构型**。对称性较低的体系可能产生多个（STRU-001、STRU-002、...），每个都需要用 ABACUS 单独计算。关于扩胞大小的选取和收敛性，见第三章 3.2 节。

---

## 2.4 第三步：ABACUS 计算原子受力

现在用 ABACUS 对 `STRU-001` 做 SCF 计算，获取 32 个原子的受力信息。

INPUT 文件完整内容如下：

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix          Al-fcc
calculation     scf
esolver_type    ksdft
symmetry        1
pseudo_dir      ./psp
orbital_dir     ./psp
cal_stress      1
cal_force       1
stru_file       STRU-001

#Parameters (2.Iteration)
ecutwfc         100
scf_thr         1e-7
scf_nmax        50

#Parameters (3.Basis)
basis_type      lcao
gamma_only      0

#Parameters (4.Smearing)
smearing_method mp
smearing_sigma  0.015

#Parameters (5.Mixing)
mixing_type     pulay
mixing_beta     0.7
mixing_gg0      1.5
```

**关键参数说明：**

| 参数 | 值 | 说明 |
|------|----|------|
| `calculation` | `scf` | 做自洽场计算，输出受力 |
| `cal_force` | `1` | **必须开启**，否则不输出原子受力 |
| `cal_stress` | `1` | 同时计算应力（可选，但建议开启） |
| `stru_file` | `STRU-001` | 指定读取 Phonopy 生成的微扰构型，而非原始 STRU |
| `basis_type` | `lcao` | 使用数值原子轨道基组 |
| `gamma_only` | `0` | 超胞计算需要使用多 k 点，不能只用 Γ 点 |
| `symmetry` | `1` | 开启对称性（Phonopy 扩胞后保留了原始结构的对称性信息） |
| `ecutwfc` | `100` | 平面波截断能，单位 Ry |
| `scf_thr` | `1e-7` | SCF 收敛阈值，需要足够小以保证受力精度 |
| `smearing_method` | `mp` | 金属体系使用 Methfessel-Paxton 展宽 |
| `smearing_sigma` | `0.015` | 展宽参数，单位 Ry |
| `mixing_gg0` | `1.5` | Kerker 预处理参数，有助于金属体系收敛 |

**提交计算：**

建议为每个微扰构型创建独立的计算子目录：

```bash
mkdir disp-001
cp INPUT KPT psp/ disp-001/
cp STRU-001 disp-001/
cd disp-001
# 单机多核运行（以8核为例）
mpirun -np 8 abacus > log 2>&1
cd ..
```

如果有多个微扰构型，则同样创建 `disp-002/`、`disp-003/` 等目录，分别提交计算。

**验证受力输出是否正确：** 计算完成后，在 `disp-001/OUT.Al-fcc/running_scf.log` 中搜索关键词 `FORCE`：

```bash
grep "FORCE" disp-001/OUT.Al-fcc/running_scf.log | head -5
```

如果看到原子受力的数值输出，说明力的计算成功。

---

## 2.5 第四步：生成 FORCE_SETS 文件

ABACUS 计算完成后，用以下命令让 Phonopy 读取受力数据：

```bash
phonopy -f ./disp-001/OUT.Al-fcc/running_scf.log
```

如果有多个微扰构型，需要列出所有构型的 `running_scf.log` 路径：

```bash
phonopy -f ./disp-001/OUT.Al-fcc/running_scf.log \
           ./disp-002/OUT.Al-fcc/running_scf.log
```

也可以使用通配符（如 OUT 目录名称不确定时）：

```bash
phonopy -f ./disp-001/OUT*/running_scf.log \
           ./disp-002/OUT*/running_scf.log
```

命令执行成功后，会在当前目录生成 `FORCE_SETS` 文件，其中包含每个位移构型的受力信息。

**常见错误排查：** 如果命令报错，优先检查：
1. ABACUS 计算是否正常结束（`running_scf.log` 最后是否有 `Total Time` 输出）
2. `running_scf.log` 中是否包含力输出（`grep FORCE disp-001/OUT.Al-fcc/running_scf.log`）
3. 路径中的 `OUT.Al-fcc` 是否与 INPUT 中的 `suffix` 参数一致

---

## 2.6 第五步：计算声子谱（band.conf）

准备 `band.conf` 文件，配置声子谱的计算参数：

```
ATOM_NAME = Al
DIM = 2 2 2
MESH = 8 8 8
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
BAND = 1 1 1  1/2 1/2 1  3/8 3/8 3/4  0 0 0   1/2 1/2 1/2
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
```

**参数逐项说明：**

| 参数 | 值 | 说明 |
|------|----|------|
| `ATOM_NAME` | `Al` | 元素名称，与 STRU 中的 `ATOMIC_SPECIES` 一致 |
| `DIM` | `2 2 2` | 扩胞大小，必须与第二步的 `--dim` 一致 |
| `MESH` | `8 8 8` | q 点采样网格，用于计算态密度等热力学量 |
| `PRIMITIVE_AXES` | 见下 | 从常规 FCC 胞到原胞的变换矩阵 |
| `BAND` | 见下 | 声子谱采样路径 |
| `BAND_POINTS` | `101` | 每段路径上的采样点数 |
| `BAND_CONNECTION` | `.TRUE.` | 在能带交叉处辅助连接声子支 |

**关于 `PRIMITIVE_AXES`：**

FCC 的常规胞（输入的 STRU）有 4 个原子，但物理上的原胞只有 1 个原子。`PRIMITIVE_AXES` 告诉 Phonopy 如何从常规胞转换到原胞：

```
0     1/2   1/2
1/2   0     1/2
1/2   1/2   0
```

这是标准的 FCC 原胞转换矩阵。设置了这个参数后，Phonopy 会以原胞的基矢量为声子谱的坐标系，`BAND` 中的坐标也对应原胞倒空间。

**关于 FCC Al 的 k 路径 Γ→X→K→Γ→L：**

`BAND` 中列出了 5 个高对称点的坐标（以原胞倒格矢为单位）：
- `1 1 1` → Γ（等价于原点）
- `1/2 1/2 1` → X
- `3/8 3/8 3/4` → K
- `0 0 0` → Γ
- `1/2 1/2 1/2` → L

这条路径覆盖了 FCC 第一布里渊区的主要高对称点。不同晶格的高对称点不同，可以使用 [SeeK-path](https://www.materialscloud.org/work/tools/seekpath) 工具自动生成对应的 k 路径。

执行计算：

```bash
phonopy -p band.conf --abacus
```

这会产生 `band.yaml` 文件，其中包含各 q 点的声子频率数据。

---

## 2.7 第六步：绘制声子谱

将 `band.yaml` 转换为 gnuplot 格式并绘图：

```bash
phonopy-bandplot --gnuplot > pho.dat
gnuplot plot_pho.gp
```

`plot_pho.gp` 文件内容如下：

```gnuplot
set terminal pngcairo size 1920, 1080 font 'Arial, 36'
set output "Al-FCC_plot.png"

set ylabel 'Frequency (THz)'
set ytics 2
unset key

x1 = 0.13115990
x2 = 0.17753200
x3 = 0.31664810
xmax = 0.43023590
ymin = 0
ymax = 12

set xrange [0:xmax]
set yrange [ymin:ymax]
set xtics ("{/Symbol G}" 0, "X" x1, "K" x2, "{/Symbol G}" x3, "L" xmax)
set arrow 1 nohead from x1,ymin to x1,ymax lt 2
set arrow 2 nohead from x2,ymin to x2,ymax lt 2
set arrow 3 nohead from x3,ymin to x3,ymax lt 2

plot 'pho.dat' using 1:($2) w l lw 3
```

**脚本说明：**

- `x1`、`x2`、`x3`、`xmax`：各高对称点在横轴上的坐标，从 `pho.dat` 第二行读取。
- `set xtics`：设置高对称点标签，`{/Symbol G}` 在 gnuplot 中显示为希腊字母 Γ。
- `set arrow`：在高对称点处画竖线分隔。
- `ymax = 12`：纵轴最大频率 12 THz，涵盖 FCC Al 声子谱的全部频率范围。

执行后生成 `Al-FCC_plot.png`，即为 FCC Al 的声子色散图。

**也可以用 Python/Origin 绘图：** `pho.dat` 的第一列是 k 点路径（横轴），第二列是声子频率（纵轴，单位 THz）。高对称点的位置可以从 `pho.dat` 的第二行读取。
