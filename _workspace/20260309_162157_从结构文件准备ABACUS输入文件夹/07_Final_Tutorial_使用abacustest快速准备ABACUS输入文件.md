---
title: "使用 abacustest 快速准备 ABACUS 输入文件"
author: "AutoTutorial 3.0"
date: "2026-03-09"
topic: "从结构文件准备ABACUS输入文件夹"
task_type: "C"
has_case: true
word_count: ~600
chapters: 6
---

# 使用 abacustest 快速准备 ABACUS 输入文件

从结构文件开始做第一性原理计算，通常需要经历几个步骤：将 CIF 或 POSCAR 转换为 ABACUS 的 STRU 格式，找到对应的赝势和轨道文件，配置 INPUT 和 KPT 文件。这些操作重复性高，在需要处理多个结构时尤为繁琐。

`abacustest model inputs` 命令将上述流程整合为一步：指定结构文件，工具自动完成结构格式转换、赝势和轨道文件匹配、INPUT 参数配置，生成可直接提交计算的输入文件夹。

本教程通过 MgO、Fe₂O₃、Pd(100) 系列等案例，演示该命令的常见用法。

## 一、安装与准备

### 安装 abacustest

```bash
# 通过 pip 安装
pip install abacustest

# 或从源码安装
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装完成后，执行 `abacustest model inputs -h` 可查看所有可用选项。

### 下载 APNS-pp-orb-v1 赝势轨道库

```bash
abacustest model inputs --download-pporb
```

命令执行完成后，当前目录下会出现两个文件夹：

```bash
$ ls
apns-orbitals-efficiency-v1  apns-pseudopotentials-v1
```

`apns-pseudopotentials-v1` 和 `apns-orbitals-efficiency-v1` 包含从 H 到 Bi 共 83 种元素推荐使用的赝势和轨道文件，均基于 APNS-pp-orb-v1。所有轨道为 DZP 水平；镧系元素的 4f 电子被视为核电子，选用截断半径 8 au 的轨道。

### 设置环境变量

将赝势和轨道路径写入环境变量，后续命令会自动读取：

```bash
export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1
export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1
```

如果不设置环境变量，也可以在每次命令中通过 `--pp` 和 `--orb` 选项显式指定路径。

## 二、基础用法：MgO SCF 计算

以 MgO 的 CIF 文件为例，生成 LCAO 基组的 SCF 计算输入文件夹：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --folder-syntax MgO
```

各参数含义：

| 参数 | 说明 |
|------|------|
| `-f MgO.cif` | 指定输入的结构文件 |
| `--ftype cif` | 结构文件格式（支持 `cif`、`poscar` 等 dpdata 支持的格式） |
| `--lcao` | 使用 LCAO 基组；不加此选项则使用 PW 基组（不需要轨道文件） |
| `--folder-syntax MgO` | 指定生成的任务文件夹名称；不设置则从 `000000` 开始编号 |

命令执行后，目录结构如下：

```
MgO
├── INPUT
├── Mg_gga_10au_100Ry_2s1p.orb -> /home/abc/apns-orbitals-efficiency-v1/Mg_gga_10au_100Ry_2s1p.orb
├── Mg.PD04.PBE.UPF -> /home/abc/apns-pseudopotentials-v1/Mg.PD04.PBE.UPF
├── O_gga_6au_100Ry_2s2p1d.orb -> /home/abc/apns-orbitals-efficiency-v1/O_gga_6au_100Ry_2s2p1d.orb
├── O.upf -> /home/abc/apns-pseudopotentials-v1/O.upf
├── STRU
└── struinfo.txt
MgO.cif
run.sh
setting.json
struinfo.json
```

几点说明：

- 赝势（`.UPF`）和轨道（`.orb`）文件是软链接，指向 `ABACUS_PP_PATH` / `ABACUS_ORB_PATH` 中的原始文件。如需复制文件到目录中（而非软链接），加 `--copy-pp-orb` 选项。
- `run.sh` 和 `setting.json` 用于通过 dpdispatcher 提交计算，不需要可以删除。
- `struinfo.txt` 和 `struinfo.json` 记录了原始结构路径，不需要也可以删除。

生成的 `INPUT` 文件内容如下：

```
calculation     scf
symmetry     1
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.8
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
#cal_force     1
#cal_stress     1
kspacing     0.14 # unit in 1/bohr
#gamma_only     0
```

默认使用 `kspacing = 0.14`（单位 1/Bohr）自动生成 K 点，这组参数对很多体系适用，可按需修改。

## 三、磁性体系：Fe₂O₃ + DFT+U

Fe₂O₃ 是磁性材料，Fe 的磁矩约为 4 μB。计算时需要开启共线自旋极化、设置初始磁矩，并对 Fe 的 d 轨道使用 DFT+U：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

新增参数说明：

| 参数 | 说明 |
|------|------|
| `--nspin 2` | 开启共线自旋极化，对应 INPUT 中的 `nspin = 2` |
| `--init_mag Fe 4.0` | 对所有 Fe 原子设置 4 μB 的初猜磁矩 |
| `--dftu` | 启用 DFT+U |
| `--dftu_param Fe 3.0` | 对 Fe 的 d 轨道设置 Ueff = 3.0 eV |

`--init_mag` 和 `--dftu_param` 均支持同时设置多个元素，依次列出即可。例如 Co₂FeAl 体系：

```bash
abacustest model inputs -f Co2FeAl.cif --ftype cif --lcao --nspin 2 --init_mag Co 1.5 Fe 2.0 --dftu --dftu_param Co 1.0 Fe 3.0
```

### 生成的 INPUT 文件

```
INPUT_PARAMETERS
calculation     scf
symmetry     0
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.4
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
#cal_force     1
#cal_stress     1
kspacing     0.14 # unit in 1/bohr
#gamma_only     0
nspin     2
onsite_radius     3
out_mul     1
dft_plus_u     1
orbital_corr     2 -1
hubbard_u     3.0 0
uramping     3.0
mixing_restart     0.001
```

相比 MgO 案例的 INPUT，主要变化：

- `symmetry 0`：磁性计算通常关闭晶体对称性
- `mixing_beta 0.4`：磁性体系适当降低 mixing_beta 以改善收敛
- `nspin 2`：共线自旋极化
- `onsite_radius 3` + `out_mul 1`：输出各原子的磁矩（结果在 `OUT.ABACUS/mulliken.txt`）
- `dft_plus_u 1`、`orbital_corr 2 -1`、`hubbard_u 3.0 0`：DFT+U 设置（Fe d 轨道 Ueff=3 eV，O 不加 U）
- `uramping 3.0` + `mixing_restart 0.001`：辅助 DFT+U 的 SCF 收敛

### 生成的 STRU 文件

```
ATOMIC_SPECIES
Fe 55.845000 Fe_ONCV_PBE-1.2.upf
O 15.999400 O.upf

NUMERICAL_ORBITAL
Fe_gga_7au_100Ry_4s2p2d1f.orb
O_gga_6au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    5.09190205000     0.00000000000     0.00000000000
   -2.54595102500     4.40971652888     0.00000000000
    0.00000000000     0.00000000000    13.77429765000

ATOMIC_POSITIONS
Cartesian

Fe
0.000000
12
    0.00000000000     0.00000000000     1.99971355156 1 1 1 mag   4.00000000
    0.00000000000     0.00000000000    11.77458409843 1 1 1 mag   4.00000000
    0.00000000000     0.00000000000     4.88743527343 1 1 1 mag   4.00000000
    0.00000000000     0.00000000000     8.88686237657 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963     6.59114610157 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963     2.59171899844 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963     9.47886782343 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963    13.47829492657 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925    11.18257865157 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925     7.18315154844 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925     0.29600272344 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925     4.29542982656 1 1 1 mag   4.00000000

O
0.000000
18
    4.31436631561     1.34673139667    10.33072323750 1 1 1 mag   0.00000000
   -1.76841529061     3.06298513222     3.44357441250 1 1 1 mag   0.00000000
    1.76841529061     3.06298513222    10.33072323750 1 1 1 mag   0.00000000
    0.77753573439     1.34673139667     3.44357441250 1 1 1 mag   0.00000000
    1.55507146878     0.00000000000    10.33072323750 1 1 1 mag   0.00000000
    0.99087955622     4.40971652888     3.44357441250 1 1 1 mag   0.00000000
    1.76841529061     2.81663690629     1.14785813750 1 1 1 mag   0.00000000
    3.32348675939     0.12317411296     8.03500696250 1 1 1 mag   0.00000000
    1.76841529061     0.12317411296     1.14785813750 1 1 1 mag   0.00000000
    3.32348675939     2.81663690629     8.03500696250 1 1 1 mag   0.00000000
    4.10102249378     1.46990550963     1.14785813750 1 1 1 mag   0.00000000
    0.99087955622     1.46990550963     8.03500696250 1 1 1 mag   0.00000000
   -0.77753573439     4.28654241592     5.73929068750 1 1 1 mag   0.00000000
    0.77753573439     1.59307962259    12.62643951250 1 1 1 mag   0.00000000
   -0.77753573439     1.59307962259     5.73929068750 1 1 1 mag   0.00000000
    0.77753573439     4.28654241592    12.62643951250 1 1 1 mag   0.00000000
    1.55507146878     2.93981101925     5.73929068750 1 1 1 mag   0.00000000
   -1.55507146878     2.93981101925    12.62643951250 1 1 1 mag   0.00000000
```

所有 12 个 Fe 原子均设置了 4 μB 的初始磁矩，O 原子磁矩为 0。

## 四、选择计算任务类型

通过 `--jtype` 参数切换计算类型，工具会自动调整 INPUT 中的相关参数。支持的任务类型：`scf`（默认）、`relax`、`cell-relax`、`md`、`band`。

以 MgO 晶胞优化为例：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --jtype cell-relax --folder-syntax MgO-cellrelax
```

生成的 `INPUT` 文件如下：

```
INPUT_PARAMETERS
calculation     cell-relax
symmetry     1
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.8
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
cal_force     1
cal_stress     1
kspacing     0.14 # unit in 1/bohr
relax_method     cg # or bfgs, bfgs_trad, cg_bfgs, sd, fire
relax_nmax     60
force_thr_ev     0.01  # unit in eV/A
stress_thr     0.5 # unit in kbar
fixed_axes     None # or volume, shape, a, b, c, ab, ac, bc; only valid for cell-relax calculation to fix some axes
#gamma_only
```

相比 SCF 的 INPUT，cell-relax 增加了：

- `cal_force 1` + `cal_stress 1`：计算力和应力
- `relax_method cg`：使用 CG 优化算法（适合变胞优化）
- `relax_nmax 60`：最大离子步数
- `force_thr_ev 0.01`：力收敛限（eV/Å）
- `stress_thr 0.5`：应力收敛限（kbar）
- `fixed_axes None`：可选固定某些晶格方向

若使用 `--jtype relax`，则 `calculation` 设为 `relax`，并自动去掉 `cal_stress` 和 `stress_thr`。

## 五、自定义 INPUT 和 K 点

默认情况下，INPUT 中使用 `kspacing = 0.14` 自动生成 K 点。如需自定义参数，可通过 `--input` 指定 INPUT 模板，通过 `--kpt` 显式设置 K 点网格。

以 Fe₂O₃ 为例，将 `smearing_sigma` 改为 0.001，`mixing_beta` 降至 0.2，并使用 5×5×5 K 点。

先准备 INPUT 模板文件（`INPUT_template`），只写需要覆盖的参数：

```
INPUT_PARAMETERS
smearing_sigma     0.001
mixing_beta     0.2
```

然后执行：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --input INPUT_template --kpt 5 5 5 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

`INPUT_template` 中的参数会替换默认值，同时生成以 Gamma 为中心、采样密度为 5×5×5 的 KPT 文件，并去掉 INPUT 中的 `kspacing` 参数。

## 六、批量结构准备

`-f` 参数支持同时传入多个结构文件，配合 `--folder-syntax` 可以批量生成各自的输入文件夹。

以 Pd(100) 表面能随层数的收敛性测试为例，结构文件如下：

```
Pd100_1layer.vasp  Pd100_3layer.vasp  Pd100_5layer.vasp  Pd100_7layer.vasp
Pd100_2layer.vasp  Pd100_4layer.vasp  Pd100_6layer.vasp  Pd100_8layer.vasp
```

使用以下命令为所有结构准备优化计算的输入文件夹：

```bash
abacustest model inputs -f Pd100_*layer.vasp --ftype poscar --lcao --jtype relax --folder-syntax "x[:-5]"
```

`--folder-syntax "x[:-5]"` 中的 `x` 代表结构文件名，`[:-5]` 是 Python 字符串切片，含义为去掉文件名末尾 5 个字符（即 `.vasp`）。执行后，每个结构生成同名的任务文件夹：

```
Pd100_1layer/
Pd100_2layer/
...
Pd100_8layer/
```

除了字符串切片，`--folder-syntax` 支持任意合法的 Python 字符串表达式，`x` 代表原文件名（含扩展名），例如：

- `"x[:-5]"`：去掉 `.vasp` 后缀
- `"x.split('.')[0]"`：取点号前的部分

## 附录

### 常用参数速查表

| 参数 | 说明 | 示例 |
|------|------|------|
| `-f` | 输入结构文件（支持多个） | `-f MgO.cif` |
| `--ftype` | 结构文件格式 | `--ftype cif` / `--ftype poscar` |
| `--lcao` | 使用 LCAO 基组（不加则用 PW） | `--lcao` |
| `--jtype` | 计算类型 | `--jtype scf` / `cell-relax` / `relax` / `md` / `band` |
| `--pp` | 赝势库路径（可用环境变量代替） | `--pp /path/to/pp` |
| `--orb` | 轨道库路径（可用环境变量代替） | `--orb /path/to/orb` |
| `--nspin` | 自旋极化 | `--nspin 2` |
| `--init_mag` | 初始磁矩（元素 + 磁矩值，支持多元素） | `--init_mag Fe 4.0` |
| `--dftu` | 启用 DFT+U | `--dftu` |
| `--dftu_param` | DFT+U 参数（元素 + Ueff，支持多元素） | `--dftu_param Fe 3.0` |
| `--input` | 自定义 INPUT 模板 | `--input INPUT_template` |
| `--kpt` | 显式 K 点网格 | `--kpt 5 5 5` |
| `--folder-syntax` | 输出文件夹命名规则 | `--folder-syntax MgO` 或 `"x[:-5]"` |
| `--copy-pp-orb` | 复制赝势轨道文件（而非软链接） | `--copy-pp-orb` |
| `--download-pporb` | 下载 APNS-pp-orb-v1 库 | `--download-pporb` |

### 参考链接

- abacustest GitHub：https://github.com/pxlxingliang/abacus-test
- ABACUS INPUT 参数文档：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html
- ABACUS STRU 文件格式：https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html
