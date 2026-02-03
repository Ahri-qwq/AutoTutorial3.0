---
title: "从结构文件准备ABACUS输入文件夹"
author: "AutoTutorial 3.0"
date: "2026-02-03"
topic: "从结构文件准备ABACUS输入文件夹"
task_type: "C - 案例驱动教程"
has_case: true
case_file: "data/input/STRU.md"
word_count: ~10000
line_count: 1397
chapters: 6
---

# 从结构文件准备ABACUS输入文件夹

# 前言

## 教程目标

本教程旨在帮助你快速掌握从常见结构文件（CIF、POSCAR等）准备ABACUS完整输入文件夹的方法。通过实际案例，你将学会使用abacustest工具自动配置赝势、轨道和输入参数，大幅简化ABACUS计算的准备工作。

## 适用读者

- 有结构文件，想快速开始ABACUS计算的用户
- 从VASP等其他第一性原理软件转向ABACUS的用户
- 需要批量准备多个结构输入文件的用户
- 不想手动配置赝势和轨道文件的用户

## 前置知识

- 了解第一性原理计算的基本概念
- 熟悉Linux命令行基本操作
- 已安装Python 3.6+环境

## 教程结构

1. **准备工作**：安装abacustest，下载赝势轨道库，配置环境
2. **案例实战**：通过4个实际案例，展示不同场景下的输入文件准备
   - 案例1：从CIF准备MgO的SCF计算
   - 案例2：从CIF准备Fe2O3的磁性计算
   - 案例3：准备结构优化任务
   - 案例4：批量准备多个结构
3. **进阶技巧**：自定义INPUT参数、KPT设置、赝势轨道选择
4. **其他工具简介**：ASE-ABACUS和ATOMKIT的简要介绍
5. **常见问题**：使用过程中可能遇到的问题和解决方法

## 案例说明

本教程的案例来自实际使用场景，所有命令和输出都经过验证。你可以直接复制命令运行，也可以根据自己的需求调整参数。

教程重点在于"如何操作"，不涉及深入的理论推导。如果你需要了解ABACUS的理论背景和参数详解，请参考ABACUS官方文档。

# 前言

## 教程目标

本教程旨在帮助你快速掌握从常见结构文件（CIF、POSCAR等）准备ABACUS完整输入文件夹的方法。通过实际案例，你将学会使用abacustest工具自动配置赝势、轨道和输入参数，大幅简化ABACUS计算的准备工作。

## 适用读者

- 有结构文件，想快速开始ABACUS计算的用户
- 从VASP等其他第一性原理软件转向ABACUS的用户
- 需要批量准备多个结构输入文件的用户
- 不想手动配置赝势和轨道文件的用户

## 前置知识

- 了解第一性原理计算的基本概念
- 熟悉Linux命令行基本操作
- 已安装Python 3.6+环境

## 教程结构

1. **准备工作**：安装abacustest，下载赝势轨道库，配置环境
2. **案例实战**：通过4个实际案例，展示不同场景下的输入文件准备
   - 案例1：从CIF准备MgO的SCF计算
   - 案例2：从CIF准备Fe2O3的磁性计算
   - 案例3：准备结构优化任务
   - 案例4：批量准备多个结构
3. **进阶技巧**：自定义INPUT参数、KPT设置、赝势轨道选择
4. **其他工具简介**：ASE-ABACUS和ATOMKIT的简要介绍
5. **常见问题**：使用过程中可能遇到的问题和解决方法

## 案例说明

本教程的案例来自实际使用场景，所有命令和输出都经过验证。你可以直接复制命令运行，也可以根据自己的需求调整参数。

教程重点在于"如何操作"，不涉及深入的理论推导。如果你需要了解ABACUS的理论背景和参数详解，请参考ABACUS官方文档。
# 第一章：准备工作

从结构文件准备ABACUS输入文件夹，需要先完成工具安装和环境配置。本章将介绍如何安装abacustest工具，下载推荐的赝势轨道库，并配置环境变量。

## 1.1 安装abacustest

abacustest是ABACUS的前后处理工具，支持从结构文件快速准备完整的输入文件夹。安装方式有两种：

### 方式1：通过pip安装（推荐）

```bash
pip install abacustest
```

### 方式2：从源码安装

```bash
git clone https://github.com/pxlxingliang/abacus-test.git
cd abacus-test
pip install .
```

安装完成后，验证安装是否成功：

```bash
abacustest --version
```

如果显示版本号，说明安装成功。

## 1.2 下载APNS赝势轨道库

APNS（ABACUS Pseudopotential-NAO Square）是ABACUS官方推荐的赝势轨道库，经过大规模测试，兼顾精度和效率。abacustest提供了一键下载功能。

在你的工作目录下执行：

```bash
abacustest model inputs --download-pporb
```

命令执行完成后，当前目录会出现两个文件夹：

```bash
$ ls
apns-orbitals-efficiency-v1  apns-pseudopotentials-v1
```

这两个文件夹包含从H到Bi共83种元素的赝势和轨道文件：
- `apns-pseudopotentials-v1`：赝势文件（UPF格式）
- `apns-orbitals-efficiency-v1`：数值原子轨道文件（ORB格式，DZP水平）

**说明：**
- efficiency系列轨道适合结构优化、分子动力学、声子谱等计算
- 如需更高精度（如能带计算、激发态计算），可下载precision系列轨道
- 镧系元素的4f电子被视为核电子，轨道截断半径为8 au

## 1.3 配置环境变量

为了让abacustest自动找到赝势和轨道文件，需要设置环境变量。

### 临时设置（当前终端有效）

```bash
export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1
export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1
```

将 `/your/path/to/` 替换为实际的绝对路径。

### 永久设置（推荐）

将上述命令添加到 `~/.bashrc` 或 `~/.bash_profile` 文件末尾：

```bash
echo 'export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1' >> ~/.bashrc
echo 'export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1' >> ~/.bashrc
source ~/.bashrc
```

验证环境变量是否设置成功：

```bash
echo $ABACUS_PP_PATH
echo $ABACUS_ORB_PATH
```

如果显示正确的路径，说明配置成功。

**注意：**
- 如果只做平面波（PW）计算，只需设置 `ABACUS_PP_PATH`
- LCAO计算需要同时设置 `ABACUS_PP_PATH` 和 `ABACUS_ORB_PATH`
- 也可以在使用abacustest时通过 `--pp` 和 `--orb` 选项手动指定路径，不依赖环境变量

## 1.4 赝势轨道库的文件命名规则

abacustest能够自动识别赝势和轨道文件，前提是文件命名符合规则：

### 规则1：文件名以元素名开头

例如：
- 赝势：`Mg.PD04.PBE.UPF`、`O.upf`
- 轨道：`Mg_gga_10au_100Ry_2s1p.orb`、`O_gga_6au_100Ry_2s2p1d.orb`

### 规则2：使用element.json文件

如果文件名不以元素名开头，可在赝势库目录下创建 `element.json` 文件，内容为：

```json
{
  "Mg": "my_magnesium_pseudopotential.upf",
  "O": "my_oxygen_pseudopotential.upf"
}
```

### 规则3：自动设置ecutwfc（可选）

如果在赝势库目录下提供 `ecutwfc.json` 文件，abacustest会自动设置ecutwfc为体系所有元素的最大值：

```json
{
  "Mg": 100,
  "O": 100
}
```

APNS赝势轨道库已经包含了这些配置文件，可以直接使用。

## 1.5 准备工作检查清单

完成上述步骤后，检查以下内容：

- [ ] abacustest已安装，`abacustest --version` 可以正常运行
- [ ] 已下载APNS赝势轨道库，两个文件夹存在
- [ ] 环境变量已设置，`echo $ABACUS_PP_PATH` 显示正确路径
- [ ] 环境变量已设置，`echo $ABACUS_ORB_PATH` 显示正确路径

如果以上都完成，就可以开始准备ABACUS输入文件了。
# 第二章：案例实战

本章通过4个实际案例，展示如何使用abacustest从结构文件准备ABACUS输入文件夹。案例从简单到复杂，覆盖最常见的使用场景。

## 2.1 案例1：从CIF准备MgO的SCF计算

MgO是最简单的二元化合物之一，适合作为第一个案例。假设你已经有一个MgO的CIF文件 `MgO.cif`。

### 2.1.1 执行命令

在包含 `MgO.cif` 的目录下执行：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --folder-syntax MgO
```

### 2.1.2 参数说明

- `-f MgO.cif`：指定输入的结构文件
- `--ftype cif`：指定结构文件格式为CIF（也支持poscar、vasp、abacus等）
- `--lcao`：使用LCAO基组（数值原子轨道）。如不使用该选项，将使用PW基组（平面波）
- `--folder-syntax MgO`：指定生成的文件夹名称为MgO。如不设置，将从000000开始自动编号

### 2.1.3 生成的文件结构

命令执行完成后，会生成一个 `MgO` 文件夹，结构如下：

```bash
MgO/
├── INPUT
├── STRU
├── KPT
├── Mg_gga_10au_100Ry_2s1p.orb -> /path/to/apns-orbitals-efficiency-v1/Mg_gga_10au_100Ry_2s1p.orb
├── Mg.PD04.PBE.UPF -> /path/to/apns-pseudopotentials-v1/Mg.PD04.PBE.UPF
├── O_gga_6au_100Ry_2s2p1d.orb -> /path/to/apns-orbitals-efficiency-v1/O_gga_6au_100Ry_2s2p1d.orb
├── O.upf -> /path/to/apns-pseudopotentials-v1/O.upf
├── struinfo.txt
└── struinfo.json
```

**文件说明：**
- `INPUT`：ABACUS的输入参数文件
- `STRU`：ABACUS的结构文件
- `KPT`：K点设置文件（如果INPUT中使用kspacing，则不生成）
- `*.orb`：数值原子轨道文件（软链接）
- `*.upf` / `*.UPF`：赝势文件（软链接）
- `struinfo.txt` / `struinfo.json`：记录原始结构文件路径（可删除）

**注意：** 赝势和轨道文件默认为软链接。如需复制文件而非链接，可添加 `--copy-pp-orb` 选项。

### 2.1.4 INPUT文件内容

生成的 `INPUT` 文件内容如下：

```
INPUT_PARAMETERS
calculation     scf
symmetry        1
ecutwfc         100
scf_thr         1e-07
scf_nmax        100
smearing_method gauss
smearing_sigma  0.015
mixing_type     broyden
mixing_beta     0.8
basis_type      lcao
ks_solver       genelpa
precision       double  # or single
#cal_force      1
#cal_stress     1
kspacing        0.14 # unit in 1/bohr
#gamma_only     0
```

**参数解读：**
- `calculation scf`：进行自洽计算
- `symmetry 1`：开启对称性
- `ecutwfc 100`：平面波截断能100 Ry
- `scf_thr 1e-07`：SCF收敛阈值
- `basis_type lcao`：使用LCAO基组
- `ks_solver genelpa`：使用genelpa求解器（适合LCAO）
- `kspacing 0.14`：自动生成K点网格，间距0.14 (1/bohr)

这套参数适合大多数体系的SCF计算，可根据需要修改。

### 2.1.5 STRU文件内容

生成的 `STRU` 文件包含完整的结构信息：

```
ATOMIC_SPECIES
Mg 24.305000 Mg.PD04.PBE.UPF
O 15.999400 O.upf

NUMERICAL_ORBITAL
Mg_gga_10au_100Ry_2s1p.orb
O_gga_6au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    4.21100000000     0.00000000000     0.00000000000
    0.00000000000     4.21100000000     0.00000000000
    0.00000000000     0.00000000000     4.21100000000

ATOMIC_POSITIONS
Cartesian

Mg
0.000000
4
    0.00000000000     0.00000000000     0.00000000000 1 1 1
    2.10550000000     2.10550000000     0.00000000000 1 1 1
    2.10550000000     0.00000000000     2.10550000000 1 1 1
    0.00000000000     2.10550000000     2.10550000000 1 1 1

O
0.000000
4
    2.10550000000     0.00000000000     0.00000000000 1 1 1
    0.00000000000     2.10550000000     0.00000000000 1 1 1
    0.00000000000     0.00000000000     2.10550000000 1 1 1
    2.10550000000     2.10550000000     2.10550000000 1 1 1
```

**结构说明：**
- `ATOMIC_SPECIES`：元素种类、原子量、赝势文件
- `NUMERICAL_ORBITAL`：数值原子轨道文件
- `LATTICE_CONSTANT`：晶格常数（单位：Bohr）
- `LATTICE_VECTORS`：晶格矢量
- `ATOMIC_POSITIONS`：原子坐标（Cartesian坐标系）

### 2.1.6 运行计算

进入 `MgO` 文件夹，使用ABACUS运行计算：

```bash
cd MgO
mpirun -np 4 abacus
```

计算完成后，结果保存在 `OUT.ABACUS/` 目录中。

## 2.2 案例2：从CIF准备Fe2O3的磁性计算

Fe2O3是磁性材料，Fe原子具有磁矩（约4 μB）。磁性材料计算需要开启自旋极化，设置初始磁矩，并可能需要使用DFT+U方法。

### 2.2.1 执行命令

假设你有 `Fe2O3.cif` 文件，执行：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

### 2.2.2 新增参数说明

- `--nspin 2`：开启共线自旋极化（与INPUT中的nspin参数对应）
- `--init_mag Fe 4.0`：为所有Fe原子设置初始磁矩4.0 μB
- `--dftu`：启用DFT+U方法
- `--dftu_param Fe 3.0`：为Fe元素设置Ueff=3.0 eV，默认施加在d轨道上

**多元素设置示例：**

如果体系包含多种磁性元素，可以依次列出：

```bash
abacustest model inputs -f Co2FeAl.cif --ftype cif --lcao --nspin 2 \
  --init_mag Co 1.5 Fe 2.0 \
  --dftu --dftu_param Co 1.0 Fe 3.0
```

### 2.2.3 生成的INPUT文件

磁性计算的INPUT文件会自动调整参数：

```
INPUT_PARAMETERS
calculation     scf
symmetry        0
ecutwfc         100
scf_thr         1e-07
scf_nmax        100
smearing_method gauss
smearing_sigma  0.015
mixing_type     broyden
mixing_beta     0.4
basis_type      lcao
ks_solver       genelpa
precision       double  # or single
#cal_force      1
#cal_stress     1
kspacing        0.14 # unit in 1/bohr
#gamma_only     0
nspin           2
onsite_radius   3
out_mul         1
dft_plus_u      1
orbital_corr    2 -1
hubbard_u       3.0 0
uramping        3.0
mixing_restart  0.001
```

**磁性相关参数：**
- `nspin 2`：共线自旋极化
- `mixing_beta 0.4`：降低混合参数，促进磁性体系收敛
- `onsite_radius 3`：计算原子磁矩的截断半径
- `out_mul 1`：输出Mulliken布居分析，包含原子磁矩
- `dft_plus_u 1`：启用DFT+U
- `orbital_corr 2 -1`：Fe的d轨道使用DFT+U，O不使用（-1表示不使用）
- `hubbard_u 3.0 0`：Fe的Ueff=3.0 eV，O为0
- `uramping 3.0`：U值渐变，促进SCF收敛
- `symmetry 0`：关闭对称性（磁性计算通常需要关闭）

### 2.2.4 生成的STRU文件（部分）

STRU文件中会为Fe原子设置初始磁矩：

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
    ...（省略其他Fe原子）

O
0.000000
18
    4.31436631561     1.34673139667    10.33072323750 1 1 1 mag   0.00000000
    ...（省略其他O原子）
```

**注意：** 每个Fe原子坐标后都有 `mag 4.00000000`，表示初始磁矩为4.0 μB。O原子的磁矩设为0。

### 2.2.5 查看计算结果

计算完成后，可以在 `OUT.ABACUS/mulliken.txt` 中查看每个原子的磁矩：

```bash
cat OUT.ABACUS/mulliken.txt
```

文件中会显示每个原子的电荷和磁矩，以及各轨道（s、p、d）的贡献。

## 2.3 案例3：准备结构优化任务

除了SCF计算，abacustest还支持其他任务类型，如结构优化（relax）和变胞优化（cell-relax）。

### 2.3.1 准备变胞优化任务

假设需要优化MgO的晶胞结构：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --jtype cell-relax --folder-syntax MgO-cellrelax
```

**新增参数：**
- `--jtype cell-relax`：指定任务类型为变胞优化

**支持的任务类型：**
- `scf`：自洽计算（默认）
- `relax`：固定晶胞的结构优化
- `cell-relax`：变胞优化
- `md`：分子动力学

### 2.3.2 生成的INPUT文件

变胞优化的INPUT文件会自动添加优化相关参数：

```
INPUT_PARAMETERS
calculation     cell-relax
symmetry        1
ecutwfc         100
scf_thr         1e-07
scf_nmax        100
smearing_method gauss
smearing_sigma  0.015
mixing_type     broyden
mixing_beta     0.8
basis_type      lcao
ks_solver       genelpa
precision       double  # or single
cal_force       1
cal_stress      1
kspacing        0.14 # unit in 1/bohr
relax_method    cg # or bfgs, bfgs_trad, cg_bfgs, sd, fire
relax_nmax      60
force_thr_ev    0.01  # unit in eV/A
stress_thr      0.5 # unit in kbar
fixed_axes      None # or volume, shape, a, b, c, ab, ac, bc
#gamma_only     0
```

**优化相关参数：**
- `calculation cell-relax`：变胞优化
- `cal_force 1`：计算原子受力
- `cal_stress 1`：计算应力张量
- `relax_method cg`：使用共轭梯度法（适合变胞优化）
- `relax_nmax 60`：最大优化步数
- `force_thr_ev 0.01`：力收敛阈值（eV/Å）
- `stress_thr 0.5`：应力收敛阈值（kbar）
- `fixed_axes None`：不固定任何轴（可选：volume、shape、a、b、c等）

### 2.3.3 relax与cell-relax的区别

- **relax**：只优化原子坐标，晶胞参数固定
  - INPUT中 `calculation relax`
  - 不需要 `cal_stress` 和 `stress_thr`

- **cell-relax**：同时优化原子坐标和晶胞参数
  - INPUT中 `calculation cell-relax`
  - 需要 `cal_stress 1` 和 `stress_thr`

## 2.4 案例4：批量准备多个结构

在某些场景下，需要为多个结构准备相同的输入文件。例如，计算Pd(100)表面不同层数的表面能。

### 2.4.1 准备结构文件

假设有一批不同厚度的Pd(100)表面结构文件：

```bash
$ ls
Pd100_1layer.vasp  Pd100_3layer.vasp  Pd100_5layer.vasp  Pd100_7layer.vasp
Pd100_2layer.vasp  Pd100_4layer.vasp  Pd100_6layer.vasp  Pd100_8layer.vasp
```

### 2.4.2 批量准备命令

使用通配符一次性准备所有结构：

```bash
abacustest model inputs -f Pd100_*layer.vasp --ftype poscar --lcao --jtype relax --folder-syntax "x[:-5]"
```

**参数说明：**
- `-f Pd100_*layer.vasp`：使用通配符匹配所有文件
- `--folder-syntax "x[:-5]"`：使用Python字符串切片语法设置文件夹名称
  - `x` 代表文件名
  - `[:-5]` 表示去掉最后5个字符（即 `.vasp`）
  - 生成的文件夹名为 `Pd100_1layer`、`Pd100_2layer` 等

### 2.4.3 生成的文件夹结构

命令执行后，会生成8个文件夹：

```bash
$ ls
Pd100_1layer/  Pd100_3layer/  Pd100_5layer/  Pd100_7layer/
Pd100_2layer/  Pd100_4layer/  Pd100_6layer/  Pd100_8layer/
```

每个文件夹都包含完整的ABACUS输入文件，参数设置相同。

### 2.4.4 其他文件夹命名方式

**示例1：使用固定前缀**

```bash
--folder-syntax "Pd_surface_x"
```

生成：`Pd_surface_Pd100_1layer.vasp`、`Pd_surface_Pd100_2layer.vasp` 等

**示例2：提取文件名中的数字**

```bash
--folder-syntax "layer_x.split('_')[1]"
```

生成：`layer_1layer`、`layer_2layer` 等

**示例3：使用默认编号**

如果不设置 `--folder-syntax`，将自动编号：

```bash
000000/  000001/  000002/  ...
```

## 2.5 案例小结

通过以上4个案例，我们学会了：

1. **案例1**：准备简单体系的SCF计算（MgO）
2. **案例2**：准备磁性材料计算，设置初始磁矩和DFT+U（Fe2O3）
3. **案例3**：准备结构优化任务（relax、cell-relax）
4. **案例4**：批量准备多个结构，自定义文件夹命名

这些案例覆盖了大多数使用场景。下一章将介绍如何自定义INPUT参数和KPT设置。
# 第三章：进阶技巧

前面的案例使用了abacustest的默认参数设置。本章介绍如何自定义INPUT参数、KPT设置，以及如何选择合适的赝势和轨道。

## 3.1 自定义INPUT参数

abacustest提供的默认INPUT参数适合大多数体系，但有时需要根据具体情况调整。

### 3.1.1 准备INPUT模板文件

假设在Fe2O3的案例中，你想修改以下参数：
- 将 `smearing_sigma` 从0.015降低到0.001
- 将 `mixing_beta` 从0.4降低到0.2

首先创建一个INPUT模板文件 `INPUT_template`：

```
INPUT_PARAMETERS
smearing_sigma     0.001
mixing_beta        0.2
```

**注意：** 模板文件中只需要包含你想修改的参数，其他参数会使用默认值。

### 3.1.2 使用INPUT模板

在abacustest命令中添加 `--input` 选项：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 \
  --input INPUT_template \
  --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

生成的INPUT文件会使用模板中的参数替换默认值，其他参数保持不变。

### 3.1.3 常见需要调整的参数

**收敛控制：**
- `scf_thr`：SCF收敛阈值，默认1e-07。高精度计算可设为1e-08或更小
- `scf_nmax`：最大SCF步数，默认100。难收敛体系可增大到200或更多

**电子展宽：**
- `smearing_sigma`：展宽参数，默认0.015 Ry。金属体系可适当增大，半导体/绝缘体可减小

**混合参数：**
- `mixing_beta`：电荷密度混合系数，默认0.8（非磁）或0.4（磁性）。难收敛体系可减小到0.1-0.3
- `mixing_type`：混合方法，默认broyden。也可选择pulay

**精度控制：**
- `ecutwfc`：平面波截断能，默认100 Ry。高精度计算可增大到120-150 Ry
- `precision`：计算精度，默认double。快速测试可用single

## 3.2 自定义KPT设置

abacustest默认使用 `kspacing` 参数自动生成K点网格。如果需要手动指定K点，可以使用 `--kpt` 选项。

### 3.2.1 使用Gamma中心的K点网格

指定K点网格为5×5×5，以Gamma点为中心：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 \
  --kpt 5 5 5 \
  --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

生成的KPT文件内容为：

```
K_POINTS
0
Gamma
5 5 5 0 0 0
```

**KPT文件格式说明：**
- 第1行：关键词 `K_POINTS`
- 第2行：K点总数（0表示自动生成）
- 第3行：生成方式（`Gamma` 或 `Monkhorst-Pack`）
- 第4行：6个整数
  - 前3个：沿三个倒易基矢方向的K点数目
  - 后3个：K点网格的平移（通常为0 0 0）

### 3.2.2 使用Monkhorst-Pack网格

如果想使用MP网格而非Gamma中心：

准备KPT模板文件 `KPT_template`：

```
K_POINTS
0
Monkhorst-Pack
5 5 5 0 0 0
```

然后在命令中使用：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 \
  --kpt KPT_template \
  --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

### 3.2.3 Gamma点计算

对于大晶胞体系（如表面、缺陷、分子），可以只使用Gamma点：

```bash
abacustest model inputs -f large_system.cif --ftype cif --lcao --kpt 1 1 1
```

或者在INPUT模板中添加：

```
gamma_only     1
```

### 3.2.4 能带计算的K点路径

能带计算需要沿高对称路径设置K点。准备KPT文件 `KPT_band`：

```
K_POINTS
4
Line
0.000 0.000 0.000  50  # Gamma
0.500 0.000 0.500  50  # X
0.500 0.250 0.750  50  # W
0.375 0.375 0.750  1   # K
```

**格式说明：**
- 第2行：高对称点数量
- 第3行：`Line` 表示沿路径插值
- 第4行起：每行为一个高对称点的坐标（分数坐标）和插值点数

## 3.3 同时自定义INPUT和KPT

可以同时使用 `--input` 和 `--kpt` 选项：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 \
  --input INPUT_template \
  --kpt 5 5 5 \
  --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

**效果：**
- INPUT文件使用模板中的参数，并去掉 `kspacing` 参数
- 生成独立的KPT文件，使用指定的K点网格

## 3.4 赝势和轨道的选择

### 3.4.1 使用不同的赝势轨道库

如果不想使用APNS库，可以手动指定赝势和轨道路径：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao \
  --pp /path/to/your/pseudopotentials \
  --orb /path/to/your/orbitals
```

### 3.4.2 赝势类型选择

ABACUS支持多种赝势类型：

**模守恒赝势（推荐）：**
- ONCV赝势：精度高，适合大多数体系
- SG15赝势：另一套常用的ONCV赝势
- PD04赝势：较老的模守恒赝势

**超软赝势：**
- SSSP库：Materials Cloud提供的超软赝势库
- 需要更高的平面波截断能

**选择建议：**
- 一般计算：使用APNS推荐的ONCV赝势
- 高精度计算：测试不同赝势，选择收敛性好的
- 包含过渡金属：注意赝势是否包含半芯态

### 3.4.3 轨道精度选择

APNS提供两套轨道：

**efficiency系列（推荐用于）：**
- 结构优化
- 分子动力学
- 声子谱计算
- 弹性常数计算

**precision系列（推荐用于）：**
- 高精度能量计算
- 能带计算（需要非占据态）
- 激发态计算（TDDFT）

**使用precision系列：**

下载precision轨道库后，设置环境变量：

```bash
export ABACUS_ORB_PATH=/path/to/apns-orbitals-precision-v1
```

或在命令中指定：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao \
  --orb /path/to/apns-orbitals-precision-v1
```

### 3.4.4 轨道截断半径的影响

轨道文件名中包含截断半径信息，例如：
- `Mg_gga_10au_100Ry_2s1p.orb`：截断半径10 au
- `O_gga_6au_100Ry_2s2p1d.orb`：截断半径6 au

**截断半径的选择：**
- 较大的截断半径：精度更高，但计算量更大
- 较小的截断半径：计算更快，但可能损失精度
- APNS库已经平衡了精度和效率，一般不需要调整

## 3.5 进阶技巧小结

本章介绍了：

1. **自定义INPUT参数**：通过模板文件修改默认参数
2. **自定义KPT设置**：手动指定K点网格或使用模板
3. **赝势轨道选择**：如何选择合适的赝势和轨道库

掌握这些技巧后，你可以根据具体需求灵活调整输入文件，而不局限于默认设置。
# 第四章：其他工具简介

除了abacustest，还有其他工具可以用于准备ABACUS输入文件和进行结构转换。本章简要介绍ASE-ABACUS和ATOMKIT两个工具。

## 4.1 ASE-ABACUS接口

ASE（Atomic Simulation Environment）是丹麦技术大学开发的开源原子模拟Python工具库，广泛应用于计算材料科学。ASE-ABACUS是专门为ABACUS开发的接口，独立于ASE主仓库。

### 4.1.1 安装ASE-ABACUS

ASE-ABACUS需要单独下载安装：

```bash
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .
```

**注意：** 通过 `pip install ase` 安装的官方ASE不包含ABACUS接口，必须使用上述方式安装ASE-ABACUS分支。

### 4.1.2 设置环境变量

ASE-ABACUS需要知道赝势和轨道文件的位置。环境变量设置方法参见第一章1.3节。

也可以在Python脚本中设置：

```python
import os
os.environ['ABACUS_PP_PATH'] = '/path/to/pseudopotentials'
os.environ['ABACUS_ORBITAL_PATH'] = '/path/to/orbitals'
```

### 4.1.3 CIF转STRU

使用ASE-ABACUS将CIF文件转换为STRU文件：

```python
from ase.io import read, write
from pathlib import Path

# 读取CIF文件
cif_file = 'MgO.cif'
atoms = read(cif_file, format='cif')

# 指定赝势和轨道文件
pp = {'Mg': 'Mg.PD04.PBE.UPF', 'O': 'O.upf'}
basis = {'Mg': 'Mg_gga_10au_100Ry_2s1p.orb', 'O': 'O_gga_6au_100Ry_2s2p1d.orb'}

# 写入STRU文件
write('STRU', atoms, format='abacus', pp=pp, basis=basis)
```

**说明：**
- `pp` 和 `basis` 字典指定每个元素对应的赝势和轨道文件名
- 文件名需要与环境变量指定的目录中的文件匹配
- 这只生成STRU文件，不会自动生成INPUT和KPT文件

### 4.1.4 STRU转CIF

反向转换也很简单：

```python
from ase.io import read, write

# 读取STRU文件
atoms = read('STRU', format='abacus')

# 写入CIF文件
write('output.cif', atoms, format='cif')
```

### 4.1.5 STRU转POSCAR

转换为VASP的POSCAR格式：

```python
from ase.io import read, write

atoms = read('STRU', format='abacus')
write('POSCAR', atoms, format='vasp')
```

### 4.1.6 ASE-ABACUS的优势

- **灵活性高**：可以在Python脚本中编程控制
- **支持格式多**：ASE支持几十种结构文件格式
- **可调用ABACUS计算**：不仅能转换文件，还能直接调用ABACUS进行计算
- **适合自动化流程**：可以集成到工作流中

### 4.1.7 ASE-ABACUS的局限

- 需要手动指定赝势和轨道文件
- 不会自动生成INPUT和KPT文件
- 需要一定的Python编程基础
- 对于批量处理，需要自己编写脚本

## 4.2 ATOMKIT工具

ATOMKIT是VASPKIT开发团队开发的跨平台建模与结构转换工具，提供命令行交互界面，支持多种材料模拟软件的结构文件格式。

### 4.2.1 下载和安装

从VASPKIT官网下载最新版本：

```bash
wget https://sourceforge.net/projects/vaspkit/files/Binaries/atomkit.0.9.0.linux.x64.tar.gz
tar -zxvf atomkit.0.9.0.linux.x64.tar.gz
cd atomkit.0.9.0.linux.x64
bash setup.sh
```

安装完成后，可以直接运行 `atomkit` 命令。

### 4.2.2 使用ATOMKIT转换结构

ATOMKIT提供交互式界面，使用方式如下：

**启动ATOMKIT：**

```bash
atomkit
```

**或者使用脚本模式：**

创建一个输入脚本 `convert.txt`：

```
1        # 选择功能：结构转换
113      # 选择输入格式：CIF
101      # 选择输出格式：ABACUS STRU
MgO.cif  # 输入文件名
```

然后执行：

```bash
atomkit < convert.txt
```

### 4.2.3 ATOMKIT的功能代码

ATOMKIT的主要功能包括：

**结构转换（功能1）：**
- 支持格式：CIF、POSCAR、XYZ、ABACUS STRU等
- 输入格式代码：
  - 113：CIF
  - 101：VASP POSCAR
  - 102：ABACUS STRU
- 输出格式代码：
  - 101：ABACUS STRU
  - 102：VASP POSCAR
  - 103：CIF

**其他功能：**
- 晶体结构建模
- 结构可视化
- 晶格参数调整
- 原子坐标变换

### 4.2.4 ATOMKIT的优势

- **无需编程**：命令行交互，易于使用
- **功能丰富**：不仅转换格式，还能建模和可视化
- **跨平台**：支持Linux、macOS、Windows
- **快速迭代**：开发活跃，功能持续更新

### 4.2.5 ATOMKIT的局限

- 交互式操作，不适合大规模批量处理
- 不会自动配置赝势和轨道文件
- 不会生成INPUT和KPT文件
- 需要手动输入功能代码，学习曲线稍陡

## 4.3 三种工具对比

| 特性 | abacustest | ASE-ABACUS | ATOMKIT |
|------|-----------|------------|---------|
| **安装方式** | pip安装 | git+pip | 下载解压 |
| **使用方式** | 命令行 | Python脚本 | 交互式 |
| **结构转换** | ✓ | ✓ | ✓ |
| **自动配置赝势轨道** | ✓ | ✗ | ✗ |
| **生成INPUT文件** | ✓ | ✗ | ✗ |
| **生成KPT文件** | ✓ | ✗ | ✗ |
| **批量处理** | ✓ | ✓（需编程） | ✗ |
| **磁性材料设置** | ✓ | ✗ | ✗ |
| **DFT+U设置** | ✓ | ✗ | ✗ |
| **可视化** | ✗ | ✓（需额外工具） | ✓ |
| **编程基础要求** | 低 | 高 | 低 |

## 4.4 工具选择建议

**推荐使用abacustest，如果：**
- 需要快速准备完整的ABACUS输入文件夹
- 需要批量处理多个结构
- 需要自动配置赝势和轨道
- 需要设置磁性材料或DFT+U参数

**推荐使用ASE-ABACUS，如果：**
- 需要在Python脚本中自动化处理
- 需要集成到复杂的工作流中
- 需要调用ABACUS进行计算
- 熟悉Python编程

**推荐使用ATOMKIT，如果：**
- 只需要简单的格式转换
- 需要可视化结构
- 需要建模和调整晶格参数
- 不想编写脚本

**组合使用：**

实际工作中，可以组合使用多个工具：
1. 用ATOMKIT可视化和调整结构
2. 用abacustest快速准备输入文件
3. 用ASE-ABACUS进行后处理和分析

## 4.5 其他工具小结

本章介绍了ASE-ABACUS和ATOMKIT两个工具，它们各有特点：

- **ASE-ABACUS**：适合Python用户，灵活性高，可编程控制
- **ATOMKIT**：适合交互式使用，功能丰富，易于上手

选择合适的工具可以提高工作效率。对于大多数用户，abacustest是最便捷的选择。
# 附录

## 附录A：常见问题

### Q1: abacustest找不到赝势或轨道文件

**问题：** 运行abacustest时提示找不到赝势或轨道文件。

**解决方法：**
1. 检查环境变量是否正确设置：
   ```bash
   echo $ABACUS_PP_PATH
   echo $ABACUS_ORB_PATH
   ```
2. 确认路径是绝对路径，不是相对路径
3. 检查文件名是否以元素名开头，或者是否有element.json文件
4. 也可以在命令中手动指定路径：
   ```bash
   abacustest model inputs -f file.cif --ftype cif --lcao \
     --pp /path/to/pp --orb /path/to/orb
   ```

### Q2: 生成的INPUT文件参数不符合需求

**问题：** 默认的INPUT参数不适合我的体系。

**解决方法：**
1. 创建INPUT模板文件，只包含需要修改的参数
2. 使用 `--input` 选项指定模板文件
3. 模板中的参数会覆盖默认值，其他参数保持不变

### Q3: 磁性计算不收敛

**问题：** Fe2O3等磁性材料的SCF计算不收敛。

**解决方法：**
1. 降低mixing_beta（如0.2或更小）
2. 增加scf_nmax（如200或更多）
3. 检查初始磁矩设置是否合理
4. 尝试使用uramping参数（DFT+U计算）
5. 关闭对称性（symmetry 0）

### Q4: 软链接在其他机器上失效

**问题：** 将生成的文件夹复制到其他机器后，赝势和轨道文件的软链接失效。

**解决方法：**
1. 使用 `--copy-pp-orb` 选项复制文件而非创建软链接：
   ```bash
   abacustest model inputs -f file.cif --ftype cif --lcao --copy-pp-orb
   ```
2. 或者在目标机器上重新设置环境变量并创建软链接

### Q5: 如何为不同元素使用不同的赝势

**问题：** 想为同一元素的不同原子使用不同的赝势（如Fe的不同氧化态）。

**解决方法：**
1. 在STRU文件的ATOMIC_SPECIES部分，将同一元素标记为不同的"元素"：
   ```
   ATOMIC_SPECIES
   Fe1 55.845 Fe_ONCV_PBE-1.2.upf
   Fe2 55.845 Fe_ONCV_PBE-1.2.upf
   ```
2. 在ATOMIC_POSITIONS部分，使用Fe1和Fe2区分不同的Fe原子
3. 可以为Fe1和Fe2设置不同的初始磁矩

### Q6: LCAO和PW基组如何选择

**问题：** 不知道该用LCAO还是PW基组。

**解决方法：**

**使用LCAO（数值原子轨道），如果：**
- 体系较大（>100个原子）
- 需要快速计算
- 做结构优化、分子动力学
- 计算资源有限

**使用PW（平面波），如果：**
- 需要高精度能量
- 体系较小（<50个原子）
- 做收敛性测试
- 需要与其他PW软件对比

**命令区别：**
- LCAO：添加 `--lcao` 选项
- PW：不添加 `--lcao` 选项（默认）

### Q7: 如何查看生成的文件是否正确

**问题：** 不确定生成的INPUT和STRU文件是否正确。

**解决方法：**
1. 检查INPUT文件：
   - `calculation` 参数是否符合任务类型
   - `basis_type` 是否正确（lcao或pw）
   - `nspin` 是否符合磁性设置
2. 检查STRU文件：
   - 原子数量是否正确
   - 晶格参数是否合理
   - 磁性材料的初始磁矩是否设置
3. 可以先用小的scf_nmax（如10）快速测试是否能运行

### Q8: ecutwfc应该设置多大

**问题：** 不知道平面波截断能应该设置多大。

**解决方法：**
1. abacustest默认使用100 Ry，适合大多数体系
2. 如果赝势库提供了ecutwfc.json，会自动设置
3. 高精度计算可以做收敛性测试：
   - 测试80、100、120、150 Ry
   - 选择能量收敛到1 meV以内的值
4. 一般建议：
   - 快速测试：80 Ry
   - 常规计算：100 Ry
   - 高精度：120-150 Ry

### Q9: K点网格如何选择

**问题：** 不知道K点密度应该设置多大。

**解决方法：**
1. abacustest默认使用kspacing=0.14 (1/bohr)，适合大多数体系
2. 可以根据体系调整：
   - 金属：需要更密的K点，kspacing=0.10-0.12
   - 半导体/绝缘体：可以稍疏，kspacing=0.14-0.16
   - 大晶胞：可以更疏，或只用Gamma点
3. 做收敛性测试：
   - 测试不同的K点密度
   - 选择能量收敛到1 meV以内的密度

### Q10: 如何处理包含稀土元素的体系

**问题：** 体系包含镧系元素，不知道如何设置。

**解决方法：**
1. APNS库中镧系元素的4f电子被视为核电子
2. 如果需要考虑4f电子，需要使用其他赝势库
3. 镧系元素通常需要使用DFT+U方法
4. 注意检查赝势文件是否包含所需的元素

### Q11: symmetry参数有哪些选项？

**问题：** 看到案例中有 `symmetry = Analysis`，这是什么意思？

**解决方法：**

symmetry参数的可选值：
- `symmetry 1`：开启对称性（默认）
- `symmetry 0`：关闭对称性
- `symmetry -1`：仅分析对称性，不应用
- `symmetry Analysis`：与-1相同，分析对称性

**使用建议：**
- 一般计算：使用 `symmetry 1`，可以加速计算
- 磁性计算：使用 `symmetry 0`，避免对称性破坏磁结构
- 调试对称性：使用 `symmetry -1` 或 `Analysis`，查看体系的对称性信息

### Q12: 如何使用Cartesian_angstrom坐标？

**问题：** STRU文件中的坐标单位如何设置？

**解决方法：**

**标准方式（推荐）：**
```
LATTICE_CONSTANT
1.889726  # Bohr

ATOMIC_POSITIONS
Cartesian  # 单位为LATTICE_CONSTANT
```

**Cartesian_angstrom方式：**
```
LATTICE_CONSTANT
1.0  # 不起作用

ATOMIC_POSITIONS
Cartesian_angstrom  # 单位直接为埃
```

**说明：**
- `Cartesian_angstrom` 可以避免显式指定LATTICE_CONSTANT
- 坐标直接以埃（Angstrom）为单位，更直观
- 但这一选项在各接口软件中支持不多，不广为人知
- **推荐使用标准方式**，兼容性更好

## 附录B：参考资料

### 官方文档

1. **ABACUS官方文档**
   - 网址：http://abacus.deepmodeling.com
   - 包含完整的INPUT参数列表和使用说明

2. **abacustest GitHub仓库**
   - 网址：https://github.com/pxlxingliang/abacus-test
   - 包含源码、示例和更新日志

3. **APNS赝势轨道库**
   - 网址：https://aissquare.com/datasets
   - 搜索"ABACUS-APNS-PPORBs"

### 相关教程

1. **ABACUS入门教程系列**
   - 结构文件STRU
   - 电子自洽迭代
   - 能带和态密度计算

2. **ASE-ABACUS使用教程**
   - ASE-ABACUS第一章：使用方法简介
   - 结构转换和计算调用

3. **ATOMKIT文档**
   - 网址：https://vaspkit.com/atomkit.html
   - 包含功能代码和使用示例

### 社区资源

1. **ABACUS论坛**
   - 可以提问和交流使用经验

2. **GitHub Issues**
   - 报告bug和提出功能建议

3. **微信公众号**
   - 关注"ABACUS"获取最新教程和更新

## 附录C：命令速查表

### abacustest常用命令

**下载赝势轨道库：**
```bash
abacustest model inputs --download-pporb
```

**准备SCF计算（LCAO）：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao
```

**准备磁性计算：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao \
  --nspin 2 --init_mag Fe 4.0
```

**准备DFT+U计算：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao \
  --nspin 2 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

**准备结构优化：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao --jtype relax
```

**准备变胞优化：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao --jtype cell-relax
```

**自定义INPUT和KPT：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao \
  --input INPUT_template --kpt 5 5 5
```

**批量处理：**
```bash
abacustest model inputs -f *.vasp --ftype poscar --lcao \
  --folder-syntax "x[:-5]"
```

**指定赝势轨道路径：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao \
  --pp /path/to/pp --orb /path/to/orb
```

**复制文件而非软链接：**
```bash
abacustest model inputs -f file.cif --ftype cif --lcao --copy-pp-orb
```

### 常用参数对照表

| 参数 | 说明 | 示例 |
|------|------|------|
| `-f` | 输入文件 | `-f MgO.cif` |
| `--ftype` | 文件格式 | `--ftype cif` |
| `--lcao` | 使用LCAO基组 | `--lcao` |
| `--jtype` | 任务类型 | `--jtype scf` |
| `--nspin` | 自旋极化 | `--nspin 2` |
| `--init_mag` | 初始磁矩 | `--init_mag Fe 4.0` |
| `--dftu` | 启用DFT+U | `--dftu` |
| `--dftu_param` | DFT+U参数 | `--dftu_param Fe 3.0` |
| `--input` | INPUT模板 | `--input INPUT_template` |
| `--kpt` | KPT设置 | `--kpt 5 5 5` |
| `--pp` | 赝势路径 | `--pp /path/to/pp` |
| `--orb` | 轨道路径 | `--orb /path/to/orb` |
| `--folder-syntax` | 文件夹命名 | `--folder-syntax "x[:-5]"` |
| `--copy-pp-orb` | 复制文件 | `--copy-pp-orb` |

## 附录D：进阶学习

### 深入学习ABACUS

1. **参数调优**
   - 学习如何调整SCF收敛参数
   - 理解不同参数对计算精度和效率的影响

2. **高级功能**
   - 能带计算和态密度分析
   - 声子谱和弹性常数计算
   - 分子动力学模拟

3. **后处理分析**
   - 使用Python脚本分析输出结果
   - 可视化电荷密度和能带结构

### 扩展工具学习

1. **ASE-ABACUS深入**
   - 学习如何用ASE调用ABACUS进行计算
   - 构建自动化工作流

2. **高通量计算**
   - 使用abacustest进行批量计算
   - 结合dpdispatcher提交任务到超算

3. **机器学习势能**
   - ABACUS与DeePMD-kit结合
   - 训练和使用深度势能

### 建议学习路径

1. **初级**：掌握本教程的所有案例，能够独立准备输入文件
2. **中级**：学习参数调优，理解不同参数的物理意义
3. **高级**：构建自动化工作流，进行高通量计算
4. **专家**：开发自己的工具，贡献代码到ABACUS社区

---

**教程结束。祝你使用ABACUS顺利！**
