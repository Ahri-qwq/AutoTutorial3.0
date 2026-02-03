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
