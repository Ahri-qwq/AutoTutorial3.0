# 第三章：MgO 态密度计算完整实战

## 3.1 案例背景与目标

### MgO 材料介绍

MgO（氧化镁）是典型的离子晶体，具有以下特点：

**晶体结构**
- 岩盐结构（NaCl 型）
- 面心立方晶格
- Mg²⁺ 和 O²⁻ 交替排列

**电子结构特征**
- 宽带隙绝缘体，实验带隙约 7.8 eV
- 强离子键：Mg 失去 2 个电子，O 获得 2 个电子
- 价带主要由 O 2p 轨道组成
- 导带主要由 Mg 3s 轨道组成

**为什么选择 MgO**
- 结构简单，只有 2 个原子
- 电子结构特征明显，便于理解
- PDOS 可清晰展示离子键特征
- 计算量适中，适合教学

### 计算目标

本案例将完成以下计算：

1. **SCF 自洽计算**：获得 MgO 的基态电荷密度
2. **NSCF 态密度计算**：计算总 DOS 和 PDOS
3. **结果分析**：
   - 确定带隙大小
   - 分析 Mg 和 O 的轨道贡献
   - 验证离子键特征

### 预期结果

**总 DOS 特征**
- 带隙约 4-5 eV（PBE 泛函低估带隙）
- 价带宽度约 5-6 eV
- 导带起始于费米能级以上约 4-5 eV

**PDOS 特征**
- 价带顶：O 2p 轨道主导
- 导带底：Mg 3s 轨道主导
- Mg 和 O 的 PDOS 分离明显（离子键特征）

## 3.2 输入文件详解

### 3.2.1 INPUT 文件

以下是完整的 INPUT 文件及逐参数说明：

```
INPUT_PARAMETERS

# System variables
calculation         scf          # 自洽场计算
symmetry            0            # 不使用对称性
kspacing            0.08         # k点间距，单位 1/Bohr
precision           double       # 双精度计算

# Plane wave related variables
ecutwfc             100          # 平面波截断能，单位 Ry

# Electronic structure
basis_type          lcao         # LCAO基组（用于输出PDOS）
ks_solver           genelpa      # Kohn-Sham方程求解器
smearing_method     gauss        # 高斯展宽
smearing_sigma      0.015        # 展宽参数，单位 Ry
mixing_type         broyden      # Broyden混合方法
mixing_beta         0.8          # 混合参数
scf_nmax            100          # 最大SCF迭代步数
scf_thr             1e-07        # SCF收敛阈值

# DOS output
init_chg            file         # 从文件读取电荷密度（NSCF用）
out_dos             2            # 输出总DOS和PDOS
```

**参数详细说明：**

| 参数 | 取值 | 说明 |
|------|------|------|
| calculation | scf | SCF自洽计算，求解基态电荷密度 |
| symmetry | 0 | 不使用对称性（NSCF需要完整k点） |
| kspacing | 0.08 | k点间距，自动生成k点网格 |
| precision | double | 双精度计算，提高数值精度 |
| ecutwfc | 100 | 平面波截断能，与轨道文件推荐值一致 |
| basis_type | lcao | LCAO基组，可输出PDOS |
| ks_solver | genelpa | 高效的本征值求解器 |
| smearing_method | gauss | 高斯展宽，适合绝缘体 |
| smearing_sigma | 0.015 | 展宽参数（约0.2 eV） |
| mixing_type | broyden | Broyden混合，加速SCF收敛 |
| mixing_beta | 0.8 | 混合系数 |
| scf_nmax | 100 | 最大迭代步数 |
| scf_thr | 1e-07 | 电荷密度收敛阈值 |
| init_chg | file | 从文件读取电荷密度（NSCF时使用） |
| out_dos | 2 | 输出总DOS和PDOS |

**注意：**
- SCF 计算时，`init_chg` 应设为 `atomic`，`out_dos` 可以不设置
- NSCF 计算时，需要将 `calculation` 改为 `nscf`，`init_chg` 改为 `file`

### 3.2.2 STRU 文件

以下是完整的 STRU 文件：

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
    2.97691954880     0.00000000000     0.00000000000
    1.48845977440     2.57808795428     0.00000000000
    1.48845977440     0.85936265143     2.43064463329

ATOMIC_POSITIONS
Cartesian

Mg
0.000000
1
    0.00000000000     0.00000000000     0.00000000000 1 1 1

O
0.000000
1
    2.97691954880     1.71872530285     1.21532231664 1 1 1
```

**结构说明：**

| 部分 | 内容 | 说明 |
|------|------|------|
| ATOMIC_SPECIES | Mg, O | 元素符号、原子量、赝势文件 |
| NUMERICAL_ORBITAL | 轨道文件 | Mg: 2s1p, O: 2s2p1d |
| LATTICE_CONSTANT | 1.889726 | Bohr转换因子（Bohr到Angstrom） |
| LATTICE_VECTORS | 3个矢量 | 晶格矢量（单位：LATTICE_CONSTANT） |
| ATOMIC_POSITIONS | Cartesian | 原子坐标（笛卡尔坐标系） |

**赝势文件：**
- `Mg.PD04.PBE.UPF`：Mg的ONCV赝势
- `O.upf`：O的赝势

**轨道文件：**
- `Mg_gga_10au_100Ry_2s1p.orb`：Mg的DZP基组（2个s轨道，1个p轨道）
- `O_gga_6au_100Ry_2s2p1d.orb`：O的TZDP基组（2个s轨道，2个p轨道，1个d轨道）

**晶格信息：**
- 晶格常数：约7.64 Bohr（4.04 Å）
- 结构类型：岩盐结构
- 原子数：2个（1个Mg，1个O）

### 3.2.3 KPT 文件

```
K_POINTS
0
Gamma
18 18 18 0 0 0
```

**参数说明：**
- 第一行：`K_POINTS` 关键字
- 第二行：`0` 表示自动生成k点
- 第三行：`Gamma` 表示Gamma中心网格
- 第四行：`18 18 18` 表示三个方向的k点数，`0 0 0` 表示不平移

**k点密度选择：**
- 18×18×18网格对应约5832个k点（考虑对称性后会减少）
- 对于MgO这样的绝缘体，该密度足以获得平滑的DOS曲线

## 3.3 第一步：SCF 自洽计算

### 准备SCF输入文件

将上述INPUT文件修改为SCF版本：

```
INPUT_PARAMETERS

calculation         scf          # SCF自洽计算
symmetry            1            # 使用对称性加速SCF
init_chg            atomic       # 从原子电荷密度开始
out_chg             1            # 输出电荷密度文件

# 其他参数与上面相同
basis_type          lcao
ecutwfc             100
ks_solver           genelpa
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-07
```

**关键修改：**
- `calculation = scf`
- `symmetry = 1`（利用对称性）
- `init_chg = atomic`
- `out_chg = 1`（输出电荷密度）

### 运行SCF计算

```bash
mpirun -np 8 abacus
```

或根据你的计算资源调整进程数。

### 检查SCF收敛性

查看 `OUT.ABACUS/running_scf.log` 文件，关注以下信息：

**收敛过程：**
```
ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)
GE1    -3.456789e+02  0.000000e+00   2.345e-01  1.234e+00
GE2    -3.457123e+02  -3.340000e-02  5.678e-02  1.123e+00
GE3    -3.457145e+02  -2.200000e-03  1.234e-03  1.098e+00
...
GE10   -3.457146e+02  -1.234567e-08  3.456e-08  1.045e+00
```

**收敛判据：**
- `EDIFF`：能量变化，应逐步减小
- `DRHO`：电荷密度变化，应小于 `scf_thr`（1e-07）
- 当 `DRHO < scf_thr` 时，SCF收敛

### 确认输出文件

检查 `OUT.ABACUS/` 目录下是否生成以下文件：

```bash
ls OUT.ABACUS/
```

**必需文件：**
- `SPIN1_CHG`：电荷密度文件（NSCF需要）
- `running_scf.log`：SCF运行日志
- `STRU_READIN_ADJUST.cif`：读入的结构文件

如果 `SPIN1_CHG` 文件存在，说明SCF计算成功，可以进行下一步。

## 3.4 第二步：NSCF 态密度计算

### 修改INPUT文件

将INPUT文件修改为NSCF版本，只需修改3个参数：

```
INPUT_PARAMETERS

calculation         nscf         # 改为NSCF
symmetry            0            # 改为0（不使用对称性）
init_chg            file         # 改为file（读取电荷密度）
out_dos             2            # 输出DOS和PDOS

# 其他参数保持不变
basis_type          lcao
ecutwfc             100
ks_solver           genelpa
smearing_method     gauss
smearing_sigma      0.015
mixing_type         broyden
mixing_beta         0.8
scf_nmax            100
scf_thr             1e-07
```

**SCF vs NSCF 参数对比：**

| 参数 | SCF | NSCF | 说明 |
|------|-----|------|------|
| calculation | scf | nscf | 计算类型 |
| symmetry | 1 | 0 | NSCF需要完整k点 |
| init_chg | atomic | file | NSCF读取电荷密度 |
| out_chg | 1 | 不需要 | NSCF不输出电荷密度 |
| out_dos | 不需要 | 2 | NSCF输出DOS/PDOS |

### 运行NSCF计算

```bash
mpirun -np 8 abacus
```

NSCF计算通常比SCF快，因为只计算一次，不需要迭代。

### 确认DOS/PDOS文件生成

检查 `OUT.ABACUS/` 目录：

```bash
ls OUT.ABACUS/
```

**应生成的文件：**
- `DOS1`：原始DOS数据
- `DOS1_smearing.dat`：展宽后的DOS数据（常用于绘图）
- `TDOS`：总DOS数据
- `PDOS`：PDOS的XML文件
- `running_nscf.log`：NSCF运行日志

**文件说明：**
- `DOS1_smearing.dat`：三列数据（能量、DOS、积分DOS）
- `TDOS`：总DOS，可用于快速绘图
- `PDOS`：XML格式，包含所有轨道的PDOS信息

## 3.5 提取费米能级

### 使用grep命令提取

费米能级记录在运行日志中，可以用grep命令提取：

```bash
grep -i 'efermi' OUT.ABACUS/running_nscf.log
```

**输出示例：**
```
EFERMI = 8.27771540465 eV
```

或者从SCF日志中提取（结果相同）：

```bash
grep -i 'efermi' OUT.ABACUS/running_scf.log
```

### 费米能级的物理意义

对于MgO这样的绝缘体：

**费米能级位置**
- 位于带隙中间
- 价带顶在费米能级以下约2-3 eV
- 导带底在费米能级以上约2-3 eV

**在DOS图中的作用**
- 通常将费米能级设为能量零点
- 负能量：占据态（价带）
- 正能量：未占据态（导带）
- 费米能级处DOS为零（带隙）

### 费米能级在后处理中的应用

绘制DOS图时，需要将能量相对于费米能级平移：

```
E_plot = E_raw - E_fermi
```

abacustest工具会自动完成这一步，无需手动处理。

至此，MgO的DOS/PDOS计算已全部完成，下一章将介绍如何使用abacustest进行结果后处理和可视化。
