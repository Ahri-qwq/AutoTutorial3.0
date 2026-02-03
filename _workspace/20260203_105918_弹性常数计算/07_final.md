# 第1章：引言

弹性常数是描述材料在弹性形变范围内应力与应变关系的物理量。它反映了材料抵抗外力形变的能力，是材料力学性质的核心参数。在材料设计、结构稳定性分析、相变研究等领域，弹性常数的准确计算至关重要。

## 计算方法

计算弹性常数的主流方法是**应力-应变法**。该方法基于广义胡克定律：对晶体施加一系列微小应变，通过第一性原理计算获得相应的应力，然后通过线性拟合得到弹性张量。相比能量-应变法，应力-应变法具有更高的数值稳定性和计算效率，是Materials Project等数据库的标准方法。

## abacustest工具

ABACUS是一款开源的第一性原理计算软件，支持平面波和数值原子轨道基组。abacustest是ABACUS的辅助工具，封装了一系列自动化工作流，可以简化弹性常数计算的全流程：
- 自动准备ABACUS输入文件
- 自动生成变形结构
- 自动提交计算任务
- 自动后处理并输出弹性张量和弹性模量

使用abacustest，用户只需提供晶体结构文件，即可完成从结构优化到弹性常数计算的全部步骤。

## 本教程内容

本教程将通过两个典型案例，介绍使用abacustest计算晶体弹性常数的完整流程：

1. **Si（立方晶系）**：具有3个独立弹性常数（C₁₁、C₁₂、C₄₄），是最简单的晶系，适合理解基本流程。

2. **金红石型TiO₂（四方晶系）**：具有6个独立弹性常数，展示如何处理更复杂的晶系。

教程将重点讲解：
- 弹性常数的理论基础和晶体对称性的影响
- abacustest工具的核心命令和参数设置
- 完整的计算流程：结构准备、优化、弹性计算、后处理
- 结果分析和验证方法
- 常见问题的排查和解决

本教程假设读者已经熟悉ABACUS的基本使用，包括输入文件格式、K点设置、赝势和轨道文件的准备等。如果对这些内容不熟悉，建议先阅读ABACUS的入门教程。

## 环境准备

在开始计算前，需要准备以下环境：

1. **ABACUS软件**：建议使用v3.10或更高版本
2. **abacustest工具**：可通过pip安装
   ```bash
   pip install abacustest
   ```
3. **环境变量**：设置赝势和轨道文件路径
   ```bash
   export ABACUS_PP_PATH=/path/to/pseudopotentials
   export ABACUS_ORB_PATH=/path/to/orbitals
   ```
4. **Python环境**：需要ase库用于结构生成
   ```bash
   pip install ase
   ```

本教程使用APNS-v1赝势库和efficiency基组，这是ABACUS推荐的高效组合。

准备好环境后，我们将从理论基础开始，逐步深入到具体的计算实践。
# 第2章：理论基础

## 2.1 广义胡克定律与弹性张量

在线弹性范围内，晶体的应力张量 σ 与应变张量 ε 满足线性关系，这就是广义胡克定律：

$$
\sigma_{ij} = C_{ijkl} \varepsilon_{kl}
$$

其中：
- σ_ij：应力张量（二阶对称张量）
- ε_kl：应变张量（二阶对称张量）
- C_ijkl：弹性张量（四阶张量），包含81个分量

由于应力和应变都是对称张量，弹性张量也具有对称性：

$$
C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij}
$$

经过对称性约化，弹性张量的独立分量从81个减少到21个。

## 2.2 Voigt记号

为了简化表示，引入Voigt记号将四阶张量降维为6×6矩阵。映射规则如下：

| 张量指标 | Voigt指标 |
|---------|----------|
| xx (11) | 1 |
| yy (22) | 2 |
| zz (33) | 3 |
| yz (23) | 4 |
| xz (13) | 5 |
| xy (12) | 6 |

使用Voigt记号后，应力和应变可表示为6维向量，弹性张量表示为6×6对称矩阵，广义胡克定律简化为：

$$
\begin{bmatrix}
\sigma_1 \\
\sigma_2 \\
\sigma_3 \\
\sigma_4 \\
\sigma_5 \\
\sigma_6
\end{bmatrix}
=
\begin{bmatrix}
C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16} \\
C_{12} & C_{22} & C_{23} & C_{24} & C_{25} & C_{26} \\
C_{13} & C_{23} & C_{33} & C_{34} & C_{35} & C_{36} \\
C_{14} & C_{24} & C_{34} & C_{44} & C_{45} & C_{46} \\
C_{15} & C_{25} & C_{35} & C_{45} & C_{55} & C_{56} \\
C_{16} & C_{26} & C_{36} & C_{46} & C_{56} & C_{66}
\end{bmatrix}
\begin{bmatrix}
\varepsilon_1 \\
\varepsilon_2 \\
\varepsilon_3 \\
\varepsilon_4 \\
\varepsilon_5 \\
\varepsilon_6
\end{bmatrix}
$$

对于一般晶体，这个6×6矩阵是对称的，包含21个独立分量。

## 2.3 晶体对称性的影响

晶体的点群对称性会进一步约束弹性张量的形式，使独立分量数量减少。不同晶系的独立弹性常数数量如下：

| 晶系 | 独立分量数 | 独立弹性常数 |
|------|-----------|-------------|
| 三斜 (Triclinic) | 21 | 全部21个 |
| 单斜 (Monoclinic) | 13 | - |
| 正交 (Orthorhombic) | 9 | - |
| 六方 (Hexagonal) | 5 | C₁₁, C₁₂, C₁₃, C₃₃, C₄₄ |
| 四方 (Tetragonal) | 6 | C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, C₆₆ |
| 立方 (Cubic) | 3 | C₁₁, C₁₂, C₄₄ |

### 立方晶系（Si）

立方晶系具有最高的对称性，弹性张量只有3个独立分量：

$$
\begin{bmatrix}
C_{11} & C_{12} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{11} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{12} & C_{11} & 0 & 0 & 0 \\
0 & 0 & 0 & C_{44} & 0 & 0 \\
0 & 0 & 0 & 0 & C_{44} & 0 \\
0 & 0 & 0 & 0 & 0 & C_{44}
\end{bmatrix}
$$

其中 C₁₁ = C₂₂ = C₃₃，C₁₂ = C₁₃ = C₂₃，C₄₄ = C₅₅ = C₆₆，其他元素为0。

### 四方晶系（TiO₂）

四方晶系（Laue类型I）有6个独立分量：

$$
\begin{bmatrix}
C_{11} & C_{12} & C_{13} & 0 & 0 & 0 \\
C_{12} & C_{11} & C_{13} & 0 & 0 & 0 \\
C_{13} & C_{13} & C_{33} & 0 & 0 & 0 \\
0 & 0 & 0 & C_{44} & 0 & 0 \\
0 & 0 & 0 & 0 & C_{44} & 0 \\
0 & 0 & 0 & 0 & 0 & C_{66}
\end{bmatrix}
$$

其中 C₁₁ = C₂₂，C₄₄ = C₅₅，但 C₁₁ ≠ C₃₃，C₄₄ ≠ C₆₆。

## 2.4 应力-应变法

应力-应变法是计算弹性常数的推荐方法，其基本原理是：

1. **施加应变**：对优化后的晶体结构施加一系列微小应变（正应变和剪切应变）
2. **计算应力**：使用第一性原理计算每个变形结构的应力张量
3. **线性拟合**：根据广义胡克定律，通过线性拟合应力-应变关系得到弹性常数

相比能量-应变法，应力-应变法具有以下优势：
- **高效**：单次DFT计算输出6个应力分量，而能量法只输出1个标量
- **稳定**：线性拟合比二次拟合对数值噪声更不敏感
- **标准**：Materials Project等数据库的标准方法

### 应变设置

对于3维晶体，需要施加6种应变状态（3种正应变 + 3种剪切应变）：
- 正应变：ε₁、ε₂、ε₃（沿x、y、z方向）
- 剪切应变：ε₄、ε₅、ε₆（yz、xz、xy平面）

对每种应变，通常施加4种不同大小：±0.5%和±1%（即 ε ∈ {-0.01, -0.005, 0.005, 0.01}），共产生24个变形结构。

### 计算流程

1. 优化晶体结构，获得无应力的平衡构型
2. 对平衡构型施加应变，生成变形结构
3. 对每个变形结构进行DFT计算（固定晶格，允许原子弛豫）
4. 提取应力张量
5. 对每种应变类型，线性拟合应力-应变关系，得到对应的弹性常数

## 2.5 晶体取向

弹性张量的具体数值与晶体的取向有关。为了便于比较，Materials Project等数据库通常将晶体旋转到IEEE 176/1987标准规定的标准取向。对于Si、Cu等立方晶系，通常使用正方体形式的惯用胞，坐标轴与晶格矢量重合。

abacustest工具在生成变形结构时，假设输入的晶体结构已经处于合适的取向。对于从Materials Project下载的结构，通常已经满足这一要求。
# 第3章：abacustest工具使用

abacustest是ABACUS的辅助工具，封装了弹性常数计算的完整工作流。本章介绍其核心命令和参数设置。

## 3.1 环境配置

使用abacustest前，需要设置环境变量指定赝势和轨道文件的路径：

```bash
export ABACUS_PP_PATH=/path/to/pseudopotentials
export ABACUS_ORB_PATH=/path/to/orbitals
```

abacustest会自动从这些路径中查找对应元素的赝势和轨道文件。本教程使用APNS-v1赝势库和efficiency基组。

## 3.2 核心命令

abacustest提供三个核心命令用于弹性常数计算：

### 3.2.1 model inputs：准备ABACUS输入文件

```bash
abacustest model inputs -f structure.cif --ftype cif --jtype cell-relax --lcao --folder-syntax output_folder
```

**功能**：根据结构文件自动生成ABACUS输入文件夹，包含INPUT、STRU、KPT等文件。

**主要参数**：

| 参数 | 说明 | 示例 |
|------|------|------|
| `-f` | 结构文件名，可指定多个文件 | `Si.cif` |
| `--ftype` | 结构文件格式 | `cif`, `vasp`, `xyz` |
| `--jtype` | ABACUS任务类型 | `cell-relax`, `scf`, `relax` |
| `--lcao` | 使用LCAO基组（数值原子轨道） | 布尔标志 |
| `--folder-syntax` | 输出文件夹名称 | `Si` |

**输出**：生成包含ABACUS输入文件的文件夹，INPUT文件中的参数适配ABACUS v3.10。

**注意事项**：
- 生成的INPUT参数在多数情况下合理，但部分参数（如nspin、磁矩、DFT+U）可能需要手动调整
- K点设置需要根据体系大小调整，确保收敛

### 3.2.2 model elastic prepare：准备弹性计算

```bash
abacustest model elastic prepare -j input_folder [--norm 0.01] [--shear 0.01] [--norelax]
```

**功能**：在已有的ABACUS输入文件夹基础上，生成一系列变形结构的输入文件夹。

**主要参数**：

| 参数 | 说明 | 默认值 | 推荐值 |
|------|------|--------|--------|
| `-j` | 包含ABACUS输入文件的文件夹 | 必需 | - |
| `--norm` | 最大正应变 | 0.01 | 0.01-0.02 |
| `--shear` | 最大剪切应变 | 0.01 | 0.01-0.02 |
| `--norelax` | 不优化变形结构的原子位置 | False | 通常不使用 |

**输出**：
- `org/`：原始结构的输入文件夹
- `deform-*/`：变形结构的输入文件夹（24个）
  - 3种正应变（x、y、z方向）
  - 3种剪切应变（yz、xz、xy平面）
  - 每种应变4个大小：±0.5%和±1%（当--norm和--shear为0.01时）

**注意事项**：
- ⚠️ **重复执行会删除已有结果**，不要在准备好后重复运行
- 应变大小的选择：过小可能导致数值噪声，过大可能超出线弹性范围。0.01（1%）是合理的默认值
- `--norelax`选项：如果不优化原子位置，计算速度更快，但精度可能降低，通常建议优化

### 3.2.3 model elastic post：后处理

```bash
abacustest model elastic post -j input_folder
```

**功能**：读取所有变形结构的计算结果，拟合应力-应变关系，输出弹性张量和弹性模量。

**输出**：
1. **屏幕输出**：
   - 弹性张量（6×6矩阵）
   - 体模量（Bulk modulus）
   - 剪切模量（Shear modulus）
   - 杨氏模量（Young's modulus）
   - 泊松比（Poisson's ratio）

2. **文件输出**：
   - `metrics.json`：所有计算指标
   - `metrics_elastic.json`：弹性性质的详细结果

## 3.3 完整工作流

使用abacustest计算弹性常数的完整流程：

```bash
# 1. 准备结构优化的输入文件
abacustest model inputs -f structure.cif --ftype cif --jtype cell-relax --lcao --folder-syntax material

# 2. 运行结构优化（需要提交ABACUS计算）
cd material
mpirun -np 4 abacus
cd ..

# 3. 提取优化后的结构，准备SCF计算
mkdir material-elastic
cp material/INPUT material/element* material-elastic/
cp material/OUT.ABACUS/STRU_ION_D material-elastic/STRU

# 4. 修改INPUT文件，将calculation改为scf

# 5. 准备弹性计算的输入文件
abacustest model elastic prepare -j material-elastic

# 6. 提交所有变形结构的计算（需要批量提交）

# 7. 计算完成后，后处理得到弹性常数
abacustest model elastic post -j material-elastic
```

## 3.4 参数调优建议

| 参数 | 建议 | 说明 |
|------|------|------|
| 应变大小 | 0.01-0.02 | 过小噪声大，过大超出线弹性范围 |
| K点密度 | 根据体系调整 | 确保应力收敛，通常需要较密的K点 |
| 能量收敛 | 1e-6 eV | 应力计算对能量收敛要求较高 |
| 力收敛 | 0.01 eV/Å | 结构优化的收敛标准 |
| 应力收敛 | 0.5 kBar | 结构优化的应力收敛标准 |

下一章将通过Si的完整案例，展示这些命令的具体使用。
# 第4章：案例1 - Si（立方晶系）

本章通过Si的完整案例，展示使用abacustest计算立方晶系弹性常数的全流程。Si具有3个独立弹性常数（C₁₁、C₁₂、C₄₄），是理解弹性常数计算的理想案例。

## 4.1 结构准备

### 4.1.1 生成晶体结构

使用ase库生成Si的惯用胞（conventional cell）：

```python
# generate_si_structure.py
from ase.build import bulk

# 生成Si的立方惯用胞
si_conv = bulk('Si', cubic=True)

# 保存为CIF文件
si_conv.write("Si_conv.cif")
```

运行脚本：
```bash
python generate_si_structure.py
```

**为什么使用惯用胞？**
- 立方晶系的弹性常数通常在正方体形式的惯用胞中定义
- 惯用胞的晶格矢量沿坐标轴，便于理解和比较
- Si的惯用胞包含8个原子，是面心立方（FCC）结构的常规表示

### 4.1.2 准备结构优化输入文件

使用abacustest生成ABACUS输入文件：

```bash
abacustest model inputs -f Si_conv.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si
```

生成的文件结构：
```
Si/
├── INPUT
├── STRU
├── KPT
├── Si_ONCV_PBE-1.0.upf  # 赝势文件（自动复制）
└── Si_gga_7au_100Ry_2s2p1d.orb  # 轨道文件（自动复制）
```

### 4.1.3 检查和调整INPUT参数

生成的INPUT文件示例：

```
INPUT_PARAMETERS
#Parameters (General)
suffix                  Si
calculation             cell-relax
ntype                   1
nbands                  auto
symmetry                1
pseudo_dir              ./
orbital_dir             ./

#Parameters (Accuracy)
ecutwfc                 100
scf_thr                 1e-6
scf_nmax                100

#Parameters (Smearing)
smearing_method         gaussian
smearing_sigma          0.002

#Parameters (Mixing)
mixing_type             broyden
mixing_beta             0.7

#Parameters (Relaxation)
relax_nmax              100
force_thr_ev            0.01
stress_thr              0.5
relax_method            cg

#Parameters (Basis)
basis_type              lcao
gamma_only              0
```

**关键参数说明**：
- `calculation = cell-relax`：同时优化原子位置和晶胞
- `ecutwfc = 100`：平面波截断能（Ry），用于数值原子轨道的积分
- `scf_thr = 1e-6`：自洽迭代收敛标准（eV）
- `force_thr_ev = 0.01`：原子受力收敛标准（eV/Å）
- `stress_thr = 0.5`：应力收敛标准（kBar）
- `gamma_only = 0`：不使用Gamma点近似，确保K点采样充分

**可能需要调整的参数**：
- K点设置：检查KPT文件，确保K点密度足够（Si通常使用8×8×8或更密）
- 收敛标准：如果需要更高精度，可以降低scf_thr和force_thr_ev

## 4.2 结构优化

### 4.2.1 提交计算

```bash
cd Si
mpirun -np 4 abacus  # 使用4个MPI进程
cd ..
```

### 4.2.2 检查优化结果

优化完成后，检查OUT.Si/running_cell-relax.log文件，确认：
- 力和应力已收敛
- 最终能量稳定
- 晶格常数合理（Si的实验值约5.43 Å）

优化后的结构保存在OUT.Si/STRU_ION_D文件中。

## 4.3 准备弹性计算

### 4.3.1 创建SCF计算文件夹

```bash
mkdir Si-elastic
cp Si/INPUT Si/Si* Si-elastic/
cp Si/OUT.Si/STRU_ION_D Si-elastic/STRU
```

### 4.3.2 修改INPUT文件

将Si-elastic/INPUT中的calculation改为scf：

```
calculation             scf
```

同时确保：
- `cal_stress = 1`：开启应力计算（abacustest会自动添加）
- `cal_force = 1`：开启力计算

### 4.3.3 生成变形结构

```bash
abacustest model elastic prepare -j Si-elastic
```

生成的文件结构：
```
Si-elastic/
├── INPUT
├── STRU
├── KPT
├── Si*  # 赝势和轨道文件
├── org/  # 原始结构
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   └── Si*
├── deform-xx-0/  # x方向正应变，-1%
├── deform-xx-1/  # x方向正应变，-0.5%
├── deform-xx-2/  # x方向正应变，+0.5%
├── deform-xx-3/  # x方向正应变，+1%
├── deform-yy-0/  # y方向正应变
├── deform-yy-1/
├── deform-yy-2/
├── deform-yy-3/
├── deform-zz-0/  # z方向正应变
├── deform-zz-1/
├── deform-zz-2/
├── deform-zz-3/
├── deform-yz-0/  # yz平面剪切应变
├── deform-yz-1/
├── deform-yz-2/
├── deform-yz-3/
├── deform-xz-0/  # xz平面剪切应变
├── deform-xz-1/
├── deform-xz-2/
├── deform-xz-3/
├── deform-xy-0/  # xy平面剪切应变
├── deform-xy-1/
├── deform-xy-2/
└── deform-xy-3/
```

共25个文件夹（1个原始 + 24个变形）。

**应变设置说明**：
- 默认应变：±0.5%和±1%（--norm和--shear默认为0.01）
- 3种正应变（xx、yy、zz）× 4个应变值 = 12个
- 3种剪切应变（yz、xz、xy）× 4个应变值 = 12个
- 总共24个变形结构

### 4.3.4 参数调整（可选）

如果需要调整应变大小：

```bash
# 使用更大的应变（2%）
abacustest model elastic prepare -j Si-elastic --norm 0.02 --shear 0.02

# 使用不同的正应变和剪切应变
abacustest model elastic prepare -j Si-elastic --norm 0.015 --shear 0.01
```

**应变大小的选择**：
- 太小（<0.5%）：数值噪声影响大
- 太大（>2%）：可能超出线弹性范围
- 推荐：0.01（1%）对大多数材料合适

## 4.4 提交计算

需要对所有25个文件夹提交ABACUS计算。可以使用循环脚本：

```bash
#!/bin/bash
# run_elastic.sh

cd Si-elastic

# 计算原始结构
cd org
mpirun -np 4 abacus > log 2>&1
cd ..

# 计算所有变形结构
for dir in deform-*; do
    cd $dir
    mpirun -np 4 abacus > log 2>&1
    cd ..
done
```

运行脚本：
```bash
bash run_elastic.sh
```

**注意**：
- 实际使用中，建议使用作业调度系统（如SLURM）批量提交
- 每个计算通常需要几分钟到几十分钟，取决于体系大小和计算资源

## 4.5 后处理与结果分析

### 4.5.1 运行后处理

所有计算完成后，运行后处理命令：

```bash
abacustest model elastic post -j Si-elastic
```

### 4.5.2 输出结果

屏幕输出示例：

```
Model: elastic
Postprocessing elastic calculation for job: Si-elastic/
             bulk_modulus  shear_modulus  young_modulus  poisson_ratio
Si-elastic/     90.705919      65.134208     157.664088       0.210302

Si-elastic/     elastic_tensor:
            0             1             2          3          4          5
0  155.464972  5.832639e+01  5.832639e+01   0.000000   0.000000   0.000000
1   58.326393  1.554650e+02  5.832639e+01   0.000000   0.000000   0.000000
2   58.326393  5.832639e+01  1.554650e+02   0.000000   0.000000   0.000000
3    0.000000  0.000000e+00  0.000000e+00  76.177486   0.000000   0.000000
4    0.000000 -2.000000e-10  0.000000e+00   0.000000  76.177486   0.000000
5    0.000000  0.000000e+00 -2.000000e-10   0.000000   0.000000  76.177486

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

### 4.5.3 结果解读

**弹性张量（单位：GPa）**：
- C₁₁ = 155.5 GPa
- C₁₂ = 58.3 GPa
- C₄₄ = 76.2 GPa
- 其他元素接近0（符合立方晶系的对称性）

**弹性模量**：
- 体模量（Bulk modulus）：90.7 GPa
- 剪切模量（Shear modulus）：65.1 GPa
- 杨氏模量（Young's modulus）：157.7 GPa
- 泊松比（Poisson's ratio）：0.21

**立方晶系特征**：
- C₁₁ = C₂₂ = C₃₃（对角元素相等）
- C₁₂ = C₁₃ = C₂₃（非对角元素相等）
- C₄₄ = C₅₅ = C₆₆（剪切模量相等）
- 其他元素为0

### 4.5.4 与Materials Project对比

| 弹性常数 | 本文结果 (GPa) | Materials Project (GPa) | 差异 |
|---------|---------------|------------------------|------|
| C₁₁ | 155.5 | 153 | +1.6% |
| C₁₂ | 58.3 | 57 | +2.3% |
| C₄₄ | 76.2 | 74 | +3.0% |

结果与Materials Project数据接近，差异在3%以内，说明计算合理。差异可能来源于：
- 赝势和基组的选择
- K点采样密度
- 收敛标准
- 应变大小

## 4.6 关键参数总结

| 参数类别 | 参数 | 推荐值 | 说明 |
|---------|------|--------|------|
| 结构优化 | force_thr_ev | 0.01 eV/Å | 力收敛标准 |
| 结构优化 | stress_thr | 0.5 kBar | 应力收敛标准 |
| SCF计算 | scf_thr | 1e-6 eV | 能量收敛标准 |
| 应变设置 | --norm | 0.01 | 正应变大小（1%） |
| 应变设置 | --shear | 0.01 | 剪切应变大小（1%） |
| K点采样 | K点网格 | 8×8×8或更密 | 确保应力收敛 |

下一章将展示四方晶系（TiO₂）的计算，流程与Si类似，但独立弹性常数增加到6个。
# 第5章：案例2 - TiO₂（四方晶系）

本章通过金红石型TiO₂的案例，展示四方晶系弹性常数的计算。四方晶系具有6个独立弹性常数，比立方晶系更复杂。

## 5.1 结构准备与优化

### 5.1.1 获取晶体结构

金红石型TiO₂在Materials Project上的编号为mp-2657。从Materials Project下载CIF文件后，使用与Si相同的流程准备输入文件：

```bash
# 准备结构优化输入文件
abacustest model inputs -f TiO2.cif --ftype cif --lcao --jtype cell-relax --folder-syntax TiO2-rutile
```

生成的文件结构：
```
TiO2-rutile/
├── INPUT
├── STRU
├── KPT
├── Ti_ONCV_PBE-1.0.upf
├── Ti_gga_8au_100Ry_2s2p2d1f.orb
├── O_ONCV_PBE-1.0.upf
└── O_gga_7au_100Ry_2s2p1d.orb
```

### 5.1.2 结构优化

流程同4.2节，提交ABACUS计算：

```bash
cd TiO2-rutile
mpirun -np 4 abacus
cd ..
```

优化完成后，提取优化后的结构：

```bash
mkdir TiO2-rutile-elastic
cp TiO2-rutile/INPUT TiO2-rutile/Ti* TiO2-rutile/O* TiO2-rutile-elastic/
cp TiO2-rutile/OUT.TiO2-rutile/STRU_ION_D TiO2-rutile-elastic/STRU
```

修改TiO2-rutile-elastic/INPUT，将calculation改为scf。

## 5.2 弹性计算

### 5.2.1 生成变形结构

流程同4.3节：

```bash
abacustest model elastic prepare -j TiO2-rutile-elastic
```

同样生成25个文件夹（1个原始 + 24个变形）。

### 5.2.2 提交计算

使用与Si相同的批量提交脚本（参考4.4节），对所有文件夹提交计算。

## 5.3 结果分析

### 5.3.1 后处理

```bash
abacustest model elastic post -j TiO2-rutile-elastic
```

### 5.3.2 输出结果

屏幕输出示例：

```
Model: elastic
Postprocessing elastic calculation for job: TiO2-rutile-elastic/
                      bulk_modulus  shear_modulus  young_modulus  poisson_ratio
TiO2-rutile-elastic/    220.705174     128.683753     323.230608       0.255911

TiO2-rutile-elastic/    elastic_tensor:
              0             1             2           3           4           5
0  2.812277e+02  1.581445e+02  1.559258e+02    0.000000    0.000000    0.000000
1  1.581445e+02  2.812277e+02  1.559258e+02    0.000000    0.000000    0.000000
2  1.558677e+02  1.558677e+02  4.840152e+02    0.000000    0.000000    0.000000
3  4.200000e-09  5.500000e-09  1.230000e-08  117.414353    0.000000    0.000000
4 -2.250000e-08 -1.770000e-08 -4.810000e-08    0.000000  117.414354    0.000000
5 -1.300000e-09  2.200002e-09  0.000000e+00    0.000000    0.000000  216.431859

The postprocess is done. The metrics are saved in 'metrics.json', and the elastic results are saved in 'metrics_elastic.json'.
```

### 5.3.3 四方晶系弹性张量

**6个独立弹性常数（单位：GPa）**：
- C₁₁ = 281.2 GPa
- C₁₂ = 158.1 GPa
- C₁₃ = 155.9 GPa
- C₃₃ = 484.0 GPa
- C₄₄ = 117.4 GPa
- C₆₆ = 216.4 GPa

**四方晶系特征**：
- C₁₁ = C₂₂（面内方向相等）
- C₁₃ = C₂₃（面内-面外耦合相等）
- C₄₄ = C₅₅（面外剪切模量相等）
- C₃₃ ≠ C₁₁（面外方向与面内方向不同）
- C₆₆ ≠ C₄₄（面内剪切与面外剪切不同）
- 其他元素接近0

**与立方晶系的差异**：
- 立方晶系：C₁₁ = C₂₂ = C₃₃，C₄₄ = C₅₅ = C₆₆
- 四方晶系：C₁₁ = C₂₂ ≠ C₃₃，C₄₄ = C₅₅ ≠ C₆₆
- 四方晶系沿c轴（z方向）的性质与面内（xy平面）不同，体现了各向异性

## 5.4 与实验和数据库对比

### 5.4.1 对比表格

| 弹性常数 | 本文结果 (GPa) | Materials Project (GPa) | 实验测量 (GPa) |
|---------|---------------|------------------------|---------------|
| C₁₁ | 281.2 | 426 | 268.0 (±1.4) |
| C₁₂ | 158.1 | 2 | 174.9 (±1.4) |
| C₁₃ | 155.9 | 149 | 147.4 (±1.5) |
| C₃₃ | 484.0 | 470 | 484.2 (±1.8) |
| C₄₄ | 117.4 | 113 | 123.8 (±0.2) |
| C₆₆ | 216.4 | 43 | 190.2 (±0.5) |

实验数据来源：D. G. Isaak et al., Physics and Chemistry of Minerals **26**, 31 (1998)

### 5.4.2 结果分析

**本文结果的特点**：
- C₁₁、C₃₃、C₁₃、C₄₄与实验值接近（误差<10%）
- C₁₂和C₆₆与实验值有一定差异，但比Materials Project更接近
- 整体上，本文结果比Materials Project更符合实验测量

**可能的差异来源**：
1. **赝势和基组**：不同的赝势和基组会影响结果
2. **K点采样**：TiO₂需要较密的K点确保收敛
3. **应变大小**：不同的应变大小可能影响拟合精度
4. **温度效应**：DFT计算是0K结果，实验通常在室温

**Materials Project结果异常的原因**：
- C₁₂ = 2 GPa和C₆₆ = 43 GPa明显偏低，可能是计算设置或后处理问题
- 这说明弹性常数计算对参数设置敏感，需要仔细验证

## 5.5 四方晶系的计算要点

### 5.5.1 与立方晶系的差异

| 特性 | 立方晶系 (Si) | 四方晶系 (TiO₂) |
|------|--------------|----------------|
| 独立弹性常数 | 3个 | 6个 |
| 对称性 | 最高 | 中等 |
| 各向异性 | 各向同性 | 面内-面外各向异性 |
| 计算复杂度 | 低 | 中等 |
| 结果验证 | 检查C₁₁=C₂₂=C₃₃ | 检查C₁₁=C₂₂, C₄₄=C₅₅ |

### 5.5.2 注意事项

1. **晶体取向**：确保c轴沿z方向，a、b轴在xy平面
2. **K点采样**：四方晶系通常需要kx=ky≠kz的K点网格
3. **对称性检查**：后处理时检查C₁₁是否等于C₂₂，C₄₄是否等于C₅₅
4. **收敛性测试**：对于各向异性材料，应力收敛要求更高

## 5.6 其他晶系的计算

abacustest支持所有晶系的弹性常数计算，流程与Si和TiO₂类似：

| 晶系 | 独立分量 | 典型材料 | 注意事项 |
|------|---------|---------|---------|
| 立方 | 3 | Si, Cu, NaCl | 最简单，各向同性 |
| 六方 | 5 | Graphite, ZnO | c轴沿z方向 |
| 四方 | 6 | TiO₂, SnO₂ | c轴沿z方向 |
| 正交 | 9 | α-U, CaSO₄ | 三个轴互相垂直 |
| 单斜 | 13 | β-S, CaSO₄·2H₂O | 一个斜角 |
| 三斜 | 21 | CuSO₄·5H₂O | 最复杂 |

对于更复杂的晶系，建议：
- 仔细检查晶体取向
- 使用更密的K点
- 降低收敛标准
- 与文献或数据库对比验证

下一章将介绍后处理方法和常见问题的排查。
# 第6章：后处理与结果分析

本章介绍如何解读abacustest的输出结果，验证计算的正确性，以及排查常见问题。

## 6.1 输出文件解读

### 6.1.1 输出文件

abacustest model elastic post命令生成两个JSON文件：

**metrics.json**：
- 包含所有计算指标
- 格式：标准JSON，可用Python读取

**metrics_elastic.json**：
- 专门存储弹性性质
- 包含弹性张量、弹性模量、泊松比等

### 6.1.2 弹性张量格式

弹性张量以6×6矩阵形式输出，单位为GPa：

```
            0             1             2          3          4          5
0  C11           C12           C13           C14        C15        C16
1  C12           C22           C23           C24        C25        C26
2  C13           C23           C33           C34        C35        C36
3  C14           C24           C34           C44        C45        C46
4  C15           C25           C35           C45        C55        C56
5  C16           C26           C36           C46        C56        C66
```

对于立方晶系，只有C₁₁、C₁₂、C₄₄非零。对于四方晶系，C₁₁、C₁₂、C₁₃、C₃₃、C₄₄、C₆₆非零。

### 6.1.3 弹性模量

**体模量（Bulk modulus, K）**：
- 描述材料抵抗均匀压缩的能力
- 单位：GPa
- 对于立方晶系：K = (C₁₁ + 2C₁₂) / 3

**剪切模量（Shear modulus, G）**：
- 描述材料抵抗剪切变形的能力
- 单位：GPa
- 计算方法：Voigt-Reuss-Hill平均

**杨氏模量（Young's modulus, E）**：
- 描述材料在拉伸或压缩时的刚度
- 单位：GPa
- 关系：E = 9KG / (3K + G)

**泊松比（Poisson's ratio, ν）**：
- 描述材料横向应变与纵向应变的比值
- 无量纲
- 关系：ν = (3K - 2G) / (6K + 2G)
- 范围：通常在0到0.5之间

## 6.2 结果验证

### 6.2.1 对称性检查

根据晶系检查弹性张量是否满足对称性要求：

**立方晶系**：
- C₁₁ = C₂₂ = C₃₃（误差<1%）
- C₁₂ = C₁₃ = C₂₃（误差<1%）
- C₄₄ = C₅₅ = C₆₆（误差<1%）
- 其他元素接近0（|C_ij| < 1 GPa）

**四方晶系**：
- C₁₁ = C₂₂（误差<1%）
- C₁₃ = C₂₃（误差<1%）
- C₄₄ = C₅₅（误差<1%）
- C₁₄ = C₁₅ = C₂₄ = C₂₅ = C₃₄ = C₃₅ = C₄₆ = C₅₆ ≈ 0

如果对称性不满足，可能原因：
- 晶体结构未充分优化
- K点采样不足
- 应变过大，超出线弹性范围
- 数值噪声

### 6.2.2 数值合理性

**弹性常数的合理范围**：
- C₁₁、C₂₂、C₃₃：通常在50-500 GPa
- C₁₂、C₁₃、C₂₃：通常小于C₁₁
- C₄₄、C₅₅、C₆₆：通常在10-200 GPa
- 所有弹性常数应为正值（力学稳定性）

**弹性模量的合理范围**：
- 体模量：10-400 GPa（金属和陶瓷）
- 剪切模量：5-200 GPa
- 杨氏模量：10-500 GPa
- 泊松比：0.1-0.5（大多数材料在0.2-0.4）

**异常值的可能原因**：
- 负的弹性常数：结构不稳定或计算错误
- 过大的弹性常数（>1000 GPa）：应变过小或拟合错误
- 过小的弹性常数（<10 GPa）：应变过大或结构问题

### 6.2.3 与文献对比

**对比来源**：
1. **Materials Project**：https://materialsproject.org/
   - 优点：数据丰富，易于查询
   - 缺点：部分数据可能有误（如TiO₂的C₁₂和C₆₆）

2. **实验测量**：
   - 优点：最可靠的参考
   - 缺点：数据较少，可能有温度效应

3. **其他DFT计算**：
   - 优点：可比较不同方法的差异
   - 缺点：需要注意计算设置的差异

**合理的误差范围**：
- 与实验值：±10%通常可接受
- 与其他DFT计算：±5%通常可接受
- 如果误差过大，需要检查计算设置

## 6.3 常见问题排查

### 6.3.1 对称性不符

**问题**：弹性张量不满足晶系的对称性要求

**可能原因**：
1. 结构优化不充分
2. K点采样不足
3. 应变过大
4. 晶体取向错误

**解决方法**：
```bash
# 1. 检查结构优化是否收敛
grep "FINAL_ETOT" Si/OUT.Si/running_cell-relax.log

# 2. 增加K点密度
# 修改KPT文件，例如从8×8×8改为12×12×12

# 3. 减小应变
abacustest model elastic prepare -j Si-elastic --norm 0.005 --shear 0.005

# 4. 检查晶体取向
# 使用ase或pymatgen检查晶格矢量方向
```

### 6.3.2 数值不稳定

**问题**：多次计算结果差异较大，或拟合R²值较低

**可能原因**：
1. 应力计算精度不足
2. K点采样不足
3. 能量收敛标准过松
4. 应变过小，数值噪声大

**解决方法**：
```bash
# 1. 提高能量收敛标准
# 在INPUT中设置：
scf_thr  1e-7  # 从1e-6降低到1e-7

# 2. 增加K点密度
# 对于应力计算，通常需要比能量计算更密的K点

# 3. 增大应变
abacustest model elastic prepare -j Si-elastic --norm 0.015 --shear 0.015

# 4. 检查应力收敛
# 查看OUT.*/running_scf.log中的应力值是否稳定
```

### 6.3.3 结果与文献差异大

**问题**：计算结果与文献或数据库差异超过20%

**可能原因**：
1. 赝势和基组不同
2. 交换关联泛函不同
3. 结构不同（原胞vs惯用胞）
4. 温度效应（DFT是0K，实验通常室温）

**解决方法**：
```bash
# 1. 检查赝势和基组
# 尝试使用与文献相同的赝势

# 2. 检查交换关联泛函
# 在INPUT中设置：
dft_functional  PBE  # 或LDA、HSE等

# 3. 检查结构
# 确保使用与文献相同的晶体结构（原胞或惯用胞）

# 4. 查阅文献的计算细节
# 尽量使用相同的计算设置
```

### 6.3.4 计算失败

**问题**：部分或全部变形结构计算失败

**可能原因**：
1. 应变过大，结构不合理
2. 内存不足
3. SCF不收敛

**解决方法**：
```bash
# 1. 减小应变
abacustest model elastic prepare -j Si-elastic --norm 0.005 --shear 0.005

# 2. 增加SCF迭代次数
# 在INPUT中设置：
scf_nmax  200  # 从100增加到200

# 3. 调整mixing参数
mixing_beta  0.4  # 从0.7降低到0.4

# 4. 检查内存使用
# 减少并行进程数或增加节点内存
```

## 6.4 收敛性测试

为了确保结果可靠，建议进行以下收敛性测试：

### 6.4.1 K点收敛

测试不同K点密度对弹性常数的影响：

| K点网格 | C₁₁ (GPa) | C₁₂ (GPa) | C₄₄ (GPa) | 计算时间 |
|---------|----------|----------|----------|---------|
| 6×6×6 | 154.2 | 57.8 | 75.5 | 1× |
| 8×8×8 | 155.5 | 58.3 | 76.2 | 2× |
| 10×10×10 | 155.6 | 58.4 | 76.3 | 4× |
| 12×12×12 | 155.6 | 58.4 | 76.3 | 8× |

当弹性常数变化<1%时，认为K点已收敛。

### 6.4.2 应变大小收敛

测试不同应变大小对结果的影响：

| 应变大小 | C₁₁ (GPa) | C₁₂ (GPa) | C₄₄ (GPa) |
|---------|----------|----------|----------|
| 0.5% | 155.8 | 58.5 | 76.4 |
| 1.0% | 155.5 | 58.3 | 76.2 |
| 1.5% | 155.2 | 58.1 | 76.0 |
| 2.0% | 154.8 | 57.8 | 75.7 |

应变过小（<0.5%）噪声大，过大（>2%）可能超出线弹性范围。1%是合理的默认值。

## 6.5 结果分析总结

**检查清单**：
- [ ] 对称性满足晶系要求（误差<1%）
- [ ] 所有弹性常数为正值
- [ ] 数值在合理范围内
- [ ] 与文献或数据库对比（误差<20%）
- [ ] K点和应变大小已收敛

**如果所有检查通过**：
- 结果可靠，可用于后续分析
- 保存metrics_elastic.json文件
- 记录计算设置（赝势、基组、K点、应变）

**如果检查未通过**：
- 参考6.3节排查问题
- 调整计算设置后重新计算
- 必要时咨询文献或专家

下一章将提供关键步骤的速查表和扩展方向。
# 第7章：总结与速查

本章提供弹性常数计算的关键步骤速查表和扩展方向，便于快速查阅。

## 7.1 完整流程速查

### 7.1.1 立方晶系（如Si）

| 步骤 | 命令 | 关键参数 | 说明 |
|------|------|---------|------|
| 1. 生成结构 | `ase` | cubic=True | 使用惯用胞 |
| 2. 准备优化 | `abacustest model inputs` | --jtype cell-relax | 生成INPUT等文件 |
| 3. 结构优化 | `mpirun -np 4 abacus` | - | 优化晶胞和原子位置 |
| 4. 提取结构 | `cp OUT.*/STRU_ION_D` | - | 获取优化后结构 |
| 5. 准备SCF | 修改INPUT | calculation=scf | 改为SCF计算 |
| 6. 准备弹性 | `abacustest model elastic prepare` | --norm 0.01 | 生成24个变形结构 |
| 7. 提交计算 | 批量运行abacus | - | 计算所有变形结构 |
| 8. 后处理 | `abacustest model elastic post` | - | 输出弹性张量和模量 |
| 9. 验证结果 | 检查对称性 | C₁₁=C₂₂=C₃₃ | 确保满足立方对称性 |

### 7.1.2 四方晶系（如TiO₂）

流程与立方晶系相同，但验证时检查：
- C₁₁ = C₂₂（面内相等）
- C₁₃ = C₂₃（面内-面外耦合相等）
- C₄₄ = C₅₅（面外剪切相等）
- C₃₃ ≠ C₁₁（面外与面内不同）

## 7.2 参数速查表

### 7.2.1 INPUT参数

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| calculation | cell-relax（优化）<br>scf（弹性计算） | 任务类型 |
| ecutwfc | 100 Ry | 平面波截断能 |
| scf_thr | 1e-6 eV | 能量收敛标准 |
| force_thr_ev | 0.01 eV/Å | 力收敛标准 |
| stress_thr | 0.5 kBar | 应力收敛标准 |
| cal_stress | 1 | 开启应力计算 |
| cal_force | 1 | 开启力计算 |
| basis_type | lcao | 使用LCAO基组 |
| gamma_only | 0 | 不使用Gamma点近似 |

### 7.2.2 abacustest参数

| 命令 | 参数 | 推荐值 | 说明 |
|------|------|--------|------|
| model inputs | --jtype | cell-relax | 结构优化 |
| model inputs | --lcao | 布尔标志 | 使用LCAO基组 |
| elastic prepare | --norm | 0.01 | 正应变大小（1%） |
| elastic prepare | --shear | 0.01 | 剪切应变大小（1%） |
| elastic prepare | --norelax | 不推荐 | 不优化原子位置 |

### 7.2.3 K点设置

| 体系大小 | K点网格 | 说明 |
|---------|---------|------|
| 小体系（<10原子） | 12×12×12 | 需要密集K点 |
| 中等体系（10-50原子） | 8×8×8 | 平衡精度和效率 |
| 大体系（>50原子） | 4×4×4 | 减少计算量 |
| 四方晶系 | kx=ky≠kz | 根据晶格常数比例调整 |

## 7.3 不同晶系的注意事项

| 晶系 | 独立分量 | 对称性检查 | 特殊注意 |
|------|---------|-----------|---------|
| 立方 | 3 | C₁₁=C₂₂=C₃₃<br>C₁₂=C₁₃=C₂₃<br>C₄₄=C₅₅=C₆₆ | 最简单，各向同性 |
| 六方 | 5 | C₁₁=C₂₂<br>C₁₃=C₂₃<br>C₄₄=C₅₅ | c轴沿z方向 |
| 四方 | 6 | C₁₁=C₂₂<br>C₁₃=C₂₃<br>C₄₄=C₅₅ | c轴沿z方向<br>C₃₃≠C₁₁ |
| 正交 | 9 | 三个轴互相垂直 | 无额外对称性 |
| 单斜 | 13 | 一个斜角 | 取向复杂 |
| 三斜 | 21 | 无对称性 | 最复杂 |

## 7.4 常见问题速查

| 问题 | 可能原因 | 解决方法 |
|------|---------|---------|
| 对称性不符 | K点不足<br>结构未优化<br>应变过大 | 增加K点<br>检查优化收敛<br>减小应变 |
| 数值不稳定 | 应力精度不足<br>应变过小 | 降低scf_thr<br>增大应变 |
| 与文献差异大 | 赝势不同<br>结构不同<br>温度效应 | 使用相同赝势<br>检查结构<br>查阅文献 |
| 计算失败 | 应变过大<br>SCF不收敛 | 减小应变<br>调整mixing参数 |
| 负弹性常数 | 结构不稳定<br>计算错误 | 检查结构<br>重新计算 |

## 7.5 扩展方向

### 7.5.1 其他晶系

本教程介绍了立方和四方晶系，其他晶系的计算流程类似：

**六方晶系（如Graphite、ZnO）**：
- 5个独立弹性常数
- 确保c轴沿z方向
- K点设置：kx=ky≠kz

**正交晶系（如α-U）**：
- 9个独立弹性常数
- 三个轴互相垂直但长度不同
- K点设置：kx≠ky≠kz

**更复杂的晶系**：
- 单斜（13个）和三斜（21个）
- 需要更仔细的晶体取向检查
- 建议与文献对比验证

### 7.5.2 温度和压力效应

**温度效应**：
- DFT计算是0K结果
- 有限温度需要考虑声子贡献
- 可使用准谐近似（QHA）或自洽声子方法

**压力效应**：
- 在不同压力下优化结构
- 计算每个压力下的弹性常数
- 研究压力诱导的相变

### 7.5.3 高级分析

**声速计算**：
- 从弹性常数计算纵波和横波速度
- 研究声学各向异性

**德拜温度**：
- 从弹性常数估算德拜温度
- 关联热力学性质

**力学稳定性**：
- 检查Born稳定性判据
- 判断晶体结构是否稳定

## 7.6 参考资料

### 7.6.1 软件文档

- **ABACUS官方文档**：http://abacus.deepmodeling.com/
- **abacustest文档**：https://github.com/deepmodeling/abacustest
- **Materials Project**：https://materialsproject.org/

### 7.6.2 理论参考

1. M. de Jong et al., "Charting the complete elastic properties of inorganic crystalline compounds", *Scientific Data* **2**, 150009 (2015)
   - Materials Project的弹性常数计算方法

2. S. Singh et al., "MechElastic: A Python library for analysis of mechanical and elastic properties", *Computer Physics Communications* **267**, 108068 (2021)
   - 弹性常数的后处理和分析

3. D. G. Isaak et al., "Elasticity of TiO₂ rutile to 1800 K", *Physics and Chemistry of Minerals* **26**, 31 (1998)
   - TiO₂弹性常数的实验测量

### 7.6.3 教程和案例

- **ABACUS用户指南**：https://mcresearch.github.io/abacus-user-guide/
- **ABACUS+pymatgen计算弹性常数**：https://mcresearch.github.io/abacus-user-guide/abacus-elastic.html

## 7.7 总结

使用abacustest计算弹性常数的关键要点：

1. **结构准备**：使用惯用胞，确保晶体取向正确
2. **结构优化**：充分优化至力和应力收敛
3. **参数设置**：K点密度、能量收敛、应变大小
4. **结果验证**：对称性检查、数值合理性、与文献对比
5. **问题排查**：参考6.3节的常见问题解决方法

abacustest大大简化了弹性常数计算的流程，但仍需要仔细设置参数和验证结果。希望本教程能帮助你顺利完成弹性常数计算。
