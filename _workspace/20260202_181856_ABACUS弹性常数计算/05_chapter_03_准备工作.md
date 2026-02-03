# 第三章：准备工作

在开始实际计算之前，我们需要安装和配置 abacustest 工具，并了解其基本使用方法。本章将介绍 abacustest 的功能特点、环境配置以及工作流程概览。

## 3.1 abacustest 简介

### 3.1.1 什么是 abacustest

**abacustest** 是 ABACUS 的辅助软件，专门用于 ABACUS 的前后处理以及高通量任务的计算。它提供了一系列自动化工作流，可以大大简化 ABACUS 计算的准备和分析过程。

abacustest 的主要功能包括：
- **输入文件准备**：从结构文件（CIF、POSCAR 等）自动生成 ABACUS 输入文件（INPUT、STRU、KPT）
- **自动化工作流**：内置多种计算工作流，如弹性常数计算、声子谱计算等
- **结果后处理**：自动提取和分析计算结果，生成可视化数据

对于弹性常数计算，abacustest 提供了完整的工作流：
1. 自动生成变形结构（24 个变形结构 + 1 个原始结构）
2. 批量提交 ABACUS 计算
3. 自动拟合应力-应变关系，输出弹性张量和弹性模量

### 3.1.2 abacustest vs ABACUS+pymatgen

除了 abacustest，还有另一种常用的方法是使用 **ABACUS+pymatgen**。两种方法的对比如下：

| 特性 | abacustest | ABACUS+pymatgen |
|-----|-----------|----------------|
| 自动化程度 | 高（一键式命令） | 中等（需要编写脚本） |
| 输入文件准备 | 自动生成 | 需要手动准备 |
| 变形结构生成 | 自动 | 需要运行 Python 脚本 |
| 结果后处理 | 自动 | 需要运行 Python 脚本 |
| 灵活性 | 中等 | 高（可自定义） |
| 学习曲线 | 平缓 | 较陡 |

对于初学者和希望快速完成计算的用户，**推荐使用 abacustest**。

### 3.1.3 核心命令

abacustest 的弹性常数计算主要使用三个命令：

1. **`abacustest model inputs`**：从结构文件准备 ABACUS 输入文件
2. **`abacustest model elastic prepare`**：生成弹性常数计算所需的变形结构
3. **`abacustest model elastic post`**：后处理计算结果，输出弹性张量

这三个命令构成了完整的计算流程。

## 3.2 环境配置

### 3.2.1 安装 abacustest

abacustest 是一个 Python 包，可以通过 pip 从 PyPI 源安装：

```bash
pip install abacustest
```

安装完成后，可以通过以下命令验证安装：

```bash
abacustest --version
```

如果显示版本号，说明安装成功。

> **注意**：abacustest 需要 Python 3.6 或更高版本。

### 3.2.2 配置环境变量

abacustest 在准备输入文件时，会自动使用环境变量中指定的赝势和数值原子轨道文件。需要配置以下两个环境变量：

```bash
export ABACUS_PP_PATH=/path/to/pseudopotentials
export ABACUS_ORB_PATH=/path/to/orbitals
```

**环境变量说明**：
- **ABACUS_PP_PATH**：赝势文件所在目录的路径
- **ABACUS_ORB_PATH**：数值原子轨道文件所在目录的路径

建议将这两行添加到 `~/.bashrc` 或 `~/.bash_profile` 文件中，以便每次登录时自动加载。

### 3.2.3 准备赝势和轨道文件

对于 LCAO 基组计算，需要准备两类文件：

**1. 赝势文件（Pseudopotential）**

赝势文件描述了原子核和内层电子对价电子的有效势。常用的赝势库包括：
- **SG15**：Optimized norm-conserving Vanderbilt (ONCV) 赝势
- **DOJO**：PseudoDojo 赝势库
- **APNS-v1**：ABACUS Pseudopotential and Numerical Orbital Set v1（推荐）

赝势文件通常以 `.upf` 为扩展名，例如 `Si_ONCV_PBE-1.0.upf`。

**2. 数值原子轨道文件（Numerical Atomic Orbital）**

数值原子轨道文件定义了 LCAO 基组。ABACUS 提供了不同精度的基组：
- **minimal**：最小基组，计算最快但精度较低
- **efficiency**：效率基组，平衡精度和速度（推荐）
- **precision**：精度基组，精度最高但计算较慢

轨道文件通常以 `.orb` 为扩展名，例如 `Si_gga_7au_100Ry_2s2p1d.orb`。

**推荐配置**：
- 赝势：APNS-v1
- 基组：efficiency

这些文件可以从 [ABACUS 官方网站](https://abacus.deepmodeling.com/) 或 [Bohrium 平台](https://bohrium.dp.tech/) 下载。

### 3.2.4 验证环境配置

配置完成后，可以通过以下命令验证环境变量是否正确设置：

```bash
echo $ABACUS_PP_PATH
echo $ABACUS_ORB_PATH
```

应该显示你设置的路径。同时，检查这些路径下是否包含赝势和轨道文件：

```bash
ls $ABACUS_PP_PATH
ls $ABACUS_ORB_PATH
```

## 3.3 工作流程概览

### 3.3.1 完整计算流程

使用 abacustest 计算弹性常数的完整流程如下：

```
结构准备 → 结构优化 → 弹性计算 → 结果分析
   ↓           ↓           ↓           ↓
 CIF文件    cell-relax   变形结构    弹性张量
            优化后STRU   应力计算    弹性模量
```

**详细步骤**：

1. **结构准备**
   - 获取或生成晶体结构文件（CIF、POSCAR 等）
   - 使用 `abacustest model inputs` 生成 ABACUS 输入文件

2. **结构优化**
   - 设置 `calculation = cell-relax`
   - 运行 ABACUS，优化晶胞和原子位置
   - 提取优化后的结构（`STRU_ION_D`）

3. **弹性计算准备**
   - 将优化后的结构复制到新文件夹
   - 修改 INPUT 为 `calculation = scf`
   - 使用 `abacustest model elastic prepare` 生成变形结构

4. **应力计算**
   - 对所有变形结构运行 ABACUS SCF 计算
   - 确保 `cal_stress = 1` 开启应力计算

5. **结果后处理**
   - 使用 `abacustest model elastic post` 提取应力数据
   - 自动拟合应力-应变关系
   - 输出弹性张量和弹性模量

### 3.3.2 变形结构的生成

`abacustest model elastic prepare` 命令会自动生成以下文件夹：

```
job_folder/
├── org/                    # 原始结构（无应变）
├── deform-xx-1/            # x方向正应变，大小1
├── deform-xx-2/            # x方向正应变，大小2
├── deform-xx-3/            # x方向负应变，大小1
├── deform-xx-4/            # x方向负应变，大小2
├── deform-yy-1/            # y方向正应变
├── deform-yy-2/
├── ...
├── deform-zz-1/            # z方向正应变
├── ...
├── deform-yz-1/            # yz剪切应变
├── ...
├── deform-xz-1/            # xz剪切应变
├── ...
├── deform-xy-1/            # xy剪切应变
└── ...
```

共 **25 个文件夹**：
- 1 个原始结构（org）
- 24 个变形结构（6 种应变类型 × 4 种应变大小）

每个文件夹包含完整的 ABACUS 输入文件（INPUT、STRU、KPT、赝势、轨道）。

### 3.3.3 应变大小的设置

默认情况下，abacustest 对每种应变类型施加 4 种不同大小的应变：

- 正应变：+0.5%, +1.0%, -0.5%, -1.0%
- 剪切应变：+0.5%, +1.0%, -0.5%, -1.0%

可以通过参数调整应变大小：

```bash
abacustest model elastic prepare -j job_folder --norm 0.01 --shear 0.01
```

参数说明：
- `--norm`：最大正应变（默认 0.01，即 ±1%）
- `--shear`：最大剪切应变（默认 0.01，即 ±1%）

> **注意**：应变大小应该在线弹性范围内，通常 ±1% 是合适的。

### 3.3.4 原子弛豫的选择

默认情况下，abacustest 会对每个变形结构进行固定晶胞的结构优化（`calculation = relax`），允许原子位置调整以达到力的平衡。

如果希望跳过原子弛豫，直接计算应力，可以使用 `--norelax` 参数：

```bash
abacustest model elastic prepare -j job_folder --norelax
```

此时，INPUT 文件中的 `calculation` 会被设置为 `scf`，原子位置固定。

**建议**：
- 对于大多数材料，**推荐允许原子弛豫**，这样得到的弹性常数更接近实验值
- 对于高对称性的简单结构（如 Si），原子弛豫的影响较小，可以考虑使用 `--norelax` 加快计算

### 3.3.5 后处理输出

`abacustest model elastic post` 命令会输出以下信息：

**1. 弹性张量矩阵**（6×6 矩阵，单位 GPa）

```
elastic_tensor:
         0          1          2       3       4       5
0  155.465  58.326  58.326   0.000   0.000   0.000
1   58.326 155.465  58.326   0.000   0.000   0.000
2   58.326  58.326 155.465   0.000   0.000   0.000
3    0.000   0.000   0.000  76.177   0.000   0.000
4    0.000   0.000   0.000   0.000  76.177   0.000
5    0.000   0.000   0.000   0.000   0.000  76.177
```

**2. 弹性模量**（单位 GPa，泊松比无量纲）

```
             bulk_modulus  shear_modulus  young_modulus  poisson_ratio
job_folder/     90.706        65.134        157.664          0.210
```

**3. 结果文件**

- `metrics.json`：包含所有计算指标
- `metrics_elastic.json`：包含弹性常数的详细数据

这些文件可以用于进一步分析或与其他数据对比。

## 3.4 命令参数详解

### 3.4.1 `abacustest model inputs` 命令

**功能**：从结构文件生成 ABACUS 输入文件

**基本用法**：
```bash
abacustest model inputs -f structure.cif --ftype cif --jtype cell-relax --lcao --folder-syntax folder_name
```

**主要参数**：

| 参数 | 说明 | 示例 |
|-----|------|------|
| `-f` | 结构文件名（可指定多个） | `Si.cif` |
| `--ftype` | 结构文件格式 | `cif`, `vasp`, `poscar` |
| `--jtype` | ABACUS 任务类型 | `scf`, `relax`, `cell-relax`, `nscf` |
| `--lcao` | 使用 LCAO 基组 | 无需参数值 |
| `--pw` | 使用平面波基组 | 无需参数值 |
| `--folder-syntax` | 输出文件夹名称 | `Si-calc` |

**任务类型说明**：
- `scf`：自洽场计算（固定结构）
- `relax`：结构优化（固定晶胞，优化原子位置）
- `cell-relax`：完全优化（优化晶胞和原子位置）
- `nscf`：非自洽计算（用于能带、态密度）

### 3.4.2 `abacustest model elastic prepare` 命令

**功能**：生成弹性常数计算所需的变形结构

**基本用法**：
```bash
abacustest model elastic prepare -j job_folder
```

**主要参数**：

| 参数 | 说明 | 默认值 |
|-----|------|-------|
| `-j` | 输入文件夹路径 | 必需 |
| `--norm` | 最大正应变 | 0.01 |
| `--shear` | 最大剪切应变 | 0.01 |
| `--norelax` | 不对变形结构做优化 | False |

> **重要提示**：不要重复执行此命令！重复执行会删除已有的 `deform-*` 文件夹和计算结果。

### 3.4.3 `abacustest model elastic post` 命令

**功能**：后处理弹性常数计算结果

**基本用法**：
```bash
abacustest model elastic post -j job_folder
```

**主要参数**：

| 参数 | 说明 |
|-----|------|
| `-j` | 计算文件夹路径（包含 deform-* 和 org 文件夹） |

此命令会自动：
1. 从每个 `deform-*` 文件夹中提取应力数据
2. 拟合应力-应变关系
3. 计算弹性张量矩阵
4. 计算弹性模量（体模量、剪切模量、杨氏模量、泊松比）
5. 将结果保存到 `metrics.json` 和 `metrics_elastic.json`

## 3.5 注意事项

### 3.5.1 INPUT 参数的调整

abacustest 自动生成的 INPUT 文件中的参数适配 ABACUS v3.10，在很多情况下是合理的。但对于特定体系，可能需要手动调整以下参数：

**必须检查的参数**：
- `nspin`：自旋极化（磁性材料需要设置为 2 或 4）
- `cal_stress`：应力计算（弹性常数计算必须设置为 1）
- `cal_force`：力计算（结构优化需要设置为 1）

**可能需要调整的参数**：
- `ecutwfc`：平面波截断能（影响精度）
- `scf_thr`：SCF 收敛阈值（影响精度）
- `smearing_method` 和 `smearing_sigma`：展宽方法和参数
- `mixing_type` 和 `mixing_beta`：电荷密度混合方法

### 3.5.2 k 点设置

k 点网格的密度直接影响计算精度。abacustest 会根据晶胞大小自动生成 KPT 文件，但建议进行 k 点收敛性测试。

对于弹性常数计算，通常需要较密的 k 点网格，以确保应力计算的精度。

### 3.5.3 计算资源

弹性常数计算需要运行 25 个 ABACUS 任务（1 个原始结构 + 24 个变形结构）。建议：
- 使用并行计算（MPI）加快单个任务
- 使用任务调度系统（如 Slurm、PBS）批量提交任务
- 合理分配计算资源，避免浪费

### 3.5.4 结果验证

计算完成后，应该验证结果的合理性：
- 检查弹性张量的对称性（Cᵢⱼ = Cⱼᵢ）
- 检查是否符合晶体对称性（例如立方晶系应该只有 3 个独立分量）
- 与已知数据对比（Materials Project、实验值）
- 检查弹性模量的合理性（体模量、剪切模量应为正值）

---

## 本章小结

本章介绍了 abacustest 工具的使用准备：

- **abacustest 简介**：ABACUS 的辅助软件，提供自动化工作流
- **环境配置**：安装 abacustest，配置环境变量，准备赝势和轨道文件
- **工作流程**：结构准备 → 优化 → 弹性计算 → 后处理
- **核心命令**：`inputs`（准备输入）、`elastic prepare`（生成变形结构）、`elastic post`（后处理）
- **注意事项**：参数调整、k 点设置、结果验证

在下一章中，我们将通过第一个案例（Si）演示完整的计算流程。
