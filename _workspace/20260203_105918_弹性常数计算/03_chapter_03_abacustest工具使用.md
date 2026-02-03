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
