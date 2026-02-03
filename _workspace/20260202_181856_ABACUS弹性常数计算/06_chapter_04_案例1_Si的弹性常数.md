# 第四章：案例1 - 计算Si的弹性常数

本章将通过一个完整的案例，演示如何使用 abacustest 计算立方晶系材料 Si（硅）的弹性常数。Si 是最常见的半导体材料，具有金刚石结构，属于立方晶系，只有 3 个独立的弹性常数。

## 4.1 Si 的晶体结构

### 4.1.1 晶体学信息

**硅（Silicon, Si）** 的基本信息：
- **晶系**：立方晶系（Cubic）
- **空间群**：Fd-3m (227)
- **晶体结构**：金刚石结构（Diamond structure）
- **晶格常数**：a ≈ 5.43 Å
- **原子数**：惯用胞包含 8 个原子

### 4.1.2 弹性常数的对称性

由于 Si 属于立方晶系，其弹性张量矩阵具有高度对称性，只有 **3 个独立的弹性常数**：

$$
C = \begin{bmatrix}
C_{11} & C_{12} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{11} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{12} & C_{11} & 0 & 0 & 0 \\
0 & 0 & 0 & C_{44} & 0 & 0 \\
0 & 0 & 0 & 0 & C_{44} & 0 \\
0 & 0 & 0 & 0 & 0 & C_{44}
\end{bmatrix}
$$

- **C₁₁**：沿主轴方向的刚度
- **C₁₂**：泊松效应
- **C₄₄**：剪切刚度

### 4.1.3 已知数据

根据 Materials Project 数据库，Si 的弹性常数实验值和计算值为：
- C₁₁ ≈ 153 GPa
- C₁₂ ≈ 57 GPa
- C₄₄ ≈ 74 GPa

计算结果应该与这些值接近。

## 4.2 生成 Si 的晶体结构

### 4.2.1 使用 ASE 生成结构文件

首先需要生成 Si 晶体的结构文件。可以使用 Python 的 ASE（Atomic Simulation Environment）库生成惯用胞的 CIF 文件：

```python
from ase.build import bulk

# 生成Si的立方惯用胞
si_conv = bulk('Si', cubic=True)

# 保存为CIF文件
si_conv.write("Si_conv.cif")
```

**代码说明**：
- `bulk('Si', cubic=True)`：生成 Si 的立方惯用胞（8 个原子）
- `cubic=True`：确保晶格矢量沿坐标轴方向，无需后续旋转
- `write("Si_conv.cif")`：保存为 CIF 格式

执行此脚本后，会在当前目录生成 `Si_conv.cif` 文件。

> **注意**：如果没有安装 ASE，可以使用 `pip install ase` 安装。

### 4.2.2 检查结构文件

可以使用文本编辑器查看生成的 CIF 文件，或使用可视化软件（如 VESTA）查看结构。

Si 惯用胞的特点：
- 晶格矢量沿 x、y、z 轴方向
- 晶格常数 a = b = c ≈ 5.43 Å
- 包含 8 个 Si 原子

## 4.3 结构优化

### 4.3.1 准备输入文件

使用 abacustest 从 CIF 文件生成 ABACUS 输入文件：

```bash
abacustest model inputs -f Si_conv.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si
```

**命令说明**：
- `-f Si_conv.cif`：指定结构文件
- `--ftype cif`：文件格式为 CIF
- `--jtype cell-relax`：任务类型为完全优化（优化晶胞和原子位置）
- `--lcao`：使用 LCAO 基组
- `--folder-syntax Si`：输出文件夹名称为 `Si`

执行后，会在当前目录生成 `Si/` 文件夹，包含以下文件：
```
Si/
├── INPUT           # ABACUS 主输入文件
├── STRU            # 结构文件
├── KPT             # k点文件
├── Si_ONCV_PBE-1.0.upf    # Si的赝势文件（自动复制）
└── Si_gga_7au_100Ry_2s2p1d.orb  # Si的轨道文件（自动复制）
```

### 4.3.2 检查 INPUT 参数

打开 `Si/INPUT` 文件，检查关键参数：

```
INPUT_PARAMETERS
calculation     cell-relax    # 完全优化
basis_type      lcao          # LCAO基组
ecutwfc         100           # 截断能（Ry）
scf_thr         1e-7          # SCF收敛阈值
cal_force       1             # 计算力
cal_stress      1             # 计算应力（重要！）
relax_nmax      100           # 最大优化步数
force_thr_ev    0.01          # 力收敛阈值（eV/Å）
stress_thr      0.5           # 应力收敛阈值（kBar）
smearing_method gaussian      # 展宽方法
smearing_sigma  0.002         # 展宽参数
mixing_type     pulay         # 混合方法
mixing_beta     0.4           # 混合参数
```

**重要参数**：
- `cal_stress = 1`：必须开启应力计算，否则无法计算弹性常数
- `cal_force = 1`：结构优化需要计算力
- `scf_thr`：SCF 收敛阈值，影响应力计算精度
- `force_thr_ev` 和 `stress_thr`：优化收敛标准

如果需要，可以手动调整这些参数。对于 Si，默认参数通常已经足够。

### 4.3.3 提交优化计算

进入 `Si/` 文件夹，提交 ABACUS 计算：

```bash
cd Si
mpirun -np 4 abacus
```

**说明**：
- `mpirun -np 4`：使用 4 个 MPI 进程并行计算
- 根据你的计算环境，可能需要使用任务调度系统（如 Slurm）提交任务

计算完成后，会在 `Si/OUT.ABACUS/` 目录生成输出文件，其中：
- `running_cell-relax.log`：计算日志
- `STRU_ION_D`：优化后的结构文件

### 4.3.4 检查优化结果

查看 `running_cell-relax.log` 文件的最后部分，确认优化已收敛：

```
RELAX IONS : max force is 0.005 eV/Å
RELAX CELL : max stress is 0.3 kBar
Relaxation converged after 15 steps
```

如果看到类似信息，说明优化成功。

### 4.3.5 提取优化后的结构

优化完成后，需要提取优化后的结构，并准备弹性常数计算的输入文件：

```bash
# 返回上级目录
cd ..

# 创建弹性计算文件夹
mkdir Si-elastic

# 复制INPUT文件和赝势/轨道文件
cp Si/INPUT Si/Si* Si-elastic/

# 复制优化后的结构
cp Si/OUT.ABACUS/STRU_ION_D Si-elastic/STRU
```

**说明**：
- `Si/Si*`：复制所有以 `Si` 开头的文件（赝势和轨道文件）
- `STRU_ION_D`：优化后的结构文件，包含优化后的晶格常数和原子位置

### 4.3.6 修改 INPUT 为 SCF 计算

编辑 `Si-elastic/INPUT` 文件，将 `calculation` 改为 `scf`：

```
calculation     scf    # 改为SCF计算
```

其他参数保持不变。确保 `cal_stress = 1` 仍然开启。

## 4.4 弹性常数计算

### 4.4.1 准备变形结构

使用 abacustest 生成弹性常数计算所需的变形结构：

```bash
abacustest model elastic prepare -j Si-elastic
```

**命令说明**：
- `-j Si-elastic`：指定输入文件夹

执行后，会在 `Si-elastic/` 目录下生成一系列文件夹：

```
Si-elastic/
├── INPUT
├── STRU
├── KPT
├── Si_ONCV_PBE-1.0.upf
├── Si_gga_7au_100Ry_2s2p1d.orb
├── org/                    # 原始结构（无应变）
├── deform-xx-1/            # x方向正应变 +0.5%
├── deform-xx-2/            # x方向正应变 +1.0%
├── deform-xx-3/            # x方向负应变 -0.5%
├── deform-xx-4/            # x方向负应变 -1.0%
├── deform-yy-1/            # y方向正应变
├── deform-yy-2/
├── deform-yy-3/
├── deform-yy-4/
├── deform-zz-1/            # z方向正应变
├── deform-zz-2/
├── deform-zz-3/
├── deform-zz-4/
├── deform-yz-1/            # yz剪切应变
├── deform-yz-2/
├── deform-yz-3/
├── deform-yz-4/
├── deform-xz-1/            # xz剪切应变
├── deform-xz-2/
├── deform-xz-3/
├── deform-xz-4/
├── deform-xy-1/            # xy剪切应变
├── deform-xy-2/
├── deform-xy-3/
└── deform-xy-4/
```

共 **25 个文件夹**（1 个原始 + 24 个变形）。

> **重要提示**：不要重复执行 `elastic prepare` 命令！重复执行会删除已有的文件夹和计算结果。

### 4.4.2 应变设置说明

默认情况下，abacustest 对每种应变类型施加 4 种不同大小的应变：

| 应变类型 | 应变值 |
|---------|-------|
| 正应变（xx, yy, zz） | +0.5%, +1.0%, -0.5%, -1.0% |
| 剪切应变（yz, xz, xy） | +0.5%, +1.0%, -0.5%, -1.0% |

这些应变值在线弹性范围内，适合大多数材料。

### 4.4.3 提交计算

需要对所有 25 个文件夹运行 ABACUS 计算。可以使用循环脚本批量提交：

**方法1：串行计算（适合小规模测试）**

```bash
cd Si-elastic

# 计算原始结构
cd org && mpirun -np 4 abacus && cd ..

# 计算所有变形结构
for dir in deform-*; do
    cd $dir
    mpirun -np 4 abacus
    cd ..
done
```

**方法2：并行提交（推荐，适合集群环境）**

如果使用 Slurm 任务调度系统，可以创建提交脚本 `submit_all.sh`：

```bash
#!/bin/bash

cd Si-elastic

# 提交原始结构
cd org
sbatch run.slurm
cd ..

# 提交所有变形结构
for dir in deform-*; do
    cd $dir
    sbatch run.slurm
    cd ..
done
```

其中 `run.slurm` 是单个任务的 Slurm 脚本。

**计算时间估计**：
- 单个任务（8 原子 Si）：约 1-5 分钟（取决于硬件）
- 总计算时间：约 25-125 分钟（如果串行）
- 如果并行提交，总时间可以大大缩短

### 4.4.4 检查计算状态

可以使用以下命令检查计算是否完成：

```bash
# 检查所有文件夹中是否有输出文件
ls Si-elastic/*/OUT.ABACUS/running_*.log

# 检查是否所有计算都正常结束
grep "Total Time" Si-elastic/*/OUT.ABACUS/running_*.log
```

如果所有文件夹都有 `running_*.log` 文件，并且包含 "Total Time"，说明计算已完成。

## 4.5 结果分析

### 4.5.1 后处理

所有计算完成后，使用 abacustest 进行后处理：

```bash
abacustest model elastic post -j Si-elastic
```

**命令说明**：
- `-j Si-elastic`：指定计算文件夹

此命令会自动：
1. 从每个 `deform-*` 和 `org` 文件夹中提取应力数据
2. 拟合应力-应变关系
3. 计算弹性张量矩阵
4. 计算弹性模量

### 4.5.2 输出结果

后处理完成后，会在屏幕上输出结果：

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

### 4.5.3 弹性张量分析

从输出的弹性张量矩阵可以看出：

**对角元素**：
- C₁₁ = C₂₂ = C₃₃ = **155.465 GPa**
- C₄₄ = C₅₅ = C₆₆ = **76.177 GPa**

**非对角元素**：
- C₁₂ = C₁₃ = C₂₃ = **58.326 GPa**
- 其他元素 ≈ 0

这完全符合立方晶系的对称性！只有 3 个独立的弹性常数：
- **C₁₁ = 155.5 GPa**
- **C₁₂ = 58.3 GPa**
- **C₄₄ = 76.2 GPa**

### 4.5.4 弹性模量分析

计算得到的弹性模量：

| 弹性模量 | 数值 | 单位 |
|---------|------|------|
| 体模量（Bulk Modulus） | 90.7 | GPa |
| 剪切模量（Shear Modulus） | 65.1 | GPa |
| 杨氏模量（Young's Modulus） | 157.7 | GPa |
| 泊松比（Poisson's Ratio） | 0.210 | 无量纲 |

**物理意义**：
- **体模量 90.7 GPa**：Si 抵抗均匀压缩的能力中等
- **剪切模量 65.1 GPa**：Si 抵抗剪切变形的能力较强
- **杨氏模量 157.7 GPa**：Si 的刚度中等
- **泊松比 0.210**：在典型范围内（0.2-0.3）

### 4.5.5 与 Materials Project 对比

将本文计算结果与 Materials Project 数据库的数据对比：

| 弹性常数 | 本文计算 | Materials Project | 差异 |
|---------|---------|-------------------|------|
| C₁₁ (GPa) | 155.5 | 153 | +2.5 (+1.6%) |
| C₁₂ (GPa) | 58.3 | 57 | +1.3 (+2.3%) |
| C₄₄ (GPa) | 76.2 | 74 | +2.2 (+3.0%) |

**结论**：
- 本文计算结果与 Materials Project 的数据接近
- 差异在 1-3% 范围内，属于正常的计算误差
- 这验证了我们的计算方法和参数设置是正确的

差异的可能来源：
- 赝势和基组的选择不同
- k 点网格密度不同
- 截断能设置不同
- 优化收敛标准不同

### 4.5.6 结果文件

后处理完成后，会在 `Si-elastic/` 目录生成两个 JSON 文件：

**1. metrics.json**
包含所有计算指标，可以用于进一步分析或与其他数据对比。

**2. metrics_elastic.json**
包含弹性常数的详细数据，格式如下：

```json
{
  "elastic_tensor": [
    [155.465, 58.326, 58.326, 0.0, 0.0, 0.0],
    [58.326, 155.465, 58.326, 0.0, 0.0, 0.0],
    [58.326, 58.326, 155.465, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 76.177, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 76.177, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 76.177]
  ],
  "bulk_modulus": 90.706,
  "shear_modulus": 65.134,
  "young_modulus": 157.664,
  "poisson_ratio": 0.210
}
```

这些文件可以用于后续的数据分析或可视化。

## 4.6 结果验证

### 4.6.1 对称性检查

检查弹性张量是否满足立方晶系的对称性：

✅ C₁₁ = C₂₂ = C₃₃ = 155.465 GPa（相等）
✅ C₁₂ = C₁₃ = C₂₃ = 58.326 GPa（相等）
✅ C₄₄ = C₅₅ = C₆₆ = 76.177 GPa（相等）
✅ 其他元素 ≈ 0（约10⁻¹⁰数量级）

对称性完全满足！

### 4.6.2 力学稳定性检查

对于立方晶系，力学稳定性的判据为：

1. C₁₁ > 0 ✅（155.5 > 0）
2. C₄₄ > 0 ✅（76.2 > 0）
3. C₁₁ > |C₁₂| ✅（155.5 > 58.3）
4. C₁₁ + 2C₁₂ > 0 ✅（155.5 + 2×58.3 = 272.1 > 0）

所有判据都满足，说明 Si 在计算的结构下是**力学稳定**的。

### 4.6.3 与实验值对比

Si 的弹性常数实验值（室温）：
- C₁₁ ≈ 165.7 GPa
- C₁₂ ≈ 63.9 GPa
- C₄₄ ≈ 79.6 GPa

我们的计算值略低于实验值（约 6-9%），这是正常的，因为：
- DFT 计算通常对应 0 K，而实验值是室温
- 不同的交换关联泛函会影响结果
- 零点振动效应未考虑

总体而言，计算结果与实验值的趋势一致。

## 4.7 本章小结

本章通过 Si 的案例，完整演示了使用 abacustest 计算弹性常数的流程：

1. **结构准备**：使用 ASE 生成 Si 惯用胞的 CIF 文件
2. **结构优化**：使用 `abacustest model inputs` 准备输入文件，运行 `cell-relax` 优化
3. **弹性计算**：使用 `abacustest model elastic prepare` 生成 25 个变形结构，运行 SCF 计算
4. **结果分析**：使用 `abacustest model elastic post` 后处理，得到弹性张量和弹性模量
5. **结果验证**：与 Materials Project 和实验值对比，验证计算的正确性

**关键结果**：
- C₁₁ = 155.5 GPa，C₁₂ = 58.3 GPa，C₄₄ = 76.2 GPa
- 与 Materials Project 数据差异在 1-3% 范围内
- 对称性和力学稳定性检查通过

在下一章中，我们将计算四方晶系材料 TiO₂ 的弹性常数，展示更复杂晶系的计算方法。
