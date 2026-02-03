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
