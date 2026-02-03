# 第五章：案例2 - 计算金红石型TiO₂的弹性常数

本章将演示如何计算四方晶系材料金红石型 TiO₂（二氧化钛）的弹性常数。与立方晶系的 Si 相比，TiO₂ 的对称性较低，有 6 个独立的弹性常数，计算流程基本相同，但结果分析更复杂。

## 5.1 TiO₂ 的晶体结构

### 5.1.1 晶体学信息

**金红石型 TiO₂（Rutile TiO₂）** 的基本信息：
- **晶系**：四方晶系（Tetragonal）
- **空间群**：P4₂/mnm (136)
- **Laue 类型**：4/mmm（Laue 类型 I）
- **晶格常数**：a = b ≈ 4.59 Å，c ≈ 2.96 Å
- **原子数**：惯用胞包含 2 个 Ti 原子和 4 个 O 原子
- **Materials Project 编号**：mp-2657

### 5.1.2 弹性常数的对称性

四方晶系（Laue 类型 I）的弹性张量矩阵形式为：

$$
C = \begin{bmatrix}
C_{11} & C_{12} & C_{13} & 0 & 0 & 0 \\
C_{12} & C_{11} & C_{13} & 0 & 0 & 0 \\
C_{13} & C_{13} & C_{33} & 0 & 0 & 0 \\
0 & 0 & 0 & C_{44} & 0 & 0 \\
0 & 0 & 0 & 0 & C_{44} & 0 \\
0 & 0 & 0 & 0 & 0 & C_{66}
\end{bmatrix}
$$

有 **6 个独立的弹性常数**：
- **C₁₁**：垂直于 c 轴方向的刚度（= C₂₂）
- **C₃₃**：沿 c 轴方向的刚度
- **C₁₂**：垂直于 c 轴平面内的泊松效应
- **C₁₃**：c 轴与垂直方向的耦合（= C₂₃）
- **C₄₄**：包含 c 轴的剪切刚度（= C₅₅）
- **C₆₆**：垂直于 c 轴平面内的剪切刚度

### 5.1.3 已知数据

根据文献和 Materials Project 数据库，金红石型 TiO₂ 的弹性常数：

**实验测量值**（D. G. Isaak et al., 1998）：
- C₁₁ = 268.0 GPa
- C₁₂ = 174.9 GPa
- C₁₃ = 147.4 GPa
- C₃₃ = 484.2 GPa
- C₄₄ = 123.8 GPa
- C₆₆ = 190.2 GPa

**Materials Project 计算值**：
- C₁₁ = 426 GPa
- C₁₂ = 2 GPa
- C₁₃ = 149 GPa
- C₃₃ = 470 GPa
- C₄₄ = 113 GPa
- C₆₆ = 43 GPa

可以看到，Materials Project 的某些值（如 C₁₂ 和 C₆₆）与实验值差异较大。

## 5.2 结构准备与优化

### 5.2.1 获取结构文件

金红石型 TiO₂ 的结构可以从 Materials Project 数据库下载：

1. 访问 [Materials Project](https://materialsproject.org/)
2. 搜索 "TiO2 rutile" 或使用编号 "mp-2657"
3. 下载 CIF 文件，保存为 `TiO2.cif`

或者，可以使用 Materials Project API 下载：

```python
from mp_api.client import MPRester

# 需要API key（在Materials Project网站注册获取）
with MPRester("your_api_key") as mpr:
    structure = mpr.get_structure_by_material_id("mp-2657")
    structure.to(filename="TiO2.cif")
```

### 5.2.2 准备输入文件

使用 abacustest 从 CIF 文件生成 ABACUS 输入文件：

```bash
abacustest model inputs -f TiO2.cif --ftype cif --lcao --jtype cell-relax --folder-syntax TiO2-rutile
```

**命令说明**：
- `-f TiO2.cif`：指定结构文件
- `--ftype cif`：文件格式为 CIF
- `--lcao`：使用 LCAO 基组
- `--jtype cell-relax`：任务类型为完全优化
- `--folder-syntax TiO2-rutile`：输出文件夹名称

执行后，会在当前目录生成 `TiO2-rutile/` 文件夹，包含：
```
TiO2-rutile/
├── INPUT
├── STRU
├── KPT
├── Ti_ONCV_PBE-1.0.upf    # Ti的赝势文件
├── Ti_gga_7au_100Ry_2s2p1d.orb  # Ti的轨道文件
├── O_ONCV_PBE-1.0.upf     # O的赝势文件
└── O_gga_7au_100Ry_2s2p1d.orb   # O的轨道文件
```

### 5.2.3 检查和调整 INPUT 参数

打开 `TiO2-rutile/INPUT` 文件，检查关键参数。对于 TiO₂，可能需要调整以下参数：

```
INPUT_PARAMETERS
calculation     cell-relax
basis_type      lcao
ecutwfc         100           # 可能需要增加到120
scf_thr         1e-7
scf_nmax        200           # TiO2可能需要更多SCF步数
cal_force       1
cal_stress      1
relax_nmax      100
force_thr_ev    0.01
stress_thr      0.5
smearing_method gaussian
smearing_sigma  0.002
mixing_type     pulay         # 如果不收敛，可以改为broyden
mixing_beta     0.4           # 如果不收敛，可以减小到0.3
```

**TiO₂ 特殊注意事项**：
- TiO₂ 是氧化物，SCF 收敛可能比 Si 困难
- 如果遇到收敛问题，可以尝试：
  - 增加 `scf_nmax`
  - 减小 `mixing_beta`
  - 改用 `mixing_type broyden`

### 5.2.4 提交优化计算

```bash
cd TiO2-rutile
mpirun -np 4 abacus
cd ..
```

等待计算完成。TiO₂ 的优化可能需要更长时间（10-30 分钟，取决于硬件）。

### 5.2.5 提取优化后的结构

优化完成后，提取优化后的结构并准备弹性计算：

```bash
# 创建弹性计算文件夹
mkdir TiO2-rutile-elastic

# 复制INPUT文件和赝势/轨道文件
cp TiO2-rutile/INPUT TiO2-rutile/Ti* TiO2-rutile/O* TiO2-rutile-elastic/

# 复制优化后的结构
cp TiO2-rutile/OUT.ABACUS/STRU_ION_D TiO2-rutile-elastic/STRU
```

**注意**：这里使用 `Ti*` 和 `O*` 复制所有 Ti 和 O 的赝势/轨道文件。

### 5.2.6 修改 INPUT 为 SCF 计算

编辑 `TiO2-rutile-elastic/INPUT` 文件：

```
calculation     scf    # 改为SCF计算
```

确保 `cal_stress = 1` 仍然开启。

## 5.3 弹性常数计算

### 5.3.1 准备变形结构

使用 abacustest 生成变形结构：

```bash
abacustest model elastic prepare -j TiO2-rutile-elastic
```

与 Si 案例相同，会生成 25 个文件夹（1 个原始 + 24 个变形）。

### 5.3.2 提交计算

对所有 25 个文件夹运行 ABACUS 计算：

```bash
cd TiO2-rutile-elastic

# 计算原始结构
cd org && mpirun -np 4 abacus && cd ..

# 计算所有变形结构
for dir in deform-*; do
    cd $dir
    mpirun -np 4 abacus
    cd ..
done
```

或者使用任务调度系统批量提交。

**计算时间估计**：
- 单个任务（6 原子 TiO₂）：约 2-10 分钟
- 总计算时间：约 50-250 分钟（如果串行）

### 5.3.3 检查计算状态

```bash
# 检查是否所有计算都完成
ls TiO2-rutile-elastic/*/OUT.ABACUS/running_*.log

# 检查是否正常结束
grep "Total Time" TiO2-rutile-elastic/*/OUT.ABACUS/running_*.log
```

## 5.4 结果分析

### 5.4.1 后处理

所有计算完成后，进行后处理：

```bash
abacustest model elastic post -j TiO2-rutile-elastic
```

### 5.4.2 输出结果

后处理完成后，会在屏幕上输出结果：

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

### 5.4.3 弹性张量分析

从输出的弹性张量矩阵可以提取 6 个独立的弹性常数：

**对角元素**：
- C₁₁ = C₂₂ = **281.2 GPa**
- C₃₃ = **484.0 GPa**
- C₄₄ = C₅₅ = **117.4 GPa**
- C₆₆ = **216.4 GPa**

**非对角元素**：
- C₁₂ = **158.1 GPa**
- C₁₃ = C₂₃ = **155.9 GPa**
- 其他元素 ≈ 0（数量级 10⁻⁸ 或更小）

这完全符合四方晶系（Laue 类型 I）的对称性！6 个独立的弹性常数为：
- **C₁₁ = 281.2 GPa**
- **C₁₂ = 158.1 GPa**
- **C₁₃ = 155.9 GPa**
- **C₃₃ = 484.0 GPa**
- **C₄₄ = 117.4 GPa**
- **C₆₆ = 216.4 GPa**

**物理意义**：
- C₃₃ (484.0 GPa) > C₁₁ (281.2 GPa)：沿 c 轴方向的刚度明显大于垂直方向，体现了四方晶系的各向异性
- C₆₆ (216.4 GPa) > C₄₄ (117.4 GPa)：垂直于 c 轴平面内的剪切刚度大于包含 c 轴的剪切刚度

### 5.4.4 弹性模量分析

计算得到的弹性模量：

| 弹性模量 | 数值 | 单位 |
|---------|------|------|
| 体模量（Bulk Modulus） | 220.7 | GPa |
| 剪切模量（Shear Modulus） | 128.7 | GPa |
| 杨氏模量（Young's Modulus） | 323.2 | GPa |
| 泊松比（Poisson's Ratio） | 0.256 | 无量纲 |

**与 Si 对比**：
- TiO₂ 的体模量（220.7 GPa）远大于 Si（90.7 GPa），说明 TiO₂ 更难压缩
- TiO₂ 的杨氏模量（323.2 GPa）也远大于 Si（157.7 GPa），说明 TiO₂ 更硬
- TiO₂ 的泊松比（0.256）略大于 Si（0.210）

## 5.5 与已知数据对比

### 5.5.1 详细对比表

将我们的计算结果与实验值和 Materials Project 数据对比：

| 弹性常数 | 本文结果 | Materials Project | 实验测量结果 | 本文 vs 实验 | MP vs 实验 |
|---------|---------|-------------------|------------|------------|-----------|
| C₁₁/GPa | 281.2   | 426               | 268.0 (1.4) | +4.9%     | +59.0%    |
| C₁₂/GPa | 158.1   | 2                 | 174.9 (1.4) | -9.6%     | -98.9%    |
| C₁₃/GPa | 155.9   | 149               | 147.4 (1.5) | +5.8%     | +1.1%     |
| C₃₃/GPa | 484.0   | 470               | 484.2 (1.8) | -0.04%    | -2.9%     |
| C₄₄/GPa | 117.4   | 113               | 123.8 (0.2) | -5.2%     | -8.7%     |
| C₆₆/GPa | 216.4   | 43                | 190.2 (0.5) | +13.8%    | -77.4%    |

**注**：实验值括号内为测量不确定度（GPa）。

### 5.5.2 结果分析

**本文计算结果的特点**：
- **C₃₃**：与实验值几乎完全一致（差异仅 0.04%）
- **C₁₁, C₁₃, C₄₄**：与实验值差异在 5-6% 范围内，吻合良好
- **C₁₂**：低估了约 10%
- **C₆₆**：高估了约 14%

总体而言，本文的计算结果**比 Materials Project 更接近实验值**，特别是：
- C₁₂：MP 严重低估（2 GPa vs 实验 174.9 GPa），本文结果合理（158.1 GPa）
- C₆₆：MP 严重低估（43 GPa vs 实验 190.2 GPa），本文结果合理（216.4 GPa）
- C₃₃：本文结果几乎完美，MP 也不错

**可能的原因**：
- 赝势和基组的选择：本文使用 APNS-v1 赝势 + efficiency 基组
- 结构优化的精度：本文进行了完全优化（cell-relax）
- 应力计算的精度：k 点网格和截断能的设置

### 5.5.3 与 Si 案例的对比

| 特征 | Si（立方晶系） | TiO₂（四方晶系） |
|-----|--------------|----------------|
| 独立弹性常数数量 | 3 | 6 |
| 对称性 | 高 | 中等 |
| 各向异性 | 低（C₁₁ ≈ C₃₃） | 高（C₃₃ >> C₁₁） |
| 体模量 | 90.7 GPa | 220.7 GPa |
| 与实验值差异 | 1-3% | 0-14% |
| 与 MP 差异 | 1-3% | 本文更接近实验 |

TiO₂ 的计算更具挑战性，但本文结果仍然可靠。

## 5.6 结果验证

### 5.6.1 对称性检查

检查弹性张量是否满足四方晶系（Laue 类型 I）的对称性：

✅ C₁₁ = C₂₂ = 281.2 GPa（相等）
✅ C₁₃ = C₂₃ = 155.9 GPa（相等）
✅ C₄₄ = C₅₅ = 117.4 GPa（相等）
✅ C₁₄ = C₁₅ = C₁₆ = C₂₄ = C₂₅ = C₂₆ = C₃₄ = C₃₅ = C₃₆ = C₄₅ = C₄₆ = C₅₆ ≈ 0

对称性完全满足！

### 5.6.2 力学稳定性检查

对于四方晶系，力学稳定性的判据为：

1. C₁₁ > 0 ✅（281.2 > 0）
2. C₃₃ > 0 ✅（484.0 > 0）
3. C₄₄ > 0 ✅（117.4 > 0）
4. C₆₆ > 0 ✅（216.4 > 0）
5. C₁₁ > |C₁₂| ✅（281.2 > 158.1）
6. C₁₁ + C₃₃ > 2C₁₃ ✅（281.2 + 484.0 = 765.2 > 2×155.9 = 311.8）

所有判据都满足，说明金红石型 TiO₂ 在计算的结构下是**力学稳定**的。

### 5.6.3 各向异性分析

计算各向异性比：

- **轴向刚度比**：C₃₃/C₁₁ = 484.0/281.2 = **1.72**
  - 说明沿 c 轴方向的刚度是垂直方向的 1.72 倍，各向异性明显

- **剪切刚度比**：C₆₆/C₄₄ = 216.4/117.4 = **1.84**
  - 说明垂直于 c 轴平面内的剪切刚度是包含 c 轴方向的 1.84 倍

这种各向异性是四方晶系的典型特征。

## 5.7 本章小结

本章通过金红石型 TiO₂ 的案例，演示了四方晶系材料弹性常数的计算：

1. **结构准备**：从 Materials Project 下载 TiO₂ 的 CIF 文件
2. **结构优化**：使用 abacustest 准备输入文件，运行 cell-relax 优化
3. **弹性计算**：生成 25 个变形结构，运行 SCF 计算
4. **结果分析**：得到 6 个独立的弹性常数和弹性模量
5. **结果验证**：与实验值和 Materials Project 对比，验证计算的准确性

**关键结果**：
- 6 个独立弹性常数：C₁₁=281.2, C₁₂=158.1, C₁₃=155.9, C₃₃=484.0, C₄₄=117.4, C₆₆=216.4 GPa
- 本文计算结果比 Materials Project 更接近实验值
- C₃₃ 与实验值几乎完全一致（差异仅 0.04%）
- 对称性和力学稳定性检查通过
- 明显的各向异性：C₃₃/C₁₁ = 1.72

**与 Si 案例的对比**：
- TiO₂ 的对称性较低（6 个独立常数 vs Si 的 3 个）
- TiO₂ 的各向异性更明显
- TiO₂ 的体模量和杨氏模量远大于 Si，更硬更难压缩
- 两个案例都验证了 abacustest 方法的可靠性

在下一章中，我们将讨论一些进阶话题和注意事项，帮助您更好地应用这些方法。
