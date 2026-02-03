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
