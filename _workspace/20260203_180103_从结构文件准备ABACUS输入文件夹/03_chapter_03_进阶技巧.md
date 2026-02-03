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
