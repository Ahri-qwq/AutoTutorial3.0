## 五、进阶：Occupation Matrix Control

### 5.1 背景与动机

DFT+U 计算存在一个已知问题：**多解性**。对于同一个体系，不同的初始磁矩或初始电荷密度可能收敛到不同的局域极小值，对应不同的 occupation matrix，从而得到不同的能量和物理性质。

Occupation Matrix Control（`omc`）功能允许用户在 SCF 迭代过程中固定 occupation matrix，可以帮助：

1. 将体系引导到特定的电子态（避免陷入非目标极小值）
2. 验证收敛后的结果是否稳定

### 5.2 三种模式

ABACUS 提供三种 `omc` 模式：

| `omc` 值 | 行为 |
|----------|------|
| 0 | 不使用 occupation matrix control，标准 DFT+U 计算（默认） |
| 1 | SCF 第一步读入 `initial_onsite.dm`，后续正常更新 occupation matrix |
| 2 | 读入 `initial_onsite.dm`，整个 SCF 过程中始终使用此 occupation matrix，不更新 |

### 5.3 使用流程

以验证 NiO 计算为例，使用 omc = 2：

**第一步：** 从上一步计算的 `OUT.NiO/onsite.dm` 复制到工作目录，改名为 `initial_onsite.dm`：

```bash
cp OUT.NiO/onsite.dm ./initial_onsite.dm
```

**第二步：** 在 INPUT 文件中添加 `omc 2`：

```
#Parameter DFT+U
dft_plus_u    1
orbital_corr  2 2 -1
hubbard_u     5.0 5.0 0.0
omc           2
```

**第三步：** 重新运行 ABACUS，在 `running_scf.log` 中可以确认整个 SCF 过程中 occupation matrix 保持不变。

**验证结果：** 使用 omc = 2 固定收敛后的 occupation matrix 重新计算，最终得到的总能量和磁矩与标准 SCF 完全一致，说明之前的 SCF 已收敛到稳定的基态。

### 5.4 initial_onsite.dm 文件格式

`initial_onsite.dm` 的格式与 `onsite.dm` 完全相同，也与日志中 `L(S)DA+U` 块的输出格式一致。用户可以手动编辑此文件，指定自定义的 occupation matrix，以引导体系进入特定的电子态。

### 5.5 进阶参数：Yukawa Potential 自动确定 U 值

除上述参数外，ABACUS 还支持通过 Yukawa 势在 SCF 中自动计算 U 值，无需手动设定：

| 参数 | 说明 |
|------|------|
| `yukawa_potential 1` | 开启 Yukawa 势方法，程序内自行计算 U 值 |
| `yukawa_lambda` | 手动指定 Yukawa 势的 screening length；不设则由电子密度 on-the-fly 计算 |

这种方法适合不确定 U 值的场合，理论细节见参考文献 [1]。
