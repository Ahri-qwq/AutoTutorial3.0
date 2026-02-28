# 四、案例二：金红石型 TiO₂ 的弹性常数

金红石型 TiO₂（Materials Project 编号 mp-2657）属于四方晶系，Laue 类型 I，有 6 个独立弹性张量元：C₁₁（= C₂₂）、C₁₂、C₁₃（= C₂₃）、C₃₃、C₄₄（= C₅₅）、C₆₆。计算流程与 Si 完全相同，这里重点展示四方晶系的结果特征及与实验数据的对比。

## 4.1 获取结构

从 Materials Project 下载金红石型 TiO₂（mp-2657）的 CIF 文件，命名为 `TiO2.cif`。

## 4.2 结构优化

用 abacustest 准备 cell-relax 输入文件并提交优化计算：

```bash
abacustest model inputs -f TiO2.cif --ftype cif --lcao --jtype cell-relax --folder-syntax TiO2-rutile
```

等待计算完成后，整理优化后的结构：

```bash
mkdir TiO2-rutile-elastic
cp TiO2-rutile/INPUT TiO2-rutile/Ti* TiO2-rutile/O* TiO2-rutile-elastic
cp TiO2-rutile/INPUT/OUT.ABACUS/STRU_ION_D TiO2-rutile-elastic/
```

> **说明：** `TiO2-rutile/Ti*` 和 `TiO2-rutile/O*` 分别复制 Ti 和 O 的赝势与轨道文件（文件名以元素符号开头）。根据实际文件名调整通配符。

## 4.3 准备并提交弹性计算

同样地，修改 `TiO2-rutile-elastic/INPUT` 中的计算类型为 `scf`，然后运行：

```bash
abacustest model elastic prepare -j TiO2-rutile-elastic
```

同样在 `TiO2-rutile-elastic/` 下生成 25 个计算任务（24 个形变 + 1 个原始结构）。提交所有任务，等待完成。

> **警告：** 与 Si 案例相同，不要在计算进行中或计算完成后重复运行 `elastic prepare`，否则已有的计算文件夹将被删除。

## 4.4 后处理

```bash
abacustest model elastic post -j TiO2-rutile-elastic
```

屏幕输出：

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

## 4.5 结果分析

### 弹性张量结构的验证

首先检查弹性张量矩阵结构是否符合四方晶系（Laue I）的预期：
- C₁₁ = C₂₂：281.2 ≈ 281.2 GPa ✓
- C₁₃ = C₂₃：155.9 ≈ 155.9 GPa ✓
- C₄₄ = C₅₅：117.4 ≈ 117.4 GPa ✓
- C₆₆ ≠ C₄₄：216.4 ≠ 117.4 GPa ✓（四方晶系区别于六方晶系的关键特征）
- 非对角非零元均为数值噪声量级（~10⁻⁸ GPa）✓

计算得到的弹性张量结构完全符合金红石型 TiO₂ 所属的四方晶系（Laue 类型 I）对 6 个独立弹性张量元的要求。

### 与 Materials Project 及实验数据的对比

下表汇总了本文计算、Materials Project 数据库和实验测量 [3] 的各独立弹性张量元：

| | 本文结果 (GPa) | Materials Project (GPa) | 实验测量 (GPa) |
|---|---:|---:|---:|
| C₁₁ | 281.2 | 426 | 268.0 ± 1.4 |
| C₁₂ | 158.1 | 2 | 174.9 ± 1.4 |
| C₁₃ | 155.9 | 149 | 147.4 ± 1.5 |
| C₃₃ | 484.0 | 470 | 484.2 ± 1.8 |
| C₄₄ | 117.4 | 113 | 123.8 ± 0.2 |
| C₆₆ | 216.4 | 43 | 190.2 ± 0.5 |

**对比分析：**

本文结果与实验数据总体吻合较好，多数分量的偏差在 10% 以内。C₃₃ 和 C₄₄ 与实验的吻合尤为出色（偏差 < 5%）。

Materials Project 数据库中部分分量与实验值差异显著，尤其是 C₁₂（数据库值 2 GPa，实验值 174.9 GPa）和 C₆₆（数据库值 43 GPa，实验值 190.2 GPa），这可能与该数据库条目的计算参数设置有关。本文的计算结果在这些分量上明显更接近实验值。

**各向同性力学性质（单位：GPa，无量纲量泊松比除外）：**

| 性质 | 数值 |
|------|------|
| 体模量（Bulk Modulus） | 220.71 |
| 剪切模量（Shear Modulus） | 128.68 |
| 杨氏模量（Young's Modulus） | 323.23 |
| 泊松比（Poisson's Ratio） | 0.256 |

金红石型 TiO₂ 的体模量（220 GPa）和杨氏模量（323 GPa）均显著高于 Si，反映了 TiO₂ 更强的化学键合和更高的力学刚度。
