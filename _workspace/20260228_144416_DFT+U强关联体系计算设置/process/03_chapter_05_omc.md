## 5. 进阶：Occupation Matrix Control (omc)

### 5.1 功能说明

在 DFT+U 计算中，SCF 迭代过程中 occupation matrix 会随着电子密度一起更新。对于复杂体系，不同的初始磁矩可能导致收敛到不同的局域极小（"多解问题"）。occupation matrix control（omc）提供了一种手动控制 occupation matrix 初始值的方法。

INPUT 中的 `omc` 参数有三种取值：

| omc 值 | 行为 |
|--------|------|
| 0 | 标准 DFT+U，不使用外部 occupation matrix，从原子初始值开始迭代 |
| 1 | 第一步读入 `initial_onsite.dm`，后续步骤照常更新 |
| 2 | 始终使用 `initial_onsite.dm` 中的 occupation matrix，全程固定，不更新 |

`initial_onsite.dm` 的文件格式与 `OUT.NiO/onsite.dm` 完全相同。

### 5.2 演示：使用 omc=2 固定 occupation matrix

#### 步骤 1：复制收敛后的 onsite.dm

将上一步计算生成的 `onsite.dm` 复制到工作目录，重命名为 `initial_onsite.dm`：

```bash
cp OUT.NiO/onsite.dm ./initial_onsite.dm
```

#### 步骤 2：修改 INPUT 文件

在 INPUT 中添加 `omc` 参数：

```
#Parameters (DFT+U)
dft_plus_u              1
orbital_corr            2 2 -1
hubbard_u               5.0 5.0 0.0
omc                     2        # 固定使用 initial_onsite.dm
```

#### 步骤 3：重新运行计算

```bash
OMP_NUM_THREADS=1 mpirun -n 8 abacus
```

#### 步骤 4：验证结果

使用 omc=2 的计算结果应与 omc=0 的标准 DFT+U 结果完全一致：

```
 !FINAL_ETOT_IS -9255.7279034240546025 eV   # 与标准计算相同
```

在 `running_scf.log` 中，每一步的 occupation matrix 应保持不变（等于 `initial_onsite.dm` 中的值）。

### 5.3 omc 的实际应用场景

- **多解问题**：当体系可能存在多个磁态（如铁磁、反铁磁、亚铁磁），可以通过指定不同的 `initial_onsite.dm` 引导收敛到目标磁态
- **重启计算**：对已知合理磁态的体系，用 omc=1 给出一个好的初始 occupation matrix，加快 SCF 收敛
- **固定轨道占据**：某些场景需要约束特定轨道的占据（如研究激发态），omc=2 可实现全程固定

### 5.4 补充：Yukawa Potential 自动确定 U 值

如果不想手动设置 U 值，ABACUS 还提供了通过 Yukawa potential 自动计算 U 值的方法：

```
yukawa_potential         1        # 开启 Yukawa potential 方法
# yukawa_lambda          X.X     # 可选：手动指定 screening length
```

开启后，程序在 SCF 过程中通过将电子相互作用近似为 Yukawa potential，自动计算各轨道的 U 值。更多细节参见 Qu et al. (2022)。
