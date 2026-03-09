## 三、案例一：SCF 计算（Si，MDFT，T ≈ 8.16 eV）

### 3.1 场景说明

`pw_Si2` 是一个 2 个原子的金刚石结构硅体系，电子温度设为 0.6 Ry（约 8.16 eV）。本例采用 MDFT 模式：4 条 KS 轨道（覆盖费米面以下的低能部分）加 64 条随机轨道，执行电子自洽迭代（SCF）计算。

这是 SDFT 功能最基础的入口案例，所涉及的核心参数在此一并详解。

### 3.2 INPUT 文件

```
INPUT_PARAMETERS
#Parameters (General)
calculation    scf
esolver_type   sdft
pseudo_dir     ../../PP_ORB
nbands         4
nbands_sto     64
nche_sto       100
method_sto     1
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  0.6
```

### 3.3 关键参数详解

**`esolver_type`**
选择系统总能量的求解方式。默认值为 `ksdft`，设为 `sdft` 后程序将根据 `nbands` 和 `nbands_sto` 的取值自动判断使用纯 SDFT 还是 MDFT。

**`nbands`**
确定性 KS 轨道数（deterministic orbitals），通过严格对角化矩阵计算得到。一般取费米面以下所有低能轨道的数目，效率最高。本例中取 4，对应 Si 的 4 个价电子轨道。

- `nbands = 0` 且 `nbands_sto > 0`：纯 SDFT
- `nbands > 0` 且 `nbands_sto > 0`：MDFT（本例）

**`nbands_sto`**
随机波函数轨道数（stochastic orbitals）。取值越大，随机误差越小，但计算量相应增加。本例取 64。

如何判断 64 是否足够？见第六章的收敛测试方法。

**`nche_sto`**
切比雪夫展开阶数。对哈密顿量进行切比雪夫多项式展开以评估费米-狄拉克函数，阶数越高精度越高、效率越低。

经验规律：
- 与温度**成反比**：温度越高，所需阶数越少
- 与 `ecutwfc`（平面波截断能）**成正比**：截断能越大，所需阶数越多

判断标准：查看输出文件 `running_scf.log`，找到 `Chebyshev Precision` 一行，确保其值小于 `1e-8`。本例 T=8.16 eV，取 100 阶。

**`method_sto`**
SDFT 内部计算方法的选择：
- `1`：内存消耗较小，速度稍慢
- `2`：速度更快，但需要更多内存（默认值）

本例取 `1`，适合内存有限的情形。

**`smearing_method`**
展宽方式。SDFT/MDFT 当前**只支持** `fd`（Fermi-Dirac）展宽，不能使用高斯或 Methfessel-Paxton 等方式。

**`smearing_sigma`**
电子温度，单位为 **Ry**（里德伯，1 Ry ≈ 13.6 eV）。本例 0.6 Ry ≈ 8.16 eV。注意这是电子温度，在 MD 计算中与离子温度（`md_tfirst`）是两个独立参数。

**常规参数**

| 参数 | 值 | 说明 |
|------|----|------|
| `ecutwfc` | 50 Ry | 平面波截断能 |
| `scf_nmax` | 20 | 最大 SCF 迭代步数 |
| `symmetry` | 1 | 开启晶体对称性 |

### 3.4 STRU 和 KPT

STRU 文件描述金刚石结构 Si 的晶胞和原子位置，格式与常规 KSDFT 计算完全相同。KPT 文件设置布里渊区 K 点采样，ABACUS 的 SDFT 支持多 K 点，在某些性质（如 DOS）计算中要注意 K 点收敛性。

### 3.5 运行方法

```bash
cd pw_Si2
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

计算完成后，结果位于 `OUT.ABACUS/` 目录，总能量和收敛信息可在 `running_scf.log` 中查看。

### 3.6 注意事项

1. **赝势高温可移植性**：本例 T=8.16 eV 时，Si 内壳层未被热激发，使用标准赝势是合理的。当温度高于 ~100 eV 时需要特别评估赝势的适用性。

2. **K 点收敛**：SDFT 支持多 K 点采样，高精度计算前应测试 K 点数目对结果的影响。

3. **`nbands_sto = 0` 的特殊行为**：ABACUS 3.2.2 版本以后，若将 `nbands_sto` 设为 0，程序会自动切换为 KSDFT 模式计算。
