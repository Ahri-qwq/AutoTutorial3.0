## 四、案例二：分子动力学模拟（Al，纯SDFT，T ≈ 100 eV）

### 4.1 场景说明

`pw_md_Al` 是一个 16 个铝原子的体系，电子温度设为 7.35 Ry（约 100 eV）。本例使用**纯 SDFT**（`nbands = 0`，没有任何 KS 轨道），执行分子动力学（MD）模拟。

极端高温（100 eV）是 SDFT 最能发挥优势的场景：此时 KSDFT 需要数百条以上的波函数，而 SDFT 只需 64 条随机轨道，且切比雪夫展开仅需 20 阶（远少于低温时的 100 阶），计算效率大幅提升。

### 4.2 INPUT 文件

```
INPUT_PARAMETERS
#Parameters (General)
calculation    md
esolver_type   sdft
pseudo_dir     ../../PP_ORB
nbands         0
nbands_sto     64
nche_sto       20
method_sto     2
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
scf_thr        1e-6
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  7.34986072
#Parameters (MD)
md_tfirst      1160400
md_dt          0.2
md_nstep       10
```

### 4.3 关键参数详解

**`calculation = md`**
设为 `md` 启动分子动力学模拟。在每个 MD 步中，程序先执行 SDFT 的 SCF 自洽迭代计算出受力，再根据受力更新原子位置。

**`nbands = 0`**
本例不使用任何 KS 轨道，全部依赖随机轨道，是纯 SDFT 模式。与案例一（`nbands = 4`，MDFT）相比，适用于极高温场景——此时费米面附近已无明显的低能 KS 轨道可供利用。

**`nche_sto = 20`**
温度 100 eV 时，费米-狄拉克函数变得非常平滑，切比雪夫展开只需 20 阶即可达到足够精度。相比案例一（T=8.16 eV，100 阶），高温大幅降低了展开阶数，这也是 SDFT 在高温下效率优异的根本原因。

**`method_sto = 2`**
高温 MD 计算通常需要更快的速度，选用方法 2（速度更快，内存更多）。

**`smearing_sigma = 7.34986072`**
电子温度，单位 Ry，约等于 100 eV（7.35 Ry × 13.6 eV/Ry ≈ 99.96 eV）。

**MD 参数**

| 参数 | 值 | 说明 |
|------|----|------|
| `md_tfirst` | 1160400 K | 离子初始温度（1160400 K ÷ 11604.5 K/eV ≈ 100 eV） |
| `md_dt` | 0.2 fs | MD 时间步长 |
| `md_nstep` | 10 | MD 模拟总步数 |

注意 `md_tfirst` 是**离子温度**（单位 K），而 `smearing_sigma` 是**电子温度**（单位 Ry）。高温 WDM 模拟中，两者通常取相同的等效温度。

### 4.4 纯 SDFT 还是 MDFT？

| 场景 | 推荐模式 | 原因 |
|------|----------|------|
| 温度 < 10 eV | MDFT（`nbands > 0`） | 低能 KS 轨道可显著减少随机误差 |
| 温度 ≥ 10 eV | 纯 SDFT（`nbands = 0`）或 MDFT | SDFT 效率已足够，MDFT 收益递减 |
| 极端高温（> 50 eV）| 纯 SDFT | KS 轨道对角化成本高，不值得 |

### 4.5 运行方法

```bash
cd pw_md_Al
OMP_NUM_THREADS=1 mpirun -np 8 abacus
```

MD 轨迹和各步能量、受力、压强等信息输出在 `OUT.ABACUS/` 目录的 `running_md.log` 和 `MD_dump` 文件中。

### 4.6 注意事项

1. **`scf_thr = 1e-6`**：MD 中每一步的 SCF 不需要收敛到静态计算的精度，`1e-6` 已足够。
2. **高温赝势**：100 eV 时铝的内壳层电子会被热激发，需要使用包含更多内壳层电子的赝势，或确认所用赝势在 100 eV 下的适用范围。
3. **SDFT 随机误差在 MD 中的影响**：每步 SCF 的受力含有随机误差，会引入额外的随机"热噪声"。通过足够多的 `nbands_sto` 可将误差控制在可接受范围。
