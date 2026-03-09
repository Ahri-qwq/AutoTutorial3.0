## 五、案例三：态密度计算（Si，MDFT+DOS，T ≈ 8.16 eV）

### 5.1 场景说明

`186_PW_SDOS_10D10S` 是一个 1 个 Si 原子的体系，电子温度 0.6 Ry（约 8.16 eV），采用 MDFT（10 条 KS 轨道 + 10 条随机轨道）计算 SCF 并输出态密度（DOS）。算例名中的 `10D10S` 正是指 10 条确定性轨道和 10 条随机轨道。

与前两个算例相比，本例的核心增量是：在 SDFT/MDFT 的 SCF 基础上，通过一组专用参数控制态密度的能量范围、展宽和精度，最终在 `OUT.autotest/` 目录下生成 `DOS1_smearing.dat` 文件。

### 5.2 INPUT 文件

```
INPUT_PARAMETERS
#Parameters (1.General)
suffix          autotest
calculation     scf
esolver_type    sdft
method_sto      2
nbands          10
nbands_sto      10
nche_sto        120
emax_sto        0
emin_sto        0
seed_sto        20000
pseudo_dir      ../../PP_ORB
symmetry        1
kpar            1
bndpar          2
#Parameters (2.Iteration)
ecutwfc         20
scf_thr         1e-6
scf_nmax        20
#Parameters (3.Basis)
basis_type      pw
#Parameters (4.Smearing)
smearing_method fd
smearing_sigma  0.6
#Parameters (5.Mixing)
mixing_type     broyden
mixing_beta     0.4
out_dos         1
dos_emin_ev     -20
dos_emax_ev     100
dos_edelta_ev   0.1
dos_sigma       4
dos_nche        240
npart_sto       2
```

### 5.3 SDFT/MDFT 控制参数（本例特有项）

**`nbands = 10` / `nbands_sto = 10`**
MDFT 模式：10 条低能 KS 轨道 + 10 条随机轨道。相比案例一（4+64），本例 KS 轨道更多、随机轨道更少，这是因为 DOS 后处理对 KS 轨道的依赖更强——更多的 KS 轨道能更好地描述低能区域的态密度细节。

**`nche_sto = 120`**
同样 T=8.16 eV，本例用 120 阶（多于案例一的 100 阶）。DOS 计算对精度要求更高，建议适当增加展开阶数，确认 `Chebyshev Precision < 1e-8`。

**`emax_sto = 0` / `emin_sto = 0`**
随机轨道的能量范围上下限。设为 0 时程序自动确定，通常无需手动设置。

**`seed_sto = 20000`**
固定随机种子，确保计算结果可复现。收敛测试时应使用不同的 `seed_sto` 值（如 1、2、3……），以评估随机误差。

**`bndpar = 2`**
将所有 MPI 进程分成 2 组，随机轨道平均分配到每组，提高并行效率。实际使用时优先用 `kpar` 进行 K 点并行，再测试不同的 `bndpar` 取值。`bndpar` 并不是越大越好。

### 5.4 DOS 专用参数

| 参数 | 值 | 说明 |
|------|----|------|
| `out_dos` | 1 | 开关：设为 1 才会输出 DOS |
| `dos_emin_ev` | -20 eV | DOS 能量范围下限 |
| `dos_emax_ev` | 100 eV | DOS 能量范围上限 |
| `dos_edelta_ev` | 0.1 eV | DOS 能量间隔（分辨率） |
| `dos_sigma` | 4 eV | 高斯展宽宽度（WDM 体系用较大展宽） |
| `dos_nche` | 240 | DOS 后处理切比雪夫展开阶数（独立于 `nche_sto`） |
| `npart_sto` | 2 | 内存控制：使用正常内存的 1/2 |

**`dos_nche`** 与 `nche_sto` 是两个独立参数：前者专用于 DOS 的后处理步骤，通常需要比 `nche_sto` 取更大的值以获得高质量的态密度曲线。本例取 240，是 `nche_sto`（120）的两倍。

**`npart_sto`** 用于处理 `method_sto = 2` 时 DOS 后处理内存不足的问题。设为 `n` 时，程序将内存使用量控制在正常的 1/n，以较慢的速度换取内存的节省。默认值为 1（不拆分），建议在内存受限时设为 2 或更大。

### 5.5 其他参数

**`mixing_type = broyden` / `mixing_beta = 0.4`**
Broyden 混合方法，`mixing_beta` 控制每步更新的比例。SDFT 的 SCF 收敛有时比 KSDFT 稍慢，可适当降低 `mixing_beta`（如从 0.7 降至 0.4）来改善收敛。

**`ecutwfc = 20`**
本例截断能比案例一（50 Ry）更低，因为单原子 Si 的超胞更简单，且电子温度相同时 `nche_sto = 120`（阶数更多）已能弥补较低截断能带来的精度损失。

### 5.6 运行与输出

```bash
cd 186_PW_SDOS_10D10S
OMP_NUM_THREADS=1 mpirun -np 4 abacus
```

输出文件位于 `OUT.autotest/`：

- `running_scf.log`：SCF 收敛信息，含 `Chebyshev Precision` 检查项
- `DOS1_smearing.dat`：态密度数据文件，两列格式（能量(eV) / DOS强度）

### 5.7 注意事项

1. **DOS 能量范围选取**：`dos_emin_ev` 和 `dos_emax_ev` 应覆盖感兴趣的能量区间。对高温 WDM，费米面以上数十 eV 仍有可观态密度，上限取 100 eV 是合理的。

2. **`dos_sigma` 的选取**：展宽过小会导致 DOS 曲线有数值噪声（源于随机误差），展宽过大会掩盖能量细节。WDM 体系通常用 1~5 eV 的展宽。

3. **`kpar` 与 `bndpar` 的优先级**：`kpar`（K 点并行）的并行效率通常高于 `bndpar`。有多个 K 点时，优先使用 `kpar`，再调整 `bndpar`。
