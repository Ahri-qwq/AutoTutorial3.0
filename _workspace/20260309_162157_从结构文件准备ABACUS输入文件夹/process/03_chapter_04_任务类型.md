## 四、选择计算任务类型

通过 `--jtype` 参数切换计算类型，工具会自动调整 INPUT 中的相关参数。支持的任务类型：`scf`（默认）、`relax`、`cell-relax`、`md`、`band`。

以 MgO 晶胞优化为例：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --jtype cell-relax --folder-syntax MgO-cellrelax
```

生成的 `INPUT` 文件如下：

```
INPUT_PARAMETERS
calculation     cell-relax
symmetry     1
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.8
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
cal_force     1
cal_stress     1
kspacing     0.14 # unit in 1/bohr
relax_method     cg # or bfgs, bfgs_trad, cg_bfgs, sd, fire
relax_nmax     60
force_thr_ev     0.01  # unit in eV/A
stress_thr     0.5 # unit in kbar
fixed_axes     None # or volume, shape, a, b, c, ab, ac, bc; only valid for cell-relax calculation to fix some axes
#gamma_only
```

相比 SCF 的 INPUT，cell-relax 增加了：

- `cal_force 1` + `cal_stress 1`：计算力和应力
- `relax_method cg`：使用 CG 优化算法（适合变胞优化）
- `relax_nmax 60`：最大离子步数
- `force_thr_ev 0.01`：力收敛限（eV/Å）
- `stress_thr 0.5`：应力收敛限（kbar）
- `fixed_axes None`：可选固定某些晶格方向

若使用 `--jtype relax`，则 `calculation` 设为 `relax`，并自动去掉 `cal_stress` 和 `stress_thr`。
