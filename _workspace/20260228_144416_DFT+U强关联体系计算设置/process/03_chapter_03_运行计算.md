## 3. 运行计算

### 3.1 运行 SCF 计算

在 `ABACUS_DFT+U/` 目录下执行：

```bash
# 进入工作目录
cd ABACUS_DFT+U

# 使用 8 个 MPI 进程，单线程运行
OMP_NUM_THREADS=1 mpirun -n 8 abacus
```

根据机器配置调整 MPI 进程数。计算在普通工作站上通常几分钟内完成。

### 3.2 检查计算收敛

计算完成后，检查 `OUT.NiO/running_scf.log`。正常收敛时，SCF 迭代的输出类似：

```
 ITER   TMAG      AMAG      ETOT(eV)          EDIFF(eV)     DRHO       TIME(s)
 GE1    0.00e+00  3.72e+00  -9.253482e+03   0.000000e+00  1.204e-01   ...
 GE2    0.00e+00  3.64e+00  -9.255098e+03  -1.616e+00     2.897e-02   ...
 ...
 GE30   0.00e+00  3.35e+00  -9.255728e+03  -1.2e-08       8.3e-09     ...
```

关键检查项：
- `ETOT` 趋于稳定：相邻步差值小于 `scf_thr`（本例为 1e-7 eV）
- `TMAG` 始终为 0：确认体系保持反铁磁态（总磁矩为零）
- `AMAG` 收敛到约 3.35：绝对磁矩反映两个 Ni 的磁矩大小

若 SCF 不收敛，可尝试减小 `mixing_beta`（如 0.2）或增大 `scf_nmax`。
