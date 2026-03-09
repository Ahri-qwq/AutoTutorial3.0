## 第三章：案例一——Li 在 Si 中的扩散

### 3.1 体系简介

本案例计算 Li 原子在 Si 晶胞中从一个间隙位扩散到相邻间隙位的过程。这是一个典型的简单扩散体系：
- 只涉及单个 Li 原子的迁移，映像间原子对应关系清晰
- 初末态对称，能垒约 0.6 eV，NEB 链收敛较快
- 适合作为 NEB 入门案例，验证工作流设置

文件准备完成后，工作目录结构如下：

```
Li-diffu-Si/
├── INIT/                    # 初态计算结果（已优化）
│   └── OUT.ABACUS/running_scf.log
├── FINAL/                   # 末态计算结果（已优化）
│   └── OUT.ABACUS/running_scf.log
├── Li_ONCV_PBE-1.2.upf      # Li 赝势
├── Si_ONCV_PBE-1.2.upf      # Si 赝势
├── Li_gga_8au_100Ry_4s1p.orb     # Li 轨道文件
└── Si_gga_8au_100Ry_2s2p1d.orb   # Si 轨道文件
```

### 3.2 生成初猜 NEB 链

用 `neb_make.py` 读取初末态，通过 IDPP 插值生成初猜：

```bash
python neb_make.py -i INIT/OUT.ABACUS/running_scf.log \
                      FINAL/OUT.ABACUS/running_scf.log \
                   -n 3
```

**参数说明：**

| 参数 | 含义 |
|------|------|
| `-i` | 初态和末态的输出文件路径（支持 abacus-out 格式） |
| `-n` | 中间映像数（不含初末态）；总映像数 = n + 2 |
| `-m` | 插值方法，默认 `IDPP`（Pymatgen 版本），优于 ASE 原生 IDPP |
| `--fix` | 固定某一高度以下的原子（`[高度]:[方向]`，用于表面计算） |
| `--mag` | 设置初始磁矩（`元素:磁矩,...`） |

执行后输出 `init_neb_chain.traj`，共 5 个映像（初态 + 3 个中间态 + 末态）。

```
Reading files: ./INIT/OUT.ABACUS/running_scf.log and ./FINAL/OUT.ABACUS/running_scf.log
Generating path, number of images: 3, sort_tol: 1.0
Optimizing path using IDPP method
Writing path: init_neb_chain.traj, Number of images: 5
```

> **选择映像数：** 可先用 `neb_dist.py` 查看初末态原子间距，以此为参考。简单扩散取 3 个映像通常足够；复杂表面反应建议 5–8 个。

### 3.3 串行 DyNEB 计算

将 `neb_run.py` 复制到工作目录并编辑，以下是本案例使用的完整脚本（`dyneb_run.py`）：

```python
# dyneb_run.py — Li-diffu-Si 串行 DyNEB

from ase.optimize import FIRE
from ase.io import read
from abacus_neb import AbacusNEB

# ===== 并行与优化设置 =====
mpi = 16           # ABACUS 每个映像使用的 MPI 进程数
omp = 1            # OpenMP 线程数
fmax = 0.05        # 收敛力阈值（eV/Å）
neb_optimizer = FIRE
neb_directory = "NEBrun"
algorism = "improvedtangent"  # IT-NEB，推荐
climb = True                  # 启用 CI-NEB
dyneb = True                  # 串行 DyNEB 模式
parallel = False              # 串行，不使用 GPAW MPI
k = 0.10                      # 弹簧力常数（eV/Å²）
init_chain = "init_neb_chain.traj"
abacus = "abacus"

# ===== 赝势与轨道 =====
lib_dir = ""          # 赝势/轨道与脚本同目录时留空
pseudo_dir = lib_dir
basis_dir  = lib_dir
pp = {
    'Li': 'Li_ONCV_PBE-1.2.upf',
    'Si': 'Si_ONCV_PBE-1.2.upf',
}
basis = {
    'Li': 'Li_gga_8au_100Ry_4s1p.orb',
    'Si': 'Si_gga_8au_100Ry_2s2p1d.orb',
}

# ===== ABACUS INPUT 参数 =====
kpts = [2, 2, 2]
parameters = {
    'calculation':    'scf',
    'nspin':          1,
    'xc':             'pbe',
    'ecutwfc':        100,
    'dft_functional': 'pbe',
    'ks_solver':      'genelpa',
    'symmetry':       0,
    'vdw_method':     'none',
    'smearing_method':'gaussian',
    'smearing_sigma': 0.001,
    'basis_type':     'lcao',
    'mixing_type':    'broyden',
    'scf_thr':        1e-6,
    'scf_nmax':       100,
    'kpts':           kpts,
    'pp':             pp,
    'basis':          basis,
    'pseudo_dir':     pseudo_dir,
    'basis_dir':      basis_dir,
    'init_wfc':       'atomic',
    'init_chg':       'atomic',
    'cal_force':      1,
    'cal_stress':     1,
    'out_stru':       1,
    'out_chg':        -1,
}

if __name__ == "__main__":
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters,
                    parallel=parallel, directory=neb_directory,
                    mpi=mpi, omp=omp, abacus=abacus,
                    algorism=algorism, k=k, dyneb=dyneb)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)
```

**关键参数说明：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `dyneb` | `True` | 启用 Dynamic NEB，串行下自动跳过已收敛映像 |
| `parallel` | `False` | 串行模式，用 `python dyneb_run.py` 运行 |
| `symmetry` | `0` | NEB 计算必须关闭对称性 |
| `cal_force` | `1` | NEB 优化需要原子受力 |
| `out_chg` | `-1` | 每步覆盖保存电荷密度（节省磁盘空间） |

串行运行：

```bash
python dyneb_run.py
```

计算运行时，ASE 会打印收敛过程：

```
Notice: Parallel calculation is not being used
Set NEB method to DyNEB automatically
----- Running Dynamic NEB -----
----- improvedtangent method is being used -----
      Step     Time          Energy          fmax
FIRE:    0 18:28:24    -7051.845345         0.733359
FIRE:    1 18:34:09    -7051.879062         0.594320
...
FIRE:   23 20:28:27    -7051.985587         0.036730
```

本案例在 `c32_m64_cpu`（32核）机器上约需 1–2 小时收敛。

### 3.4 并行 NEB 对比

并行 NEB 将每个映像分配给一个独立进程，多映像同时计算，适合核数充足时使用。修改 `neb_run.py` 中的关键参数：

```python
mpi      = 5       # 每个进程内 ABACUS 使用的核数（无强制要求）
omp      = 1
dyneb    = False   # 并行模式不使用 DyNEB
parallel = True    # 启用并行
```

并行启动需要通过 GPAW 的 MPI 环境：

```bash
mpirun -np 3 gpaw python neb_run.py
```

> **注意：** `mpirun -np` 的进程数应与中间映像数一致（本例为 3）。ABACUS 和 GPAW 需依赖同一套 MPI 环境，镜像中通常已配置好。

并行 NEB 与串行 DyNEB 的计算结果基本一致，但挂钟时间（walltime）显著缩短。

### 3.5 后处理与结果分析

计算完成后，用 `neb_post.py` 处理轨迹文件：

```bash
# NEB 后处理
python neb_post.py neb.traj
```

输出：
- `nebplots_chain.pdf`：NEB 链的势能-反应坐标曲线（含各映像切线力投影）
- `neb_latest.traj`：收敛时刻的 NEB 链轨迹

终端还会打印每个映像的绝对能量与能垒：

```
num: 0; Energy: -7052.6037618 (eV)
num: 1; Energy: -7052.2767038 (eV)
num: 2; Energy: -7051.9855868 (eV)   ← 过渡态（能量最高）
num: 3; Energy: -7052.2797724 (eV)
num: 4; Energy: -7052.6035978 (eV)
Reaction Barrier and Energy Difference: (0.618 eV, 0.0002 eV)
```

| 结果量 | 数值 |
|--------|------|
| 扩散能垒 | 0.618 eV |
| 反应能差 | 0.0002 eV（近似对称路径） |

可用 ASE 可视化 NEB 链：

```bash
ase -T gui neb_latest.traj
```

或在 Jupyter/Bohrium 中：

```python
from ase.io import read
from ase.visualize import view

atoms = read("neb_latest.traj", ":")
view(atoms, viewer='ngl')
```
