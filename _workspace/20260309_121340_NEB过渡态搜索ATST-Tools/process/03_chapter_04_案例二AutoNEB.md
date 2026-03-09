## 第四章：案例二——环己烷在 Pt@石墨烯上的 AutoNEB

### 4.1 体系简介与方法选择

本案例计算环己烷（C₆H₁₂）在 Pt 原子负载石墨烯（Cy-Pt@graphene）表面上的 C-H 键断裂反应：

$$\text{C}_6\text{H}_{12}\text{@Pt-graphene} \rightarrow \text{C}_6\text{H}_{11}\text{@Pt-graphene} + \text{H}$$

这一体系比案例一复杂得多：
- 原子数多，势能面较"软"，NEB 映像间的原子匹配难度更高
- 涉及 Pt 的 d 轨道和石墨烯的 π 电子，电子结构自洽收敛较慢
- 直接用 IT-NEB 计算（8 个映像）验证表明会收敛到**错误的过渡态**

AutoNEB 通过动态调整映像分布，在过渡态附近集中计算资源，使该体系得到正确的过渡态（能垒 1.328 eV）。

### 4.2 初末态准备与插值

工作目录准备：

```
Cy-Pt@graphene/
├── IS/                      # 初态（已优化）
│   └── OUT.ABACUS/running_relax.log
├── FS/                      # 末态（已优化）
│   └── OUT.ABACUS/running_relax.log
├── C_ONCV_PBE-1.0.upf
├── H_ONCV_PBE-1.0.upf
├── Pt_ONCV_PBE-1.0.upf
├── C_gga_7au_100Ry_2s2p1d.orb
├── H_gga_6au_100Ry_2s1p.orb
└── Pt_gga_7au_100Ry_4s2p2d1f.orb
```

生成初猜（4 个中间映像，用于 AutoNEB 的初始映像数 = `n_simul`）：

```bash
python neb_make.py -i IS/OUT.ABACUS/running_relax.log \
                      FS/OUT.ABACUS/running_relax.log \
                   -n 4
```

```
Reading files: IS/OUT.ABACUS/running_relax.log and FS/OUT.ABACUS/running_relax.log
Generating path, number of images: 4, sort_tol: 1.0
Optimizing path using IDPP method
Writing path: init_neb_chain.traj, Number of images: 6
```

> Pymatgen 的 IDPP 在匹配失败时会发出警告（`Auto sorting is turned off`），通常不影响计算，AutoNEB 会在优化过程中自动修正映像位置。

### 4.3 AutoNEB 脚本配置

将 `autoneb_run.py` 复制到工作目录并编辑：

```python
# autoneb_run.py — Cy-Pt@graphene AutoNEB

from ase.optimize import FIRE
from ase.io import read
from ase.parallel import world
from abacus_autoneb import AbacusAutoNEB

# ===== AutoNEB 特有设置 =====
mpi          = 16
omp          = 4
neb_optimizer = FIRE
neb_directory = "AutoNEBrun"
algorism     = "improvedtangent"
init_chain   = "init_neb_chain.traj"
climb        = True
fmax         = [0.20, 0.05]    # 两阶段收敛阈值：粗→精
n_simul      = world.size      # 同时计算的映像数 = MPI 进程数
n_images     = 10              # 最终映像总数上限
smooth_curve = False
k            = 0.10

abacus     = "abacus"
lib_dir    = ""
pseudo_dir = f"{lib_dir}/"
basis_dir  = f"{lib_dir}/"

# ===== 赝势与轨道 =====
pp = {
    'C':  'C_ONCV_PBE-1.0.upf',
    'H':  'H_ONCV_PBE-1.0.upf',
    'Pt': 'Pt_ONCV_PBE-1.0.upf',
}
basis = {
    'C':  'C_gga_7au_100Ry_2s2p1d.orb',
    'H':  'H_gga_6au_100Ry_2s1p.orb',
    'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb',
}

# ===== ABACUS INPUT 参数 =====
kpts = [2, 1, 2]
parameters = {
    'calculation':    'scf',
    'nspin':          2,           # 有 Pt 原子，开启自旋
    'xc':             'pbe',
    'ecutwfc':        100,
    'dft_functional': 'pbe',
    'ks_solver':      'genelpa',
    'symmetry':       0,
    'vdw_method':     'd3_bj',     # 范德华校正（表面吸附体系推荐）
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
    'out_chg':        0,
    'out_mul':        0,
    'out_wfc_lcao':   0,
    'out_bandgap':    0,
    'efield_flag':    1,           # 电场修正（表面偶极校正）
    'dip_cor_flag':   1,
    'efield_dir':     1,           # 修正方向：y（真空层方向）
}

if __name__ == "__main__":
    init_chain = read(init_chain, index=':')
    neb = AbacusAutoNEB(init_chain, parameters,
                        algorism=algorism,
                        directory=neb_directory,
                        k=k, n_simul=n_simul, n_max=n_images,
                        abacus=abacus, mpi=mpi, omp=omp)
    neb.run(optimizer=neb_optimizer, climb=climb,
            fmax=fmax, smooth_curve=smooth_curve)
```

**与案例一（neb_run.py）的关键差异：**

| 参数 | 案例一 NEB | 案例二 AutoNEB | 说明 |
|------|-----------|---------------|------|
| `fmax` | `0.05` | `[0.20, 0.05]` | AutoNEB 分两阶段收敛 |
| `n_images` | — | `10` | AutoNEB 最终映像上限 |
| `n_simul` | — | `world.size` | 同时优化的映像数 |
| `nspin` | `1` | `2` | 含 Pt，需要自旋极化 |
| `vdw_method` | `none` | `d3_bj` | 表面吸附需要色散修正 |
| `efield_flag` | `0` | `1` | 表面偶极修正 |

并行启动（4 进程对应 4 个同时计算的映像）：

```bash
mpirun -np 4 gpaw python autoneb_run.py
```

### 4.4 AutoNEB 迭代过程

AutoNEB 的运行日志（`running_autoneb.out`）记录了动态映像添加过程：

```
NSIMUL is 4
===== AutoNEB Job Starting =====
----- Running AutoNEB -----
----- improvedtangent method is being used -----
The NEB initially has 6 images (including the end-points)

Start of evaluation of the initial images
Now starting iteration 1 on [0, 1, 2, 3, 4, 5]
Finished initialisation phase.

****Now adding another image until n_max is reached (6/10)****
Adding image between 2 and 3.
Now starting iteration 2 on [1, 2, 3, 4, 5, 6]

****Now adding another image until n_max is reached (7/10)****
Adding image between 1 and 2.
Now starting iteration 3 on [1, 2, 3, 4, 5, 6]

****Now adding another image until n_max is reached (8/10)****
Adding image between 5 and 6.
Now starting iteration 4 on [2, 3, 4, 5, 6, 7]

****Now adding another image until n_max is reached (9/10)****
Adding image between 7 and 8.
Now starting iteration 5 on [4, 5, 6, 7, 8, 9]

n_max images has been reached
****Now doing the CI-NEB calculation****
Now starting iteration 6 on [2, 3, 4, 5, 6, 7]
----- AutoNEB calculation finished -----
```

**6 次迭代的含义：**

| 迭代 | 内容 | 映像数 |
|------|------|--------|
| iter 1 | 评估初始 6 个映像 | 6 |
| iter 2–5 | 动态添加映像（每次添加 1 个，聚焦于差异最大处） | 7 → 10 |
| iter 6 | 最终 CI-NEB，精细收敛过渡态 | 10 |

最终 CI-NEB（iter 6）的收敛过程：

```
      Step     Time          Energy          fmax
FIRE:    0 01:57:54   -11866.084572         0.652487
...
FIRE:   88 03:06:42   -11865.497000         0.040039   ← 收敛（fmax < 0.05）
```

### 4.5 后处理与结果

AutoNEB 后处理使用 `--autoneb` 模式：

```bash
python neb_post.py --autoneb run_autoneb???.traj
```

Shell 会自动按数字顺序展开 `run_autoneb000.traj` 到 `run_autoneb009.traj`。

各映像能量与结果：

```
num: 0; Energy: -11866.824886 (eV)
num: 1; Energy: -11866.827220 (eV)
num: 2; Energy: -11866.801706 (eV)
num: 3; Energy: -11866.782620 (eV)
num: 4; Energy: -11866.587954 (eV)
num: 5; Energy: -11865.497000 (eV)   ← 过渡态（CI 映像）
num: 6; Energy: -11866.330196 (eV)
num: 7; Energy: -11866.396832 (eV)
num: 8; Energy: -11866.429487 (eV)
num: 9; Energy: -11866.435419 (eV)
Reaction Barrier and Energy Difference: (1.328 eV, 0.389 eV)
```

| 结果量 | 数值 |
|--------|------|
| C-H 断裂能垒 | 1.328 eV |
| 反应能差（H 解离吸热量） | 0.389 eV |

**判断过渡态是否正确：** 在 `nebplots_all.pdf` 中，过渡态处的切线应基本水平（切线力接近零）。最终通过振动分析确认唯一虚频（见第五章）。

可视化 AutoNEB 过程（实时监控）：

```bash
# 查看最新 NEB 路径
ase -T gui -n -1 run_autoneb???.traj

# 查看最终收敛链
ase -T gui neb_latest.traj
```
