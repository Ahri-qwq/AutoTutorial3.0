# 使用 ABACUS 和 ATST-Tools 进行 NEB/AutoNEB 过渡态搜索

## 引言

过渡态是化学反应机理分析的核心对象。无论是固体中的原子扩散、表面催化反应，还是均相有机反应，微观机理研究都离不开对过渡态的定位与分析。

本教程介绍如何用 **ATST-Tools** 工具包调用 **ABACUS** 进行过渡态搜索计算。ATST-Tools 是基于 ASE 构建的工作流脚本集，封装了 NEB、AutoNEB、Dimer、Sella 等多种过渡态方法，并通过 ASE-ABACUS 接口将 ABACUS 作为电子结构计算引擎。

**学完本教程，你将能够：**
- 理解 NEB 与 AutoNEB 方法的原理与区别
- 配置 ATST-Tools 环境，组织 NEB 计算目录
- 完成从初末态准备、初猜插值、NEB 计算到后处理的全流程操作
- 用振动分析验证过渡态的正确性

**前置条件：**
- 已安装 ABACUS 及 ASE-ABACUS 接口
- 已安装 ATST-Tools 及其依赖（pymatgen、pymatgen-analysis-diffusion、GPAW）
- 熟悉 ABACUS 结构优化计算（初末态需要提前优化完毕）

本教程涵盖两个计算案例：

| 案例 | 体系 | 方法 | 能垒 |
|------|------|------|------|
| 案例一 | Li 在 Si 中的扩散 | 串行 DyNEB / 并行 NEB | 0.618 eV |
| 案例二 | 环己烷在 Pt@石墨烯上的 C-H 解离 | AutoNEB | 1.328 eV |
## 第一章：过渡态与 NEB 方法

### 1.1 过渡态与最小能量路径

过渡态理论将反应过程描述为势能面（PES）上的运动：反应物从一个局域极小值出发，越过一个鞍点（过渡态），到达另一个局域极小值（产物）。这条连接初末态、经过鞍点的路径称为**最小能量路径（MEP）**。

鞍点具有如下特征：
- 能量对所有坐标的一阶梯度为零
- 在反应坐标方向的二阶梯度为负，其他方向为正
- 振动频率分析中，**只在反应坐标方向存在唯一虚频**

Eyring 过渡态理论将反应速率常数与活化自由能联系起来：

$$k = \kappa \frac{k_B T}{h} \exp\!\left(\frac{-\Delta G^{\ddagger}}{RT}\right)$$

因此，找到过渡态结构、计算能垒，就能估算基元反应速率。

### 1.2 IT-NEB 与 CI-NEB

NEB（Nudged Elastic Band）方法是双端过渡态搜索中最常用的方法。其基本思路：

1. 在初态和末态之间生成若干**映像（images）**，映像之间通过弹簧力耦合，形成一条映像链
2. 对映像链进行优化，使其收敛到 MEP；链上能量最高的映像即为过渡态

在目前最常用的**IT-NEB（Improved Tangent NEB）**中，第 $i$ 个映像的切线方向定义为指向能量更高邻近映像的方向：

$$\boldsymbol{\tau}_i = \begin{cases} \boldsymbol{\tau}_i^+ & \text{if } E_{i+1} > E_i > E_{i-1} \\ \boldsymbol{\tau}_i^- & \text{if } E_{i+1} < E_i < E_{i-1} \end{cases}$$

在能量极值点附近做加权平均，保证切线的连续性。

**CI-NEB（Climbing Image NEB）** 在 IT-NEB 基础上，将能量最高点的映像的受力修改为：

$$\mathbf{F}_{i_{\max}} = -\nabla E(\mathbf{R}_{i_{\max}}) + 2\nabla E(\mathbf{R}_{i_{\max}})\big|_{\|}$$

即沿切线方向的力分量取反，驱使该映像爬升至鞍点。CI-NEB 不增加计算量，但显著提升过渡态定位精度。实际计算中，IT-NEB + CI-NEB 是标准组合。

**ASE 中对应的参数：**

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `algorism` | NEB 切线方法 | `"improvedtangent"`（IT-NEB） |
| `climb` | 是否启用 CI-NEB | `True` |
| `k` | 弹簧力常数（eV/Å²） | `0.10` |
| `fmax` | 收敛力阈值（eV/Å） | `0.05` |

优化器推荐使用 `FIRE`，它对 CI-NEB 的收敛表现优于 BFGS。

### 1.3 AutoNEB：动态映像机制

传统 NEB 计算从一开始就固定映像数，存在两个问题：
1. 各映像难度不均，靠近过渡态的映像电子结构难以自洽收敛，导致并行计算资源分配不均
2. 均匀插值的初始映像对过渡态附近的分辨率不足

AutoNEB 方法（Kolsbjerg et al., J. Chem. Phys. 145, 094107, 2016）通过**动态添加映像**解决上述问题，算法流程如下：

1. 从初末态 + 少量初猜出发，完成粗精度 NEB 计算
2. 定位 NEB 链上几何/能量差异最大的相邻映像对，插值添加一个新映像
3. 以新映像为中心继续局部 NEB 优化
4. 重复步骤 2–3，直到映像总数达到设定上限 `n_max`
5. 以能量最高点为中心，进行最终的 CI-NEB 计算

**AutoNEB 相对 NEB 的优势：**
- 过渡态附近的映像密度更高，定位更精确
- 每次只优化局部子链，计算量更小
- 对复杂体系（多原子、软势能面）的鲁棒性更好

**AutoNEB 特有参数：**

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `n_images` / `n_max` | 最终映像总数上限 | 8–12 |
| `n_simul` | 同时计算的映像数 | 等于并行进程数 |
| `fmax` | 两阶段收敛阈值 `[粗, 精]` | `[0.20, 0.05]` |
| `smooth_curve` | 最后额外平滑 NEB 链 | `False`（可选） |
## 第二章：ATST-Tools 工作流

### 2.1 环境配置

ATST-Tools 可从 GitHub 获取：

```bash
git clone https://github.com/QuantumMisaka/ATST-Tools
```

**依赖项：**

| 依赖 | 用途 |
|------|------|
| ASE-ABACUS 接口 | ABACUS 计算后端接入 |
| pymatgen + pymatgen-analysis-diffusion | IDPP 插值生成初猜 |
| GPAW | 并行 NEB / AutoNEB 所需的 MPI 环境 |
| Sella（可选） | 单端过渡态搜索 |

配置 Python 路径，使 source 目录下的库可被识别：

```bash
export PYTHONPATH=/path/to/ATST-Tools/source:$PYTHONPATH
```

将上述命令加入 `~/.bashrc` 即可永久生效。Bohrium 平台的推荐镜像（`atst-tools:1.5.0`）已预置上述环境。

### 2.2 目录结构

```
ATST-Tools/
├── neb/          # NEB 和 AutoNEB 计算脚本
│   ├── neb_make.py        # 生成初猜 NEB 链
│   ├── neb_run.py         # 串行/并行 NEB 计算
│   ├── autoneb_run.py     # AutoNEB 计算
│   ├── neb_post.py        # NEB 后处理
│   ├── neb_submit.sh      # Slurm 提交脚本（NEB）
│   └── autoneb_submit.sh  # Slurm 提交脚本（AutoNEB）
├── dimer/        # Dimer 方法脚本
├── sella/        # Sella 方法脚本
├── vibration/    # 振动分析脚本
│   └── vib_analysis.py
├── relax/        # ASE 结构优化脚本
└── source/       # 工作流核心库
    ├── abacus_neb.py      # AbacusNEB 类
    ├── abacus_autoneb.py  # AbacusAutoNEB 类
    └── neb2vib.py         # NEB 结果转振动分析
```

### 2.3 通用工作流

NEB 计算的标准三步流程：

```
步骤 1：neb_make.py   → 读取初末态，IDPP 插值 → init_neb_chain.traj
步骤 2：neb_run.py    → NEB 迭代优化           → neb.traj
        autoneb_run.py → AutoNEB 动态优化       → run_autoneb???.traj
步骤 3：neb_post.py   → 后处理                 → nebplots.pdf + neb_latest.traj
```

计算过程中，每个映像的 ABACUS 计算在独立子目录中运行：
- NEB：`NEBrun/NEB-rank{i}/`
- AutoNEB：`AutoNEBrun/AutoNEBrun_rank{i}/`

> **注意：** NEB 计算的初末态必须是已收敛的优化结构。建议用 ABACUS 的 `calculation = relax` 提前完成初末态弛豫，或使用 ATST-Tools 的 `relax/` 目录下的脚本。

> **表面计算注意：** 若体系为表面或六方体系，建议将真空层和 c 轴沿 y 方向设置，ABACUS 的 LCAO 计算在此方向效率更高。
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
## 第五章：振动分析与过渡态验证

### 5.1 为什么需要振动分析

NEB 计算给出能量最高的映像作为过渡态候选，但这只是能量上的判断。严格意义上的过渡态必须满足：**在反应坐标方向存在且仅存在一个虚频**。

振动分析（频率计算）的作用：
1. **验证过渡态正确性**：唯一虚频对应反应坐标方向的不稳定模式
2. **排查错误**：若虚频不止一个，说明过渡态结构不对
3. **零点能（ZPE）校正**：将电子能量转换为 0 K 下的振动包含能量
4. **有限温度自由能校正**：通过谐振近似计算振动熵和自由能

### 5.2 vib_analysis.py 配置

ATST-Tools 的 `vibration/vib_analysis.py` 使用 ASE 的有限差分方法，对指定原子施加微小位移，通过力的变化计算 Hessian 矩阵，进而得到振动模式。

**从 NEB 结果中提取过渡态和活跃原子：**

```python
from ase.io import read
from neb2vib import neb2vib

# 从 NEB 轨迹中自动识别过渡态结构和活跃原子
neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)
```

`neb2vib` 函数自动选取 NEB 链中位移最大的原子作为振动计算对象，无需手动指定原子序号。也可以手动指定：

```python
# 手动指定：只计算 H 原子（序号 0, 1）和 C 原子（序号 37）
vib_indices = [0, 1, 37]
```

完整脚本（针对案例二配置）：

```python
# vib_analysis.py — Cy-Pt@graphene 过渡态振动分析

from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world
from neb2vib import neb2vib
import os

# ===== 读取 NEB 结果 =====
neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)

T = 523.15  # K（振动自由能校正温度）

# ===== ABACUS 计算设置（与 NEB 保持一致）=====
abacus = "abacus"
mpi = 16
omp = 4
# ... pp, basis, parameters 与 autoneb_run.py 相同 ...

# ===== 振动计算参数 =====
vib_name = 'vib'
delta    = 0.01   # 有限位移步长（Å）
nfree    = 2      # 每个方向位移次数（正+负）

def set_calculator(abacus, parameters, mpi=1, omp=1):
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(f"mpirun -np {mpi} {abacus}")
    out_directory = f"SCF-rank{world.rank}"
    return Abacus(profile=profile, directory=out_directory, **parameters)

if __name__ == "__main__":
    atoms.calc = set_calculator(abacus, parameters, mpi=mpi, omp=omp)
    vib = Vibrations(atoms, indices=vib_indices,
                     name=vib_name, delta=delta, nfree=nfree)
    vib.run()
    vib.summary()
    vib.write_mode()

    # 热力学校正
    vib_energies = vib.get_energies()
    thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True)
    entropy     = thermo.get_entropy(T)
    free_energy = thermo.get_helmholtz_energy(T)
    print(f"Entropy:     {entropy:.6e} eV/K")
    print(f"Free Energy: {free_energy:.6f} eV")
```

运行：

```bash
python vib_analysis.py
```

### 5.3 计算结果解读

案例二过渡态的振动分析结果：

```
  #    meV     cm⁻¹
---------------------
  0   87.4i   705.0i   ← 唯一虚频：反应坐标方向（C-H 断裂）
  1    2.4     19.3
  2    3.5     28.3
  ...
 53  377.5   3044.9
---------------------
Zero-point energy: 4.416 eV
```

| 结果 | 数值 | 含义 |
|------|------|------|
| 虚频（imaginary） | 87.4i meV（705.0i cm⁻¹） | 对应 C-H 键断裂的反应坐标 |
| 虚频个数 | 1 | ✓ 过渡态正确（唯一虚频） |
| 零点能（ZPE） | 4.416 eV | 包含所有实振动模式的零点贡献 |
| 振动熵（523 K） | 见输出 | 用于 Helmholtz 自由能校正 |

**结果判读：**
- 虚频 705.0i cm⁻¹，对应标准 C-H 伸缩频率区间，物理上合理
- 只有 1 个虚频，确认该映像为真正的过渡态
- 查看 `vib.0.traj` 可在 ASE-GUI 中动画展示该虚频振动模式

```bash
ase -T gui vib.0.traj
```

计算完成后，`vib/` 目录中存储了所有位移结构的 DFT 结果（JSON 格式），下次修改温度重新分析时无需重新调用 ABACUS。
## 附录

### A. 辅助脚本

ATST-Tools 提供若干辅助脚本，简化日常操作：

| 脚本 | 用途 | 典型命令 |
|------|------|---------|
| `neb_dist.py` | 计算初末态原子间距，辅助确定映像数 | `python neb_dist.py INIT FINAL` |
| `traj_transform.py` | 将 traj 文件转换为 cif/extxyz/STRU 等格式 | `python traj_transform.py neb_latest.traj cif` |
| `traj_collect.py` | 从多个轨迹文件合并为一个 traj（续算用） | `python traj_collect.py ./AutoNEB_iter/run_autoneb???iter005.traj` |

### B. 续算方法

NEB 或 AutoNEB 因机器故障等中断后，可基于已有轨迹继续计算：

**NEB 续算：**
```bash
# 方法 1：从 neb.traj 提取最新状态
python neb_post.py neb.traj       # 生成 neb_latest.traj
python neb_make.py -i neb_latest.traj -n 3

# 方法 2：直接从中断的轨迹续算
python neb_make.py -i neb.traj -n 3
```

**AutoNEB 续算：**
```bash
# 收集某迭代阶段的轨迹（如第 5 步）
python traj_collect.py ./AutoNEB_iter/run_autoneb???iter005.traj
# 或收集最新的所有轨迹
python traj_collect.py ./run_autoneb???.traj
# 生成续算初猜
python neb_make.py -i collection.traj -n 4
```

> `n_max` 的值等于映像总数减 2（不含初末态）。

### C. 参考资料

1. ATST-Tools 主页：https://github.com/QuantumMisaka/ATST-Tools
2. ASE NEB 文档：https://wiki.fysik.dtu.dk/ase/ase/neb.html
3. AutoNEB 原始论文：E. L. Kolsbjerg, M. N. Groves, B. Hammer, *J. Chem. Phys.* 145, 094107 (2016)
4. ABACUS 官方文档：https://abacus.deepmodeling.com
5. ASE-ABACUS 接口：https://gitlab.com/1041176461/ase-abacus
