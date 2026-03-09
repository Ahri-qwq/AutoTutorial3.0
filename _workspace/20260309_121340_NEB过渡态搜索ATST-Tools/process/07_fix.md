---
title: "使用 ABACUS 和 ATST-Tools 进行 NEB/AutoNEB 过渡态搜索"
author: "AutoTutorial 3.0"
date: "2026-03-09"
topic: "NEB AutoNEB 过渡态搜索 ATST-Tools"
task_type: "C"
has_case: true
word_count: ~5500
chapters: 5
---

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

---

## 第一章：过渡态与 NEB 方法

### 1.1 过渡态与最小能量路径

过渡态理论将反应过程描述为势能面（PES）上的运动：反应物从一个局域极小值出发，越过一个鞍点（过渡态），到达另一个局域极小值（产物）。这条连接初末态、经过鞍点的路径称为**最小能量路径（MEP）**。

鞍点具有如下特征：
- 能量对所有坐标的一阶梯度为零
- 在反应坐标方向的二阶梯度为负，其他方向为正
- 振动频率分析中，**只在反应坐标方向存在唯一虚频**

Eyring 过渡态理论将反应速率常数与活化自由能联系起来：

$$k = \kappa \frac{k_B T}{h} \exp\!\left(\frac{-\Delta G^{\ddagger}}{RT}\right)$$

找到过渡态结构、计算能垒，就能估算基元反应速率。

### 1.2 IT-NEB 与 CI-NEB

NEB（Nudged Elastic Band）方法是双端过渡态搜索中最常用的方法：

1. 在初态和末态之间生成若干**映像（images）**，映像之间通过弹簧力耦合，形成一条映像链
2. 对映像链进行优化，使其收敛到 MEP；链上能量最高的映像即为过渡态

在 **IT-NEB（Improved Tangent NEB）** 中，第 $i$ 个映像的切线方向定义为指向能量更高邻近映像的方向：

$$\boldsymbol{\tau}_i = \begin{cases} \boldsymbol{\tau}_i^+ & \text{if } E_{i+1} > E_i > E_{i-1} \\ \boldsymbol{\tau}_i^- & \text{if } E_{i+1} < E_i < E_{i-1} \end{cases}$$

在能量极值点附近做加权平均，保证切线的连续性。

**CI-NEB（Climbing Image NEB）** 在 IT-NEB 基础上，将能量最高点映像的受力修改为：

$$\mathbf{F}_{i_{\max}} = -\nabla E(\mathbf{R}_{i_{\max}}) + 2\nabla E(\mathbf{R}_{i_{\max}})\big|_{\|}$$

即沿切线方向的力分量取反，驱使该映像爬升至鞍点。CI-NEB 不增加计算量，但显著提升过渡态定位精度。实际计算中，IT-NEB + CI-NEB 是标准组合，优化器推荐使用 `FIRE`。

**NEB 关键参数：**

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `algorism` | NEB 切线方法 | `"improvedtangent"`（IT-NEB） |
| `climb` | 是否启用 CI-NEB | `True` |
| `k` | 弹簧力常数（eV/Å²） | `0.10` |
| `fmax` | 收敛力阈值（eV/Å） | `0.05` |

### 1.3 AutoNEB：动态映像机制

传统 NEB 从一开始就固定映像数，存在两个问题：
1. 靠近过渡态的映像电子结构难以自洽收敛，并行计算资源分配不均
2. 均匀插值的初始映像对过渡态附近的分辨率不足

AutoNEB 方法（Kolsbjerg et al., *J. Chem. Phys.* 145, 094107, 2016）通过**动态添加映像**解决上述问题：

1. 从初末态 + 少量初猜出发，完成粗精度 NEB 计算
2. 定位链上几何/能量差异最大的相邻映像对，插值添加新映像
3. 以新映像为中心继续局部 NEB 优化
4. 重复步骤 2–3，直到映像总数达到上限 `n_max`
5. 以能量最高点为中心，进行最终 CI-NEB 计算

**AutoNEB 特有参数：**

| 参数 | 含义 | 推荐值 |
|------|------|--------|
| `n_images` / `n_max` | 最终映像总数上限 | 8–12 |
| `n_simul` | 同时计算的映像数 | 等于并行进程数 |
| `fmax` | 两阶段收敛阈值 `[粗精度, 最终CI]` | `[0.20, 0.05]` |
| `smooth_curve` | 最终额外平滑 NEB 链 | `False`（可选） |

---

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

配置 Python 路径：

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
步骤 1：neb_make.py    → 读取初末态，IDPP 插值   → init_neb_chain.traj
步骤 2：neb_run.py     → NEB 迭代优化             → neb.traj
        autoneb_run.py → AutoNEB 动态优化          → run_autoneb???.traj
步骤 3：neb_post.py    → 后处理，输出曲线与轨迹   → nebplots.pdf + neb_latest.traj
```

计算过程中，每个映像的 ABACUS 计算在独立子目录中运行：
- NEB：`NEBrun/NEB-rank{i}/`
- AutoNEB：`AutoNEBrun/AutoNEBrun_rank{i}/`

> **注意：** NEB 计算的初末态必须是已收敛的优化结构，建议用 ABACUS 的 `calculation = relax` 提前弛豫，或使用 ATST-Tools `relax/` 目录下的脚本。

> **表面计算：** 若体系含真空层或为六方晶系，建议将真空层和 c 轴沿 y 方向设置，ABACUS LCAO 计算在此方向效率更高。

---

## 第三章：案例一——Li 在 Si 中的扩散

### 3.1 体系简介

本案例计算 Li 原子在 Si 晶胞中从一个间隙位扩散到相邻间隙位的过程：
- 只涉及单个 Li 原子的迁移，映像间原子对应关系清晰
- 初末态对称，能垒约 0.6 eV，NEB 链收敛快
- 适合作为 NEB 入门案例

工作目录结构：

```
Li-diffu-Si/
├── INIT/                         # 初态计算结果（已优化）
│   └── OUT.ABACUS/running_scf.log
├── FINAL/                        # 末态计算结果（已优化）
│   └── OUT.ABACUS/running_scf.log
├── Li_ONCV_PBE-1.2.upf
├── Si_ONCV_PBE-1.2.upf
├── Li_gga_8au_100Ry_4s1p.orb
└── Si_gga_8au_100Ry_2s2p1d.orb
```

### 3.2 生成初猜 NEB 链

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
| `-m` | 插值方法，默认 `IDPP`（Pymatgen 版本，处理复杂体系更稳健） |
| `--fix` | 固定某一高度以下的原子（`[高度]:[方向]`，用于表面计算） |
| `--mag` | 设置初始磁矩（`元素:磁矩,...`） |

执行后输出 `init_neb_chain.traj`，共 5 个映像（初态 + 3 个中间态 + 末态）：

```
Reading files: ./INIT/OUT.ABACUS/running_scf.log and ./FINAL/OUT.ABACUS/running_scf.log
Generating path, number of images: 3, sort_tol: 1.0
Optimizing path using IDPP method
Writing path: init_neb_chain.traj, Number of images: 5
```

> 可先用 `neb_dist.py` 查看初末态原子间距来估计合适的映像数。简单扩散取 3 个映像通常足够；复杂表面反应建议 5–8 个。

### 3.3 串行 DyNEB 计算

将 `neb_run.py` 复制到工作目录，按以下配置编辑（`dyneb_run.py`）：

```python
# dyneb_run.py — Li-diffu-Si 串行 DyNEB

from ase.optimize import FIRE
from ase.io import read
from abacus_neb import AbacusNEB

# ===== NEB 优化设置 =====
mpi           = 16           # ABACUS 每个映像使用的 MPI 进程数
omp           = 1            # OpenMP 线程数
fmax          = 0.05         # 收敛力阈值（eV/Å）
neb_optimizer = FIRE
neb_directory = "NEBrun"
algorism      = "improvedtangent"  # IT-NEB
climb         = True               # 启用 CI-NEB
dyneb         = True               # 串行 DyNEB 模式
parallel      = False              # 串行，不使用 GPAW MPI
k             = 0.10               # 弹簧力常数（eV/Å²）
init_chain    = "init_neb_chain.traj"
abacus        = "abacus"

# ===== 赝势与轨道（与脚本同目录时 lib_dir 留空）=====
lib_dir    = ""
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
    'calculation':     'scf',
    'nspin':           1,
    'xc':              'pbe',
    'ecutwfc':         100,
    'dft_functional':  'pbe',
    'ks_solver':       'genelpa',
    'symmetry':        0,          # NEB 必须关闭对称性
    'vdw_method':      'none',
    'smearing_method': 'gaussian',
    'smearing_sigma':  0.001,
    'basis_type':      'lcao',
    'mixing_type':     'broyden',
    'scf_thr':         1e-6,
    'scf_nmax':        100,
    'kpts':            kpts,
    'pp':              pp,
    'basis':           basis,
    'pseudo_dir':      pseudo_dir,
    'basis_dir':       basis_dir,
    'init_wfc':        'atomic',
    'init_chg':        'atomic',
    'cal_force':       1,          # NEB 优化需要原子受力
    'cal_stress':      1,
    'out_stru':        1,
    'out_chg':         -1,         # 覆盖保存，节省磁盘空间
}

if __name__ == "__main__":
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters,
                    parallel=parallel, directory=neb_directory,
                    mpi=mpi, omp=omp, abacus=abacus,
                    algorism=algorism, k=k, dyneb=dyneb)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)
```

运行：

```bash
python dyneb_run.py
```

ASE 打印 FIRE 优化过程（每步约 6 分钟/映像，16 核机器）：

```
Notice: Parallel calculation is not being used
Set NEB method to DyNEB automatically
----- Running Dynamic NEB -----
----- improvedtangent method is being used -----
      Step     Time          Energy          fmax
FIRE:    0 18:28:24    -7051.845345         0.733359
FIRE:    1 18:34:09    -7051.879062         0.594320
FIRE:    2 18:39:55    -7051.922971         0.350617
...
FIRE:   23 20:28:27    -7051.985587         0.036730
```

本案例在 32 核机器上约需 1–2 小时收敛。

### 3.4 并行 NEB

并行 NEB 将每个映像分配给独立进程，多映像同时计算。修改以下参数：

```python
mpi      = 5      # 每个映像内 ABACUS 使用的核数
dyneb    = False  # 并行模式不使用 DyNEB
parallel = True   # 启用并行
```

并行启动：

```bash
mpirun -np 3 gpaw python neb_run.py
```

> `-np` 进程数应等于中间映像数（本例为 3）。ABACUS 和 GPAW 需依赖同一套 MPI 环境。

并行 NEB 与串行 DyNEB 结果一致，但挂钟时间显著缩短。

### 3.5 后处理与结果

```bash
python neb_post.py neb.traj
```

输出：
- `nebplots_chain.pdf`：势能–反应坐标曲线（含各映像切线力投影）
- `neb_latest.traj`：收敛时的 NEB 链

各映像能量与结果：

```
num: 0; Energy: -7052.6037618 (eV)
num: 1; Energy: -7052.2767038 (eV)
num: 2; Energy: -7051.9855868 (eV)    ← 过渡态（能量最高）
num: 3; Energy: -7052.2797724 (eV)
num: 4; Energy: -7052.6035978 (eV)
Reaction Barrier and Energy Difference: (0.618 eV, 0.0002 eV)
```

| 结果 | 数值 |
|------|------|
| Li 扩散能垒 | 0.618 eV |
| 反应能差 | 0.0002 eV（近似对称路径） |

可视化：

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

---

## 第四章：案例二——环己烷在 Pt@石墨烯上的 AutoNEB

### 4.1 体系简介与方法选择

本案例计算环己烷（C₆H₁₂）在 Pt 原子负载石墨烯表面的 C-H 键断裂反应。

这一体系比案例一复杂得多：原子数多，势能面较"软"，且涉及 Pt 的 d 轨道和石墨烯 π 电子，电子结构自洽收敛较慢。实验表明，直接用 IT-NEB（8 个固定映像）计算会**收敛到错误的过渡态**（能量路径异常，切线力在最高点不接近零）。

AutoNEB 通过动态调整映像分布，使过渡态附近的映像密度更高，从而得到正确的过渡态（能垒 1.328 eV）。

### 4.2 初末态准备与插值

工作目录结构：

```
Cy-Pt@graphene/
├── IS/                      # 初态（已优化）
│   └── OUT.ABACUS/running_relax.log
├── FS/                      # 末态（已优化）
│   └── OUT.ABACUS/running_relax.log
├── C_ONCV_PBE-1.0.upf
├── H_ONCV_PBE-1.0.upf
├── Pt_ONCV_PBE-1.0.upf
├── C_gga_8au_100Ry_2s2p1d.orb
├── H_gga_8au_100Ry_2s1p.orb
└── Pt_gga_7au_100Ry_4s2p2d1f.orb
```

生成初猜（取 4 个中间映像，与 AutoNEB 的 `n_simul` 对应）：

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

> Pymatgen IDPP 在无法自动匹配原子时会发出 `Auto sorting is turned off` 警告，通常不影响计算，AutoNEB 会在优化过程中自动修正映像位置。

### 4.3 AutoNEB 脚本配置

将 `autoneb_run.py` 复制到工作目录并编辑：

```python
# autoneb_run.py — Cy-Pt@graphene AutoNEB

from ase.optimize import FIRE
from ase.io import read
from ase.parallel import world
from abacus_autoneb import AbacusAutoNEB

# ===== AutoNEB 设置 =====
mpi           = 16
omp           = 4
neb_optimizer = FIRE
neb_directory = "AutoNEBrun"
algorism      = "improvedtangent"
init_chain    = "init_neb_chain.traj"
climb         = True
fmax          = [0.20, 0.05]   # 两阶段收敛阈值：粗精度 → 最终 CI-NEB
n_simul       = world.size     # 同时计算的映像数 = MPI 进程数
n_images      = 10             # 最终映像总数上限
smooth_curve  = False
k             = 0.10

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
    'C':  'C_gga_8au_100Ry_2s2p1d.orb',
    'H':  'H_gga_8au_100Ry_2s1p.orb',
    'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb',
}

# ===== ABACUS INPUT 参数 =====
kpts = [2, 1, 2]
parameters = {
    'calculation':     'scf',
    'nspin':           2,           # 含 Pt，开启自旋极化
    'xc':              'pbe',
    'ecutwfc':         100,
    'dft_functional':  'pbe',
    'ks_solver':       'genelpa',
    'symmetry':        0,
    'vdw_method':      'd3_bj',     # 表面吸附体系推荐色散修正
    'smearing_method': 'gaussian',
    'smearing_sigma':  0.001,
    'basis_type':      'lcao',
    'mixing_type':     'broyden',
    'scf_thr':         1e-6,
    'scf_nmax':        100,
    'kpts':            kpts,
    'pp':              pp,
    'basis':           basis,
    'pseudo_dir':      pseudo_dir,
    'basis_dir':       basis_dir,
    'init_wfc':        'atomic',
    'init_chg':        'atomic',
    'cal_force':       1,
    'cal_stress':      1,
    'out_stru':        1,
    'out_chg':         0,
    'out_mul':         0,
    'out_wfc_lcao':    0,
    'out_bandgap':     0,
    'efield_flag':     1,           # 表面偶极修正
    'dip_cor_flag':    1,
    'efield_dir':      1,           # 修正方向 y（真空层方向）
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

**与案例一 neb_run.py 的关键差异：**

| 参数 | 案例一（NEB） | 案例二（AutoNEB） | 说明 |
|------|-------------|-----------------|------|
| `fmax` | `0.05` | `[0.20, 0.05]` | AutoNEB 分两阶段收敛 |
| `n_images` | — | `10` | 最终映像总数上限 |
| `n_simul` | — | `world.size` | 同时优化的映像数 |
| `nspin` | `1` | `2` | 含 Pt，需自旋极化 |
| `vdw_method` | `none` | `d3_bj` | 表面吸附需要色散修正 |
| `efield_flag` | `0` | `1` | 表面偶极修正 |

并行启动（4 进程对应 4 个同时计算的映像）：

```bash
mpirun -np 4 gpaw python autoneb_run.py
```

### 4.4 AutoNEB 迭代过程

AutoNEB 的运行日志（`running_autoneb.out`）记录动态映像添加过程：

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

| 迭代 | 内容 | 映像总数 |
|------|------|---------|
| iter 1 | 评估初始 6 个映像 | 6 |
| iter 2–5 | 动态添加映像（聚焦于差异最大处） | 7 → 10 |
| iter 6 | 最终 CI-NEB，精细收敛过渡态 | 10 |

最终 CI-NEB（iter 6）收敛过程：

```
      Step     Time          Energy          fmax
FIRE:    0 01:57:54   -11866.084572         0.652487
...
FIRE:   88 03:06:42   -11865.497000         0.040039
```

### 4.5 后处理与结果

```bash
# AutoNEB 后处理（Shell 自动按序展开文件名）
python neb_post.py --autoneb run_autoneb???.traj
```

各映像能量：

```
num: 0; Energy: -11866.824886 (eV)
num: 1; Energy: -11866.822720 (eV)
num: 2; Energy: -11866.801706 (eV)
num: 3; Energy: -11866.782620 (eV)
num: 4; Energy: -11866.587954 (eV)
num: 5; Energy: -11865.497000 (eV)    ← 过渡态（CI 映像）
num: 6; Energy: -11866.330196 (eV)
num: 7; Energy: -11866.396832 (eV)
num: 8; Energy: -11866.429487 (eV)
num: 9; Energy: -11866.435419 (eV)
Reaction Barrier and Energy Difference: (1.328 eV, 0.389 eV)
```

| 结果 | 数值 |
|------|------|
| C-H 断裂能垒 | 1.328 eV |
| 反应能差（H 解离吸热量） | 0.389 eV |

**判断过渡态是否正确：** 在 `nebplots_all.pdf` 的势能曲线中，过渡态处的切线力投影应接近零（切线近乎水平）。最终通过第五章的振动分析确认唯一虚频。

实时监控 AutoNEB 进度：

```bash
ase -T gui -n -1 run_autoneb???.traj
```

---

## 第五章：振动分析与过渡态验证

### 5.1 唯一虚频判据

NEB 给出能量最高的映像作为过渡态候选，严格验证还需要振动分析：**真正的过渡态只在反应坐标方向有一个虚频**。

振动分析的作用：
1. **验证过渡态**：唯一虚频对应反应坐标方向的不稳定振动模式
2. **零点能（ZPE）校正**：将电子能量修正为包含零点振动的能量
3. **有限温度自由能校正**：通过谐振近似计算振动自由能

### 5.2 vib_analysis.py 配置

ATST-Tools 使用 ASE 的有限差分方法：对选定原子施加微小位移（默认 0.01 Å），通过力的变化计算 Hessian，进而得到振动模式。

从 NEB 轨迹中自动提取过渡态和活跃原子：

```python
from ase.io import read
from neb2vib import neb2vib

neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)
# neb2vib 自动选取 NEB 链中位移最大的原子作为振动计算对象
```

也可手动指定：

```python
# 只对参与反应的原子做振动计算（大幅减少计算量）
vib_indices = [atom.index for atom in atoms if atom.symbol in ('C', 'H')]
```

完整脚本（案例二配置，ABACUS 参数与 autoneb_run.py 保持一致）：

```python
# vib_analysis.py — Cy-Pt@graphene 过渡态振动分析

import os
import numpy as np
from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world
from neb2vib import neb2vib

# ===== 读取 NEB 结果 =====
neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)

T = 523.15  # K

# ===== ABACUS 设置（与 autoneb_run.py 保持一致）=====
abacus = "abacus"
mpi, omp = 16, 4
pp    = {'H': 'H_ONCV_PBE-1.0.upf',
         'C': 'C_ONCV_PBE-1.0.upf',
         'Pt': 'Pt_ONCV_PBE-1.0.upf'}
basis = {'H': 'H_gga_8au_100Ry_2s1p.orb',
         'C': 'C_gga_8au_100Ry_2s2p1d.orb',
         'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb'}
parameters = {
    'calculation': 'scf', 'nspin': 2, 'xc': 'pbe', 'ecutwfc': 100,
    'ks_solver': 'genelpa', 'symmetry': 0, 'vdw_method': 'd3_bj',
    'smearing_method': 'gaussian', 'smearing_sigma': 0.001,
    'basis_type': 'lcao', 'mixing_type': 'broyden',
    'scf_thr': 1e-6, 'scf_nmax': 100,
    'kpts': [2, 1, 2], 'pp': pp, 'basis': basis,
    'pseudo_dir': '', 'basis_dir': '',
    'init_wfc': 'atomic', 'init_chg': 'file',
    'cal_force': 1, 'cal_stress': 1,
    'out_stru': 1, 'out_chg': 1,
    'efield_flag': 1, 'dip_cor_flag': 1, 'efield_dir': 1,
}

# ===== 振动计算参数 =====
vib_name = 'vib'
delta    = 0.01   # 有限位移步长（Å）
nfree    = 2      # 每个方向正负各一次

def set_calculator(abacus, parameters, mpi=1, omp=1):
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(f"mpirun -np {mpi} {abacus}")
    return Abacus(profile=profile,
                  directory=f"SCF-rank{world.rank}",
                  **parameters)

if __name__ == "__main__":
    atoms.calc = set_calculator(abacus, parameters, mpi=mpi, omp=omp)
    vib = Vibrations(atoms, indices=vib_indices,
                     name=vib_name, delta=delta, nfree=nfree)
    vib.run()
    vib.summary()
    vib.write_mode()   # 输出各振动模式的轨迹动画

    # 热力学校正
    vib_energies = vib.get_energies()
    thermo       = HarmonicThermo(vib_energies, ignore_imag_modes=True)
    entropy      = thermo.get_entropy(T)
    free_energy  = thermo.get_helmholtz_energy(T)
    print(f"Entropy:     {entropy:.6e} eV/K")
    print(f"Free Energy: {free_energy:.6f} eV")
```

> `init_chg = 'file'` 使振动计算复用 NEB 最后一步的电荷密度，加速 SCF 收敛。

运行：

```bash
python vib_analysis.py
```

### 5.3 结果解读

案例二过渡态的振动分析结果（`running_vib.out`）：

```
  #    meV     cm⁻¹
---------------------
  0   87.4i   705.0i    ← 唯一虚频
  1    2.4     19.3
  2    3.5     28.3
  3    6.1     48.9
  ...
 53  377.5   3044.9
---------------------
Zero-point energy: 4.416 eV
```

| 结果 | 数值 | 含义 |
|------|------|------|
| 虚频 | 87.4i meV（705.0i cm⁻¹） | 对应 C-H 键断裂的反应坐标 |
| 虚频个数 | 1 | ✓ 唯一虚频，过渡态正确 |
| ZPE | 4.416 eV | 零点振动能（含所有实振动模式） |

705 cm⁻¹ 处于 C-H 弯曲振动区间，物理上合理。

查看该虚频振动模式的动画：

```bash
ase -T gui vib.0.traj
```

计算完成后，`vib/` 目录中以 JSON 格式缓存了所有位移结构的 DFT 结果，修改温度重新进行热力学分析时无需重复调用 ABACUS。

---

## 附录

### A. 辅助脚本

| 脚本 | 用途 | 典型用法 |
|------|------|---------|
| `neb_dist.py` | 查看初末态原子间距，辅助估计映像数 | `python neb_dist.py INIT FINAL` |
| `traj_transform.py` | 将 traj 转换为 cif/extxyz/STRU 等格式 | `python traj_transform.py neb_latest.traj cif` |
| `traj_collect.py` | 从多个轨迹合并为单个 traj（续算用） | `python traj_collect.py ./AutoNEB_iter/run_autoneb???iter005.traj` |

### B. 续算方法

NEB 或 AutoNEB 因中断需要续算时：

**NEB 续算：**
```bash
python neb_post.py neb.traj          # 生成 neb_latest.traj
python neb_make.py -i neb_latest.traj -n 3   # 重新生成初猜
```

**AutoNEB 续算：**
```bash
# 收集某迭代阶段的轨迹
python traj_collect.py ./AutoNEB_iter/run_autoneb???iter005.traj
# 或收集最新状态
python traj_collect.py ./run_autoneb???.traj
# 生成续算初猜（n = 映像总数 - 2）
python neb_make.py -i collection.traj -n 4
```

### C. 参考资料

1. ATST-Tools 主页：https://github.com/QuantumMisaka/ATST-Tools
2. ASE NEB 文档：https://wiki.fysik.dtu.dk/ase/ase/neb.html
3. AutoNEB 原始论文：E. L. Kolsbjerg, M. N. Groves, B. Hammer, *J. Chem. Phys.* 145, 094107 (2016)
4. ABACUS 官方文档：https://abacus.deepmodeling.com
5. ASE-ABACUS 接口：https://gitlab.com/1041176461/ase-abacus
