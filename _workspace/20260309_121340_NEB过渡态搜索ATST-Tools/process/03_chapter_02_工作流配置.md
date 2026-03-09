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
