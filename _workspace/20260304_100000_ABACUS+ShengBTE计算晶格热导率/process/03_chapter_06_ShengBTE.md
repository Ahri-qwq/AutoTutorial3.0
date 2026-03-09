# 第六章：Step 3 — 运行 ShengBTE

进入 `shengbte/` 目录，确认以下三个文件已就位：
- `FORCE_CONSTANTS_2ND`（来自 Step 1）
- `FORCE_CONSTANTS_3RD`（来自 Step 2）
- `CONTROL`（ShengBTE 参数文件，已提供）

## 6.1 CONTROL 文件

`CONTROL` 是 ShengBTE 的核心参数文件，使用 Fortran namelist 格式，包含四个段：

```fortran
&allocations
    nelements=1
    natoms=2
    ngrid(:)=10 10 10
&end
&crystal
    lfactor=0.100000
    lattvec(:,1)=0 2.81594778072 2.81594778072
    lattvec(:,2)=2.81594778072 0 2.81594778072
    lattvec(:,3)=2.81594778072 2.81594778072 0
    elements="Si"
    types=1 1
    positions(:,1)=0.8750000000000000  0.8750000000000000  0.8750000000000000
    positions(:,2)=0.1250000000000000  0.1250000000000000  0.1250000000000000
    scell(:)=2 2 2
&end
&parameters
    !T=300,
    T_min=200
    T_max=500
    T_step=50
    scalebroad=1.0
&end
&flags
    !espresso=.true.
    nonanalytic=.true.,
    isotopes=.true.
&end
```

各段参数说明：

**&allocations**
| 参数 | 值 | 说明 |
|------|----|------|
| `nelements` | 1 | 元素种类数（Si 只有 1 种） |
| `natoms` | 2 | 原胞中原子数 |
| `ngrid(:)` | 10 10 10 | 布里渊区 q 点积分网格，越大精度越高 |

**&crystal**
| 参数 | 说明 |
|------|------|
| `lfactor` | 格矢量长度单位因子（0.1 nm = 1 Å） |
| `lattvec(:,i)` | 三个格矢量，单位为 lfactor × nm = Å，与 STRU 中 LATTICE_VECTORS 对应 |
| `elements` | 元素名称 |
| `types` | 每个原子的元素类型编号（两个 Si 均为类型 1） |
| `positions(:,i)` | 每个原子的分数坐标，与 STRU 中的原子坐标一致 |
| `scell(:)` | 计算力常数时使用的超胞大小，**必须与 Step 2 保持一致** |

**&parameters**
| 参数 | 值 | 说明 |
|------|----|------|
| `T_min` | 200 | 计算温度下限（K） |
| `T_max` | 500 | 计算温度上限（K） |
| `T_step` | 50 | 温度步长（K），本例计算 200、250、300、350、400、450、500 K |
| `scalebroad` | 1.0 | 声子线宽展宽因子，控制散射积分的展宽精度 |

（`!T=300` 行以感叹号开头表示注释，ShengBTE 不会读取该行。）

**&flags**
| 参数 | 值 | 说明 |
|------|----|------|
| `nonanalytic` | .true. | 开启非解析修正（对极性材料有较大影响，Si 中效果较小） |
| `isotopes` | .true. | 考虑天然同位素散射 |

## 6.2 运行 ShengBTE

```bash
mpirun -n 10 ShengBTE
```

ShengBTE 运行期间会输出声子对称性信息、q 点数量等日志，例如：

```
Info: symmetry group Fd-3m detected
Info: 48 symmetry operations
Info: Ntot = 1000
Info: Nlist = 47
Info: about to obtain the spectrum
Info: expecting Phonopy 2nd-order format
```

计算完成后，热导率结果保存在 `BTE.kappa_*` 系列文件中。

## 6.3 计算结果

本教程参数（2×2×2 超胞，2×2×2 k 点，ngrid=10×10×10）下，Si 的计算结果如下：

| 基组 | 300 K 热导率 | 实验值（300 K） |
|------|-------------|----------------|
| LCAO | ~100 W/(m·K) | ~150 W/(m·K) |
| PW   | ~100 W/(m·K) | ~150 W/(m·K) |

LCAO 和 PW 基组的结果高度一致，说明两种基组对 Si 体系均能给出可靠的力常数。

计算值低于实验值约 30%，原因是本教程为演示目的使用了较小的超胞（2×2×2）和较稀的 k 点（2×2×2），导致力常数截断误差和 k 点采样误差同时存在。实际科研中需要系统地测试这两个参数（见第八章）。

对于 PW 基组的计算，流程与 LCAO 完全相同，唯一差别是三阶力常数 SCF 中的 `scf_thr` 需要设为至少 `1e-12`。
