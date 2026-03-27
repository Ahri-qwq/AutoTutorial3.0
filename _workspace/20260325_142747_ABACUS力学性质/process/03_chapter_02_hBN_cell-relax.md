# 第二章：案例一——h-BN 变胞优化

六方氮化硼（h-BN）是典型的层状材料，属六方晶系，结构类似石墨。层内 B-N 键强，层间为弱 van der Waals 相互作用。这使 a/b 方向与 c 方向的弹性性质差异显著。

层间距离对 c 方向应力非常敏感，cell-relax 需要同时收紧力和应力才能得到可靠的平衡结构。本案例以此体系演示变胞优化的完整参数配置。

本案例体系：**h-BN 192原子超胞**（B₉₆N₉₆），LCAO 基组，PBE 泛函。

---

## 2.1 INPUT 文件

```
# h-BN cell-relax INPUT
INPUT_PARAMETERS
suffix                  BN
pseudo_dir              ./
orbital_dir             ./
calculation             cell-relax
symmetry                0
basis_type              lcao
ecutwfc                 100
scf_thr                 1e-07
scf_nmax                100
smearing_method         gauss
smearing_sigma          0.002
mixing_type             pulay
mixing_beta             0.3
cal_force               1
cal_stress              1
force_thr_ev            0.01
stress_thr              0.5
relax_nmax              100
out_stru                1
kspacing                0.08
```

参数说明：

| 参数 | 值 | 说明 |
|------|------|------|
| `calculation` | `cell-relax` | 同时优化晶胞和原子位置 |
| `symmetry` | `0` | 关闭对称性。层状结构在弛豫过程中可能降低对称性，关闭后更安全 |
| `cal_force` | `1` | 开启受力计算，cell-relax 必须 |
| `cal_stress` | `1` | 开启应力张量计算，cell-relax 必须 |
| `force_thr_ev` | `0.01` | 力收敛阈值，0.01 eV/Å 适用于大多数材料 |
| `stress_thr` | `0.5` | 应力收敛阈值，0.5 kBar 较严格，弹性计算前建议使用 |
| `kspacing` | `0.08` | 自动生成 K 点网格，代替手写 KPT 文件，单位 Å⁻¹ |
| `smearing_method` | `gauss` | 高斯展宽，处理绝缘体时展宽量可小 |
| `smearing_sigma` | `0.002` | 展宽参数（Ry），h-BN 是绝缘体，取小值 |
| `mixing_type` | `pulay` | Pulay 混合，适合 LCAO 计算 |
| `mixing_beta` | `0.3` | 混合系数，大体系取 0.3–0.4 合适 |

> **注意：** `scf_thr` 设为 `1e-7` 比默认值严格，确保每一步 SCF 充分收敛，减少应力计算误差。结构优化的精度上限由 SCF 精度决定。

---

## 2.2 STRU 文件

h-BN 192原子体系的 STRU 关键部分如下：

```
ATOMIC_SPECIES
B   10.81   B_ONCV_PBE-1.0.upf
N   14.007  N_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
B_gga_8au_100Ry_2s2p1d.orb
N_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
2.6665300000  0.0000000000  0.0000000000
0.0000000000  16.970664000  0.0000000000
0.0000000000  0.0000000000  28.0218

ATOMIC_POSITIONS
Direct
B
0.0
96
...（96个B原子坐标）
N
0.0
96
...（96个N原子坐标）
```

晶格常数单位为 Bohr（`LATTICE_CONSTANT` 给出 Bohr 换算因子），`LATTICE_VECTORS` 给出无量纲晶格矢量（乘以 `LATTICE_CONSTANT` 后为 Bohr 单位）。

轨道文件说明：
- `B_gga_8au_100Ry_2s2p1d.orb`：B 的 DZP 轨道，截断半径 8 a.u.，截断能 100 Ry
- `N_gga_8au_100Ry_2s2p1d.orb`：N 的 DZP 轨道，截断半径 8 a.u.，截断能 100 Ry
- 两者均与 `ONCV_PBE-1.0` 系列赝势匹配

> 完整的 192原子坐标文件可从 ABACUS 官方 GitHub 仓库的 examples 目录获取。

---

## 2.3 运行与收敛过程

提交计算后，ABACUS 的输出日志（`running_cell-relax.log`）中会依次输出每轮 RELAX CELL 和 RELAX IONS 的信息（以下为示意性输出，实际数值因体系而异）：

```
-------------------------------------------
RELAX CELL : 1
RELAX IONS : 1 (in total: ...)
-------------------------------------------
...
  LARGEST GRAD (eV/A)  :      0.312
  TOTAL-PRESSURE: -2.070e+00 KBAR     <-- 绝对值 2.07 > 0.5，未收敛
-------------------------------------------
RELAX CELL : 2
RELAX IONS : 1 (in total: ...)
-------------------------------------------
...
  LARGEST GRAD (eV/A)  :      0.054
  TOTAL-PRESSURE: -8.500e-01 KBAR     <-- 绝对值 0.85 > 0.5，未收敛
-------------------------------------------
RELAX CELL : 3
RELAX IONS : 1 (in total: ...)
-------------------------------------------
...
  LARGEST GRAD (eV/A)  :      0.008   <-- 小于 force_thr_ev=0.01，满足
  TOTAL-PRESSURE: -4.350e-01 KBAR     <-- 绝对值 0.435 < 0.5，满足 → 收敛！
```

判断收敛的两个关键量：

- **LARGEST GRAD**：所有原子受力分量中的最大值（对应 `force_thr_ev`）
- **TOTAL-PRESSURE**：应力张量对角元均值（对应 `stress_thr`）

当两者**同时**低于设定阈值时，计算收敛，输出最终优化结构 `OUT.BN/STRU_ION_D`（通过 `out_stru 1` 开启）。

---

## 2.4 常见问题

**Q：cell-relax 跑了很多步但应力一直不降？**

可能原因：初始结构残余应力太大，或 `kspacing` 取值偏大（K 点密度不够）。建议检查 K 点收敛性，或适当减小 `kspacing`（如从 0.1 改为 0.08）。

**Q：`symmetry 0` 和 `symmetry 1` 有什么区别？**

开启对称性（`symmetry 1`）会约束原子位置使其满足对称性，减少独立自由度，加速收敛，但有时会因程序识别对称性不准确导致弛豫出错。大体系或低对称性材料建议用 `symmetry 0`。

**Q：弛豫完成后如何继续进行弹性计算？**

将 `OUT.BN/STRU_ION_D` 复制为新目录的 `STRU`，修改 `INPUT` 中 `calculation scf`（或交给 abacustest 自动处理），即可在优化后结构上计算弹性常数。
