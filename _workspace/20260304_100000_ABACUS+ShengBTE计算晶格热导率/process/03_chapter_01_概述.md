# 第一章：概述

本教程以金刚石结构 Si 为例，演示从 ABACUS 第一性原理计算出发，
结合 Phonopy、thirdorder 和 ShengBTE，完整地计算晶格热导率。
计算采用数值原子轨道基组（LCAO），并在结尾给出平面波（PW）基组的参数差异。

整体流程如下：

```
  ABACUS (relax)
        │ 优化后 STRU
        ▼
  Phonopy -d                        thirdorder sow
        │ 超胞构型 ×1                      │ 微扰构型 ×40
        ▼                                 ▼
  ABACUS SCF                     ABACUS SCF ×40
        │ 原子受力                         │ 原子受力
        ▼                                 ▼
  phonopy -f                       aba2vasp.py
        │ FORCE_SETS                       │ vasprun.xml ×40
        ▼                                 ▼
  phonopy band.conf                thirdorder reap
        │ FORCE_CONSTANTS                  │
      au2si.py                            │
        │                                 │
        ▼                                 ▼
  FORCE_CONSTANTS_2ND         FORCE_CONSTANTS_3RD
                    │                 │
                    └────────┬────────┘
                             ▼
                         ShengBTE
                             │
                             ▼
                       κ(T)  [W/(m·K)]
```

体系信息：
- 材料：Si，金刚石结构（Fd-3m），2 原子原胞
- 基组：LCAO（主），PW（对比）
- 赝势：`Si_ONCV_PBE-1.0.upf`（模守恒，GGA-PBE）
- 轨道文件：`Si_gga_7au_100Ry_2s2p1d.orb`（DZP，7 au 截断，100 Ry 能量截断）
