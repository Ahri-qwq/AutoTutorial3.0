# 任务简报

## 基本信息
- 日期：2026-03-04
- 任务类型：C（案例驱动教程）
- 主题：ABACUS+ShengBTE 计算晶格热导率
- 案例文件：data/input/ShengBTE计算晶格热导率.md

## 案例核心信息
- 材料体系：Si（金刚石结构，2原子原胞）
- 基组方案：LCAO（主）+ PW（对比）
- 赝势：Si_ONCV_PBE-1.0.upf
- 轨道：Si_gga_7au_100Ry_2s2p1d.orb

## 计算流程（三步）
1. 二阶力常数：ABACUS(relax → SCF) + Phonopy + ASE → FORCE_CONSTANTS_2ND
2. 三阶力常数：ABACUS(40个SCF) + thirdorder + aba2vasp.py → FORCE_CONSTANTS_3RD
3. 晶格热导率：ShengBTE → κ(T)

## 关键参数记录
- 结构优化：k=2×2×2，ecut=100 Ry
- 超胞：DIM=2 2 2
- Phonopy命令：phonopy setting.conf --abacus -d
- thirdorder命令：thirdorder_vasp.py sow 2 2 2 -2（产生40个构型）
- scf_thr：LCAO 1e-8，PW 1e-12
- ShengBTE命令：mpirun -n 10 ShengBTE
- 结果：300K下~100 W/(m·K)（实验值~150 W/(m·K)）

## 特殊要求
- 围绕案例展开（Type C）
- 标题需包含"ABACUS"

## 执行计划
Step 1 → 检索知识 + 解析案例 + 学习风格
Step 2 → 讨论大纲（等待用户选择）
Step 3-7 → 撰写、审查、输出
