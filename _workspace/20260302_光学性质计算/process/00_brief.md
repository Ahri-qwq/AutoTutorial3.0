# 任务简报

## 任务信息
- **主题：** 光学性质（介电函数/吸收谱）计算
- **任务类型：** C（案例驱动教程）
- **案例文件：** data/input/介电函数与线性光学的性质的计算.md
- **特殊要求：** 围绕案例展开
- **日期：** 2026-03-02

## 案例核心信息

### 材料
晶态 SiO₂（单胞，立方结构，8 Si + 16 O）

### 计算流程
1. **ABACUS SCF**：LCAO基组自洽计算，输出 HR/SR/rR 稀疏矩阵
2. **PYATB**：读取矩阵文件，计算光学电导率/介电函数
3. **Python 后处理**：由介电函数计算折射率、消光系数、吸收系数、能量损失函数

### ABACUS 关键参数
- basis_type = lcao
- ecutwfc = 100 Ry
- smearing_method = gaussian, smearing_sigma = 0.01
- out_mat_hs2 = 1（输出 HR/SR 矩阵）
- out_mat_r = 1（输出偶极矩阵）
- 赝势：Si/O ONCV PBE；轨道：2s2p1d-7au
- KPT：6×6×6 Gamma
- 运行：mpirun -np 16 abacus

### PYATB 关键参数
- fermi_energy = 5.5385382545 eV, occ_band = 64
- omega 0~30 eV, domega = 0.01, eta = 0.1
- grid = 20×20×20
- 运行：mpirun -np 16 pyatb

### 计算环境
- 镜像：registry.dp.tech/dptech/prod-19853/abacus-pyatb-open:v0.0.1
- 推荐配置：c16_m32_cpu

## 执行计划
Step 1 → RAG检索 + 风格学习
Step 2 → 提供3个大纲方案，等待用户选择
Step 3-7 → 撰写、审查、输出
