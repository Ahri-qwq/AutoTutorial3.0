# 任务简报

## 基本信息
- **任务类型：** C（案例驱动教程）
- **主题：** 从结构文件准备ABACUS输入文件夹
- **工具：** abacustest model inputs
- **案例文件：** data/input/从结构文件准备ABACUS输入文件夹.md
- **特殊要求：** 围绕案例展开

## 案例核心内容
1. 下载 APNS-pp-orb-v1 赝势轨道库（`--download-pporb`）
2. MgO SCF 基础案例（CIF → 完整输入目录，LCAO 基组）
3. Fe₂O₃ 磁性 + DFT+U 案例（nspin=2，init_mag，dftu_param）
4. 不同任务类型：scf / relax / cell-relax（MgO cell-relax 案例）
5. 自定义 INPUT 和 KPT（--input / --kpt 选项）
6. 批量结构准备（Pd(100) 系列，--folder-syntax Python 切片语法）

## 执行计划
- Step 1：RAG 检索 + 案例解析 + 风格学习
- Step 2：大纲讨论（3 方案）→ 等待用户选择
- Step 3-7：撰写、审查、输出
