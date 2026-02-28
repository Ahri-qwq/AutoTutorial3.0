# 任务简报

## 基本信息
- **创建时间**：2026-02-28 14:44:16
- **任务类型**：类型C（案例驱动教程）
- **主题**：ABACUS DFT+U强关联体系计算设置
- **案例文件**：data/input/DFT+U计算.md
- **特殊要求**：围绕案例展开；完成后进行计算测试

## 案例概述
- **体系**：NiO 反铁磁体（antiferromagnetic NiO）
- **基组**：LCAO
- **计算类型**：SCF + DFT+U + occupation matrix control (omc)
- **特点**：
  - 两个Ni原子设为不同atomic species（Ni1磁矩+2、Ni2磁矩-2）
  - 对Ni的d轨道施加Hubbard U=5 eV
  - 演示omc=2固定occupation matrix的用法

## DFT+U核心参数（来自案例）
```
dft_plus_u    1
orbital_corr  2 2 -1
hubbard_u     5.0 5.0 0.0
```

## 预期结果（来自案例）
- 总能量：-9255.7279034240546025 eV
- 能隙：spin-up 0.205 eV，spin-down 2.794 eV
- 总磁矩：0（反铁磁）
- 绝对磁矩：3.35321634 Bohr mag/cell
- Ni1磁矩：+1.827，Ni2磁矩：-1.827

## 执行计划
1. Step 1：RAG知识检索 + 风格参考学习
2. Step 2：提供3个大纲方案，等待用户确认
3. Step 3：撰写完整初稿
4. Step 4-6：多轮审查
5. Step 7：最终输出
6. Step 8：计算测试
