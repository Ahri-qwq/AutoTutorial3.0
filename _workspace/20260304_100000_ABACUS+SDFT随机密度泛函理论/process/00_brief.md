# 任务简报

## 任务信息
- **主题：** 随机密度泛函理论（SDFT）在 ABACUS 中的使用
- **案例文件：** data/input/abacus_sdft.md
- **任务类型：** C（案例驱动教程）
- **特殊要求：** 标题需包含"ABACUS"
- **日期：** 2026-03-04

## 案例概述
案例文件包含三个算例：
1. **pw_Si2** — 2原子Si，T=0.6 Ry（≈8.16 eV），SCF 计算
2. **pw_md_Al** — 16原子Al，T=7.35 Ry（≈100 eV），MD 模拟
3. **186_PW_SDOS_10D10S** — 1原子Si，T=0.6 Ry，DOS 计算（MDFT，10个KS轨道+10个随机轨道）

## 核心概念
- SDFT：随机波函数密度泛函理论，用随机轨道替代对角化，适用于高温高压计算
- MDFT：混合KS轨道与随机轨道，加速收敛
- 目标场景：温稠密物质（Warm Dense Matter，WDM），温度数十至上千eV

## 教程目标
- 介绍 SDFT/MDFT 的物理原理和适用场景
- 详细说明关键参数（nbands, nbands_sto, nche_sto, method_sto等）
- 通过三个案例展示 SCF、MD、DOS 的完整使用流程
