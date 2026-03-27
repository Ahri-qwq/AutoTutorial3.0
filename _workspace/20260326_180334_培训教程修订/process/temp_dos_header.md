---
title: "用 ABACUS 计算能带与态密度：从硅的带隙到铝的费米面"
author: "AutoTutorial 3.0"
date: "2026-03-25"
topic: "电子结构"
task_type: "A"
has_case: false
word_count: ~3500
chapters: 4
---

# 用 ABACUS 计算能带与态密度：从硅的带隙到铝的费米面

> 本文由 AutoTutorial 3.0 提供。

## 前言

ABACUS 的电子结构计算包含两类核心任务：电子能带结构（band structure）和态密度（density of states，DOS）。两类计算共享同一套流程骨架：先做自洽计算（SCF）获得收敛的电子密度，再做非自洽计算（NSCF）读入该密度，在任意 k 点上求解 Kohn-Sham 方程。本教程通过两个典型案例演示完整流程：

- **案例一：硅（Si）** — 计算能带结构，验证 PBE 间接带隙约 0.57 eV
- **案例二：铝（Al）** — 计算 TDOS 和 PDOS，分析金属特征电子结构

后处理统一使用 `abacus-plot`（ABACUS 官方 Python 包）完成绘图。

**适用读者：** 已安装 ABACUS 并能运行基础 SCF 计算的用户

**前置知识：** ABACUS 的 INPUT/STRU/KPT 文件格式，DFT 基本概念

**计算环境：** 案例一使用 LCAO 基组，案例二使用平面波（PW）基组，均基于 ABACUS v3.10.1

---

