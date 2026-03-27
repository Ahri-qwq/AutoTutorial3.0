#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
正确重新生成教程2：只保留弹性张量理论和Si案例，移除h-BN案例
"""

# 读取原始教程2
with open('process/教程2_原始.md', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# 找到关键分界点
start_elastic_theory = None  # 1.2 弹性张量理论开始
start_si_case = None  # Si案例开始
end_file = len(lines)

for i, line in enumerate(lines):
    if '## 1.2 弹性张量与广义胡克定律' in line:
        start_elastic_theory = i
    elif '# 三、案例二：Si 弹性常数' in line:
        start_si_case = i

print(f"Found elastic theory at line {start_elastic_theory}")
print(f"Found Si case at line {start_si_case}")

# 构建新教程2
new_content = []

# 添加元数据
new_content.append('---\n')
new_content.append('title: "ABACUS 弹性常数计算"\n')
new_content.append('author: "AutoTutorial 3.0"\n')
new_content.append('date: "2026-03-26"\n')
new_content.append('topic: "弹性常数"\n')
new_content.append('task_type: "D"\n')
new_content.append('has_case: true\n')
new_content.append('word_count: ~1500\n')
new_content.append('chapters: 2\n')
new_content.append('---\n\n')

# 添加标题和前言
new_content.append('# ABACUS 弹性常数计算\n\n')
new_content.append('> 本文由 AutoTutorial 3.0 提供。\n\n')
new_content.append('## 前言\n\n')
new_content.append('本教程介绍如何使用 ABACUS 和 abacustest 工具计算材料的弹性常数。\n\n')
new_content.append('**适用读者：** 已掌握 ABACUS 基本使用和结构优化方法（参见教程1）\n\n')
new_content.append('**前置知识：** ABACUS 的 INPUT/STRU/KPT 文件格式，结构优化原理\n\n')
new_content.append('**教程结构：**\n')
new_content.append('- 第一章：弹性张量理论\n')
new_content.append('- 第二章：Si 弹性常数计算案例（完整流程）\n\n')
new_content.append('**计算环境：** LCAO 基组，ABACUS v3.10.1，abacustest 工具\n\n')
new_content.append('---\n\n')

# 添加第一章：弹性张量理论（从1.2开始到Si案例之前）
new_content.append('## 第一章 弹性张量理论\n\n')
new_content.extend(lines[start_elastic_theory:start_si_case])

# 添加第二章：Si案例（从Si案例开始到文件结束）
new_content.append('\n## 第二章 案例：Si 弹性常数计算\n\n')
# 跳过原来的"# 三、案例二：Si 弹性常数"标题
new_content.extend(lines[start_si_case+1:end_file])

# 写入新文件
with open('教程2_弹性常数计算.md', 'w', encoding='utf-8') as f:
    f.writelines(new_content)

print("Tutorial 2 regenerated successfully!")
