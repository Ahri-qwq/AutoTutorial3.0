#!/usr/bin/env python3
"""组装教程2：移除结构优化章节，只保留弹性常数计算"""

# 读取原始教程2
with open("教程2_原始.md", "r", encoding="utf-8") as f:
    lines = f.readlines()

# 找到关键位置
chapter2_start = 0  # h-BN 变胞优化章节开始
chapter3_start = 0  # Si 弹性常数章节开始

for i, line in enumerate(lines):
    if '# 二、案例一：h-BN 变胞优化' in line or '## 2.' in line:
        chapter2_start = i
    elif '# 三、案例二：Si 弹性常数' in line or '## 3.' in line:
        chapter3_start = i

# 组装新教程
new_tutorial = []

# 1. 修改元数据
new_tutorial.append('---\n')
new_tutorial.append('title: "ABACUS 弹性常数计算"\n')
new_tutorial.append('author: "AutoTutorial 3.0"\n')
new_tutorial.append('date: "2026-03-26"\n')
new_tutorial.append('topic: "弹性常数"\n')
new_tutorial.append('task_type: "D"\n')
new_tutorial.append('has_case: true\n')
new_tutorial.append('word_count: ~2000\n')
new_tutorial.append('chapters: 2\n')
new_tutorial.append('---\n\n')
new_tutorial.append('# ABACUS 弹性常数计算\n\n')
new_tutorial.append('> 本文由 AutoTutorial 3.0 提供。\n\n')
new_tutorial.append('---\n\n')

# 2. 修改前言
new_tutorial.append('## 前言\n\n')
new_tutorial.append('本教程介绍如何使用 ABACUS 和 abacustest 计算材料的弹性常数。\n\n')
new_tutorial.append('**前置知识：** 已掌握 ABACUS 基本使用和结构优化方法（参见教程1）\n\n')
new_tutorial.append('**教程结构：**\n')
new_tutorial.append('- 第一章：弹性张量理论背景\n')
new_tutorial.append('- 第二章：Si 弹性常数计算（abacustest 完整流程）\n\n')
new_tutorial.append('---\n\n')

# 3. 保留弹性张量理论部分（从1.2开始）
new_tutorial.append('## 第一章 弹性张量理论\n\n')
for i in range(len(lines)):
    if '## 1.2 弹性张量与广义胡克定律' in lines[i]:
        # 从1.2开始到第二章之前
        for j in range(i, chapter2_start):
            new_tutorial.append(lines[j])
        break

new_tutorial.append('\n---\n\n')

# 4. 保留Si弹性常数章节（第三章改为第二章）
new_tutorial.append('## 第二章 案例：Si 弹性常数计算\n\n')
for i in range(chapter3_start + 1, len(lines)):
    line = lines[i]
    # 跳过原第四章
    if '# 四、' in line or '## 4.' in line:
        break
    # 调整章节编号
    line = line.replace('## 3.', '## 2.')
    new_tutorial.append(line)

# 写入文件
with open("教程2_弹性常数计算.md", "w", encoding="utf-8") as f:
    f.writelines(new_tutorial)

print("Tutorial 2 assembled successfully!")
print(f"Total lines: {len(new_tutorial)}")
