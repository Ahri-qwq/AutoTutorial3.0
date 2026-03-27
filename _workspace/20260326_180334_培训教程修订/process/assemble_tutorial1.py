#!/usr/bin/env python3
"""组装教程1：插入新章节到原始教程"""

# 读取原始教程
with open("教程1_原始.md", "r", encoding="utf-8") as f:
    original = f.read()

# 读取新章节
with open("process/01_教程1_软件背景章节.md", "r", encoding="utf-8") as f:
    software_bg = f.read()

with open("process/05_教程1_结构优化章节.md", "r", encoding="utf-8") as f:
    relax_chapter = f.read()

with open("process/03_教程1_abacustest准备输入.md", "r", encoding="utf-8") as f:
    abacustest_prepare = f.read()

with open("process/04_教程1_abacustest抓取结果.md", "r", encoding="utf-8") as f:
    abacustest_extract = f.read()

with open("process/02_教程1_输出文件解读章节.md", "r", encoding="utf-8") as f:
    output_files = f.read()

# 分割原始教程
lines = original.split('\n')

# 找到关键位置
header_end = 0
chapter1_start = 0
appendix_start = 0

for i, line in enumerate(lines):
    if line.startswith('## 前言'):
        header_end = i
    elif line.startswith('## 第一章'):
        chapter1_start = i
    elif line.startswith('## 附录'):
        appendix_start = i

# 组装新教程
new_tutorial = []

# 1. 元数据和标题（修改）
new_tutorial.append('---')
new_tutorial.append('title: "ABACUS 基本介绍：输入输出体系与结构优化"')
new_tutorial.append('author: "AutoTutorial 3.0"')
new_tutorial.append('date: "2026-03-26"')
new_tutorial.append('topic: "ABACUS 基本介绍"')
new_tutorial.append('task_type: "D"')
new_tutorial.append('has_case: true')
new_tutorial.append('word_count: ~6500')
new_tutorial.append('chapters: 7')
new_tutorial.append('---')
new_tutorial.append('')
new_tutorial.append('# ABACUS 基本介绍：输入输出体系与结构优化')
new_tutorial.append('')
new_tutorial.append('> 本文由 AutoTutorial 3.0 提供。')
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 2. 插入软件背景
new_tutorial.append('## 第零章 ABACUS 软件背景')
new_tutorial.append('')
new_tutorial.append(software_bg)
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 3. 前言（修改教程结构表）
new_tutorial.append('## 前言')
new_tutorial.append('')
new_tutorial.append('本教程介绍 ABACUS 的基本使用，包括软件背景、三个核心输入文件（STRU、KPT、INPUT）、结构优化、输出文件解读，以及 abacustest 工具的使用。')
new_tutorial.append('')
new_tutorial.append('**教程目标**')
new_tutorial.append('')
new_tutorial.append('- 了解 ABACUS 软件的背景和特点')
new_tutorial.append('- 理解 STRU、KPT、INPUT 各自承担什么职责')
new_tutorial.append('- 学会为平面波（PW）和数值原子轨道（LCAO）两种基组写 STRU')
new_tutorial.append('- 掌握 k 点采样的两种格式：MP 网格（SCF 用）和路径 k 点（能带用）')
new_tutorial.append('- 理解 INPUT 中的关键参数分组，重点掌握 mixing 参数的物理含义')
new_tutorial.append('- 掌握结构优化（relax 和 cell-relax）的原理和参数设置')
new_tutorial.append('- 学会解读 ABACUS 输出文件（running_scf.log）')
new_tutorial.append('- 掌握 abacustest 工具准备输入和抓取结果的方法')
new_tutorial.append('')
new_tutorial.append('**适用读者**')
new_tutorial.append('')
new_tutorial.append('初次使用 ABACUS 的用户，或熟悉其他 DFT 软件（VASP、QE）需要迁移的用户。建议具备基本的 DFT 概念（自洽场计算、布里渊区采样）。')
new_tutorial.append('')
new_tutorial.append('**前置知识**')
new_tutorial.append('')
new_tutorial.append('- 晶体结构基本概念（晶格常数、分数坐标）')
new_tutorial.append('- 自洽场（SCF）计算的基本流程')
new_tutorial.append('')
new_tutorial.append('**教程结构**')
new_tutorial.append('')
new_tutorial.append('| 章节 | 主题 |')
new_tutorial.append('|------|------|')
new_tutorial.append('| 第零章 | ABACUS 软件背景 |')
new_tutorial.append('| 第一章 | STRU 结构文件 |')
new_tutorial.append('| 第二章 | KPT k 点文件 |')
new_tutorial.append('| 第三章 | INPUT 计算参数 |')
new_tutorial.append('| 第四章 | 结构优化 |')
new_tutorial.append('| 第五章 | 输出文件解读 |')
new_tutorial.append('| 第六章 | abacustest 准备输入 |')
new_tutorial.append('| 第七章 | abacustest 抓取结果 |')
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 4. 第一、二、三章（原始内容）
new_tutorial.extend(lines[chapter1_start:appendix_start])
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 5. 第四章：结构优化
new_tutorial.append('## 第四章 结构优化')
new_tutorial.append('')
new_tutorial.append(relax_chapter)
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 6. 第五章：输出文件解读（替换原附录）
new_tutorial.append('## 第五章 输出文件解读')
new_tutorial.append('')
new_tutorial.append(output_files)
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 7. 第六章：abacustest 准备输入
new_tutorial.append('## 第六章 使用 abacustest 准备输入文件')
new_tutorial.append('')
new_tutorial.append(abacustest_prepare)
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 8. 第七章：abacustest 抓取结果
new_tutorial.append('## 第七章 使用 abacustest 抓取计算结果')
new_tutorial.append('')
new_tutorial.append(abacustest_extract)
new_tutorial.append('')
new_tutorial.append('---')
new_tutorial.append('')

# 9. 附录（简化版）
new_tutorial.append('## 附录：进阶学习方向')
new_tutorial.append('')
new_tutorial.append('本教程覆盖了 ABACUS 的基本使用。后续可以进一步学习：')
new_tutorial.append('')
new_tutorial.append('- **弹性常数计算**：在结构优化基础上，计算材料的力学性质')
new_tutorial.append('- **能带与 DOS**：计算电子结构，分析材料的电学性质')
new_tutorial.append('- **LCAO 基组与收敛**：LCAO 计算中基组大小对精度的影响')
new_tutorial.append('- **高级功能**：DFT+U、隐式溶剂模型、RT-TDDFT 等')
new_tutorial.append('')
new_tutorial.append('### 参考资料')
new_tutorial.append('')
new_tutorial.append('1. ABACUS 线上文档：https://abacus.deepmodeling.com')
new_tutorial.append('2. ABACUS GitHub：https://github.com/deepmodeling/abacus-develop')
new_tutorial.append('3. abacustest GitHub：https://github.com/pxlxingliang/abacus-test')

# 写入文件
output = '\n'.join(new_tutorial)
with open("教程1_ABACUS基本介绍.md", "w", encoding="utf-8") as f:
    f.write(output)

print("Tutorial 1 assembled successfully!")
print(f"Total lines: {len(new_tutorial)}")
print(f"Output file: tutorial1_ABACUS_basics.md")
