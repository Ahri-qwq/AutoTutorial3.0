# Step 1 研究结果

## 1.1 RAG 检索结果评估

### 主要查询：abacustest model inputs 工具
检索质量：**高度相关**。知识库中有多条关于 abacustest model inputs 的文档，涵盖了：
- 完整命令行参数列表（-f, --ftype, --jtype, --pp, --orb, --input, --kpt, --lcao 等）
- 安装方式（pip 和 git+pip）
- 使用 ABACUS_PP_PATH / ABACUS_ORB_PATH 环境变量
- 批量准备输入文件（param.json 方式和 model inputs 方式）

注意：知识库中还存在旧版方式（ASE-ABACUS 接口），本教程案例介绍的是更新的 abacustest model inputs 方式，无需混淆。

### 补充查询：磁性/DFT+U
检索结果与案例高度匹配：知识库有 nspin=2、初始磁矩设置的完整文档，可作为知识背景。

### 补充查询：批量/自定义 INPUT
检索到 kspacing 参数说明和批量使用案例，与案例内容一致。

---

## 1.2 案例解析

### 案例核心：5 个递进场景

| 场景 | 体系 | 关键特性 |
|------|------|----------|
| 基础 SCF | MgO（CIF） | LCAO，自动赝势轨道，生成完整目录 |
| 磁性+DFT+U | Fe₂O₃（CIF） | nspin=2, init_mag, dftu, dftu_param |
| 多任务类型 | MgO（cell-relax） | jtype, relax 参数自动设置 |
| 自定义 INPUT/KPT | Fe₂O₃ | --input 模板，--kpt 5 5 5 |
| 批量结构 | Pd(100) 系列 | --folder-syntax Python 切片 |

### 关键命令汇总
```bash
# 下载赝势轨道库
abacustest model inputs --download-pporb

# MgO SCF（基础）
abacustest model inputs -f MgO.cif --ftype cif --lcao --folder-syntax MgO

# Fe₂O₃ 磁性+DFT+U
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0

# MgO cell-relax
abacustest model inputs -f MgO.cif --ftype cif --lcao --jtype cell-relax --folder-syntax MgO-cellrelax

# 自定义 INPUT 和 KPT
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --input INPUT_template --kpt 5 5 5 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0

# 批量处理
abacustest model inputs -f Pd100_*layer.vasp --ftype poscar --lcao --jtype relax --folder-syntax "x[:-5]"
```

### 案例完整性检查清单
- [x] MgO 目录结构（含符号链接说明）
- [x] MgO INPUT 文件内容
- [x] Fe₂O₃ INPUT 文件内容（含 DFT+U 参数）
- [x] Fe₂O₃ STRU 文件内容（含 12 个 Fe + 18 个 O 原子坐标）
- [x] MgO cell-relax INPUT 文件内容
- [x] Co₂FeAl 多元素磁矩+DFT+U 示例命令
- [x] Pd(100) 批量结构名称列表 + folder-syntax 命令

---

## 1.3 风格参考总结

**磁性材料教程**
- 开头：直接介绍概念，不用虚词
- 使用表格比较参数
- 分小节细化，层级清晰（###/####）
- 语言直接、技术性强

**ABACUS+DeePMD-kit 教程**
- 编号章节（一、二、三）
- 每步给出完整命令
- 较多代码块，少理论展开
- 适合工具使用类教程

**PySCF 教程**
- 学术风格，理论展开较多
- 不适合本篇（工具使用类）

**本教程适用风格：** 接近 ABACUS+DeePMD-kit 教程风格——步骤驱动，代码为主，适量说明参数含义，不过度展开理论。

**参考文章平均长度：** ~345 行、~883 words
**本教程目标长度：** 500-800 行（工具教程，以案例代码为主，适中）
