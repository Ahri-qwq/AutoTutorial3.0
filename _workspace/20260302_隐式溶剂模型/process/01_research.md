# Step 1: 知识调研与案例解析

## 1.1 RAG检索结果

### 查询1：隐式溶剂模型原理
- **主要来源：** abacus_user_guide__abacus_sol.html.md（与案例文件内容一致）
- **评估：** 知识库中最主要的资料就是案例文件本身，无额外扩展内容
- **有效信息：**
  - Mathew et al. J. Chem. Phys. 140, 084106 (2014) 是 ABACUS 实现的理论依据
  - 溶剂作为连续介质处理（非显式溶剂分子）
  - 计算成本远低于显式溶剂方法

### 查询2：implicit solvation 物理原理
- **主要来源：** en__latest__advanced__scf__advanced.html.md（ABACUS 官方文档）
- **有效信息：**
  - 官方文档位置确认：advanced/scf/advanced.html#implicit-solvation-model
  - Pt-slab 是官方标准算例

### 查询3：ABACUS SCF PW基础
- **有效信息：** 平面波基组 SCF 计算基础步骤，参数如 ecutwfc, scf_thr, basis_type 用法

## 1.2 案例解析

### 案例来源
- 文件：`data/input/ABACUS隐式溶剂模型使用教程.md`
- 原始来源：https://mcresearch.github.io/abacus-user-guide/abacus-sol.html
- 作者：刘裕、孙梦琳；审核：许审镇、陈默涵

### 案例体系
- **物理体系：** H₂分子置于超胞中，模拟水溶液环境
- **计算方法：** PW基组，SCF

### 案例文件结构（完整性检查清单）
- [x] INPUT — 完整，见案例原文
- [ ] STRU — **案例中未列出**（需后续注明为算例目录中的文件）
- [ ] KPT — **案例中未列出**（同上）
- [ ] 赝势文件 — 未明确列出（H 的赝势）

### INPUT 文件（完整内容）
```
INPUT_PARAMETERS
#Parameters (1.General)
suffix              H2
calculation         scf
ntype               1
nbands              2
symmetry            0
pseudo_dir          ./

#Parameters (2.Iteration)
ecutwfc             60
scf_thr             1e-6
scf_nmax            100

#Parameters (3.Basis)
basis_type          pw

#Parameters (Solvation Model)
imp_sol             1
eb_k                80
tau                 0.000010798
sigma_k             0.6
nc_k                0.00037
```

### 关键参数说明（来自案例）
| 参数 | 类型 | 值 | 说明 |
|------|------|-----|------|
| `imp_sol` | Bool | 1 | 隐式溶剂模型开关，1=开 |
| `eb_k` | Real | 80 | 溶剂相对介电常数（水=80） |
| `tau` | Real | 0.000010798 | 有效表面张力参数，单位 Ry/Bohr² |
| `sigma_k` | Real | 0.6 | 扩散腔宽度（无量纲） |
| `nc_k` | Real | 0.00037 | 介电腔形成时的电子密度，单位 Bohr⁻³ |

### 预期输出
- 输出文件：`OUT.ABACUS/running_scf.log`
- 关键量：`E_sol_el`（静电作用，通常为负值）
- 关键量：`E_sol_cav`（空腔作用，通常为正值）
- 溶剂化能计算：E(imp_sol=1) - E(imp_sol=0)（结构优化后）

### 案例完整性评估
- ✅ INPUT 文件完整
- ⚠️ STRU/KPT 未在案例文件中列出，需引导读者下载算例
- ✅ 预期结果有定性描述

## 1.3 风格参考总结

### 参考文章特点
1. **磁性材料计算教程：**
   - 直接定义概念（不绕弯子），使用表格对比参数
   - 代码块带语言标签，分节清晰
   - 段落短，信息密度高

2. **ABACUS+DeePMD-kit 教程：**
   - "本教程旨在介绍..." 直接说明目的
   - 有"准备"章节说明算例下载方式
   - 分节用二级+三级标题

### 写作风格关键点
- 开篇直接说明目的和适用场景
- 有专门的"算例获取"说明
- 参数用表格+代码双呈现
- 预期结果有具体输出描述
- 不用"在当今"、"综上所述"等 AI 腔表达
