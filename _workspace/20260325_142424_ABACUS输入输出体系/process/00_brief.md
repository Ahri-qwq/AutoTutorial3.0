# 任务简报

- **创建时间：** 2026-03-25 14:24:24
- **主题：** ABACUS 输入输出体系
- **任务类型：** 类型A（新建教程，无外部案例文件，自行通过RAG检索获取案例数据）
- **标题方向：** ABACUS 输入输出体系详解

## 必须覆盖的子主题
1. INPUT 参数详解（分类、语法、常用参数）
2. SCF 收敛算法参数（mixing_type / mixing_beta）
3. KPT 采样策略（MP 网格 + 路径 KPT，Al 收敛测试案例）
4. STRU 结构文件定义（PW + LCAO 两种写法）
5. 平面波截断能收敛测试（ecutwfc，Al 案例）

## 知识简报摘要
- INPUT/STRU/KPT 资料均充足（RAG ✅）
- 核心案例：Si SCF（平面波，8原子）贯穿全文
- Al 收敛测试案例（k点 + ecutwfc）有完整文档支持
- 运行输出解读（running_scf.log）降级为附录
- 基组/赝势匹配、基组收敛测试本版不做
- 参考来源：
  - `ABACUS 使用教程｜电子自洽迭代.md`
  - `ABACUS 的平面波计算与收敛性测试.md`
  - `Al元素晶体...k点收敛性测试.md`
  - `en__latest__advanced__input_files__kpt.html.md`
  - `ABACUS 使用教程｜结构优化.md`

## 特殊要求
- 无外部 docx 案例文件，RAG + 联网检索获取案例数据
- 风格参考：data/reference_materials/style_references/ 中的参考文章
