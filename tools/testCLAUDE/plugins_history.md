# 插件开发历史与扩展参考

> 本文件由 testCLAUDE.md 按需 read。
> **新增插件时，只需在本文件"插件开发历史"表末尾追加一行，无需修改主文件。**

---

## 附录D: 未来计算类型扩展参考

### 已知待扩展场景

| 教程主题 | 核心关键词 | 主要挑战 | 参考 |
|----------|-----------|----------|------|
| HSE06 杂化泛函 | `exx_hybrid_type = hse` | 计算量大，需增大 maxRunTime | — |
| 声子谱 Phonopy（已完成）| `phonopy`、`FORCE_SETS` | 多步骤：SCF → 位移 → 后处理 | phonopy_plugin.py |
| 实时 TDDFT | `esolver_type = tddft` | 结果是轨迹，难以定量对比 | — |
| 分子动力学 MD | `calculation = md` | 轨迹验证，不用对比单一数值 | — |
| DFT+U（已完成）| `dft_plus_u = 1` | 结构不完整，磁矩容差 | dftu_plugin.py |

### 插件开发历史

| 插件 | 计算类型 | 添加日期 | 首次测试教程 |
|------|----------|----------|--------------|
| relax_plugin.py | relax | 2026-02 | 弹性常数教程 |
| elastic_plugin.py | elastic | 2026-02 | 弹性常数教程 |
| band_plugin.py | band | 2026-02 | 能带结构教程 |
| dos_plugin.py | dos | 2026-02 | 态密度教程 |
| dftu_plugin.py | dftu | 2026-02-28 | DFT+U 强关联体系教程 |
| optic_plugin.py | optic | 2026-03-02 | 光学性质（介电函数/吸收谱）计算教程（SiO₂） |
| solvation_plugin.py | solvation | 2026-03-02 | 隐式溶剂模型使用教程（H₂@水溶液） |
| phonopy_plugin.py | phonopy | 2026-03-03 | ABACUS+Phonopy 声子谱计算教程（FCC Al） |
| elf_plugin.py | elf | 2026-03-09 | ELF 电子局域函数计算与可视化教程（H₂O PW/LCAO + Fe BCC） |

> 每次新增插件时，在此表格追加一行。

---

### AutoTutorial 3.0相关
- 教程生成指南：`CLAUDE.md`
- 测试功能总结：`docs/CALCULATION_TESTING_SUMMARY.md`
- 成功报告：`docs/2SUCCESS_REPORT.md`
- 综合解决方案：`docs/1COMPREHENSIVE_SOLUTION.md`

### Bohrium平台
- 官网：https://bohrium.dp.tech/
- CLI文档：https://docs.bohrium.com/docs/bohrctl/about
- ABACUS软件案例：https://bohrium-doc.dp.tech/docs/software/ABACUS/

### ABACUS相关
- 官网：http://abacus.deepmodeling.com/
- 赝势库：http://abacus.deepmodeling.com/upf/
- 轨道文件：http://abacus.deepmodeling.com/orbital/

---

**文档版本：** 1.0
**创建日期：** 2026-02-26
**适用范围：** AutoTutorial 3.0 教程测试
**维护者：** AutoTutorial开发团队
