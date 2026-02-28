# AutoTutorial 3.0 - 测试插件注册表

本文件记录所有已注册的测试插件。**每次新增插件后必须更新此表。**

当 `analysis.json` 为空（`detected_types = []`）时，Claude 应首先查阅本表，判断是已有插件未匹配还是需要创建新插件。

---

## 已注册插件

| 插件文件 | 计算类型 | can_handle 检测关键词 | 添加日期 | 首次测试教程 |
|----------|----------|----------------------|----------|--------------|
| relax_plugin.py | relax | `calculation\s*=\s*(relax\|cell-relax)` | 2026-02 | 弹性常数教程（relax步骤） |
| elastic_plugin.py | elastic | `弹性常数`、`elastic` | 2026-02 | 弹性常数教程 |
| band_plugin.py | band | `能带结构`、`band structure` | 2026-02 | 能带结构教程 |
| dos_plugin.py | dos | `态密度`、`density of states` | 2026-02 | 态密度教程 |
| dftu_plugin.py | dftu | `dft_plus_u\s*[=\s]+1`、`DFT\+U` | 2026-02-28 | DFT+U强关联体系教程（NiO） |

---

## 诊断流程

当 `analysis.json` 为空时：

```
1. 阅读教程，找到核心 INPUT 关键词
2. 查表：已注册插件中是否有匹配的 can_handle 关键词？
   ├── 有匹配 → 关键词格式问题（空格/大小写），检查教程写法
   └── 无匹配 → 需要创建新插件，执行 testCLAUDE.md Step 1.5
```

---

## 新增插件时的更新步骤

1. 创建 `tools/test_plugins/<type>_plugin.py`
2. 在 `tools/test_framework_integrated.py` 中导入并注册
3. 在本表追加一行（格式同上）
4. 在 `testCLAUDE.md` 正文"已支持计算类型"表中追加一行
5. 在 `testCLAUDE.md` 附录D"插件开发历史"表中追加一行

---

## 已知待扩展的计算类型

以下计算类型尚无插件，遇到时按 testCLAUDE.md Step 1.5 创建：

| 计算类型 | 核心关键词 | 主要挑战 |
|----------|-----------|----------|
| HSE06 杂化泛函 | `exx_hybrid_type = hse` | 计算量大，需增大 maxRunTime |
| 声子谱 Phonopy | Phonopy 后处理脚本 | 多步骤：SCF → 位移 → 后处理 |
| 实时 TDDFT | `esolver_type = tddft` | 结果是轨迹，难以定量对比 |
| 分子动力学 MD | `calculation = md` | 轨迹验证，不用对比单一数值 |
| Yukawa DFT+U | `yukawa_potential = 1` | 参数格式与标准 DFT+U 略有不同 |

---

**文档版本：** 1.0
**创建日期：** 2026-02-28
**维护者：** AutoTutorial 3.0
