# testCLAUDE.md 模块化拆分报告

**日期：** 2026-03-09
**涉及文件：** `testCLAUDE.md`、`CLAUDE.md`、`tools/testCLAUDE/`（新建目录）

---

## 背景

`testCLAUDE.md` 原有 1934 行，单次读取约占 40K tokens。结合 CLAUDE.md（~15K）和对话历史（~30–60K），实际测试中已出现自动 compact，导致关键状态（job ID、路径变量）丢失。

---

## 变更内容

### 1. 新建模块目录 `tools/testCLAUDE/`

将主文件中低频使用的大块内容拆分为 5 个独立模块，按需读取：

| 模块文件 | 来源 | 行数 | 触发条件 |
|---------|------|------|---------|
| `bohrium_setup.md` | Step 0.2 + 附录A | 180 | Step 0 环境未配置时 |
| `plugin_dev_guide.md` | Step 1.5 | 157 | `detected_types = []` 时 |
| `plugins_history.md` | 附录D | 59 | 创建新插件时 |
| `troubleshooting.md` | 附录B | 136 | Step 4.4 失败 / Step 6.3 偏差时 |
| `file_formats.md` | 附录C | 109 | Step 2 生成输入文件时 |

### 2. 修改 `testCLAUDE.md`

- Step 0.2（113 行）→ 3 行桩代码，指向 `bohrium_setup.md`
- Step 1.5（148 行）→ 4 行桩代码，指向 `plugin_dev_guide.md` + `plugins_history.md`
- 四个附录（344 行）→ 12 行参考表，列出所有模块的触发条件
- 内部 2 处"附录B"引用 → `tools/testCLAUDE/troubleshooting.md`

### 3. 修改 `CLAUDE.md` Step 8

在用户同意执行测试后、读取 `testCLAUDE.md` 之前，增加暂停点：

- 明确输出 `tutorial_path` 和 `test_dir` 的值（确保 compact 后摘要保留关键信息）
- 建议用户运行 `/compact` 释放上下文空间
- 收到"继续测试"后再进入测试流程

---

## 效果

| 指标 | 改造前 | 改造后 |
|------|--------|--------|
| testCLAUDE.md 行数 | 1934 | **1346**（↓588 行，30%）|
| 正常流程读入 tokens | ~40K | ~26K |
| 新增插件需修改的文件 | testCLAUDE.md 主文件 | 仅 `plugins_history.md` |

---

## 相关文件

- 设计文档：`docs/plans/2026-03-09-testclaude-modularization-design.md`
- 实施计划：`docs/plans/2026-03-09-testclaude-modularization-plan.md`
