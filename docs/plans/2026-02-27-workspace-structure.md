# Workspace Structure Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 在 CLAUDE.md 中将所有过程文件路径加上 `process/` 前缀，使 `_workspace/YYYYMMDD_主题/` 文件夹打开后只显示最终产物 `07_Final_Tutorial_*.md` 和 `08_summary.md`。

**Architecture:** 纯文本替换，修改 `CLAUDE.md` 中约 20 处文件路径引用。Step 0 新增创建 `process/` 子目录的指令；Steps 1-7.1b 所有输出文件加 `process/` 前缀；Step 7.3、7.4 最终文件不变。无代码逻辑改动。

**Tech Stack:** 直接 Edit CLAUDE.md，grep 验证，单次 commit。

---

## Task 1: 更新 Step 0（目录创建 + 00_brief.md）

**Files:**
- Modify: `CLAUDE.md`（Step 0 区域，约第 62-71 行）

**Step 1: 修改 Step 0 的执行步骤**

使用 Edit 工具，将：

```
1. 询问任务信息（主题、是否有案例、是否案例驱动、特殊要求）
2. 判断任务类型
3. 创建工作目录：`_workspace/YYYYMMDD_HHMMSS_主题/`
4. 保存任务简报到 `00_brief.md`
```

替换为：

```
1. 询问任务信息（主题、是否有案例、是否案例驱动、特殊要求）
2. 判断任务类型
3. 创建工作目录：`_workspace/YYYYMMDD_HHMMSS_主题/`，同时在其中创建 `process/` 子目录
4. 保存任务简报到 `process/00_brief.md`
```

**Step 2: 验证修改**

```bash
grep -n "process/\|00_brief" CLAUDE.md | head -5
```

预期输出包含：
```
69:3. 创建工作目录：`_workspace/YYYYMMDD_HHMMSS_主题/`，同时在其中创建 `process/` 子目录
70:4. 保存任务简报到 `process/00_brief.md`
```

---

## Task 2: 更新 Step 1（01_research.md，4处）

**Files:**
- Modify: `CLAUDE.md`（Step 1 区域，约第 75-108 行）

**Step 1: 修改「保存到 \`01_research.md\`」**

将：
```
- 保存到 `01_research.md`
```
替换为：
```
- 保存到 `process/01_research.md`
```

**Step 2: 修改两处「追加到 \`01_research.md\`」**

将：
```
- 追加到 `01_research.md`
```
替换为：
```
- 追加到 `process/01_research.md`
```

注意：有两处，需使用 `replace_all: true` 参数。

**Step 3: 修改完成标志**

将：
```
**完成标志：** `01_research.md` 包含检索结果、案例解析（如有）、风格总结
```
替换为：
```
**完成标志：** `process/01_research.md` 包含检索结果、案例解析（如有）、风格总结
```

**Step 4: 验证**

```bash
grep -n "01_research" CLAUDE.md
```

预期：所有出现均含 `process/01_research.md`，无裸露的 `01_research.md`。

---

## Task 3: 更新 Step 2（01_research.md 引用 + 02_outline.md，3处）

**Files:**
- Modify: `CLAUDE.md`（Step 2 区域，约第 111-145 行）

**Step 1: 修改「回顾 \`01_research.md\`」**

将：
```
- 回顾 `01_research.md`
```
替换为：
```
- 回顾 `process/01_research.md`
```

**Step 2: 修改「保存3个方案到 \`02_outline.md\`」**

将：
```
- 保存3个方案到 `02_outline.md`
```
替换为：
```
- 保存3个方案到 `process/02_outline.md`
```

**Step 3: 修改完成标志**

将：
```
**完成标志：** `02_outline.md` 包含3个方案，用户已选择并确认
```
替换为：
```
**完成标志：** `process/02_outline.md` 包含3个方案，用户已选择并确认
```

**Step 4: 验证**

```bash
grep -n "02_outline\|回顾.*01_research" CLAUDE.md
```

预期：所有出现均含 `process/` 前缀。

---

## Task 4: 更新 Step 3（03_*.md，6处）

**Files:**
- Modify: `CLAUDE.md`（Step 3 区域，约第 148-198 行）

**Step 1: 修改前言保存路径**

将：
```
- 保存为 `03_00_前言.md`
```
替换为：
```
- 保存为 `process/03_00_前言.md`
```

**Step 2: 修改章节文件路径说明**

将：
```
- 保存为独立章节文件：`03_chapter_XX_标题.md`
```
替换为：
```
- 保存为独立章节文件：`process/03_chapter_XX_标题.md`
```

**Step 3: 修改文件命名规范示例（4个示例路径）**

将：
```
  - `03_00_前言.md`
  - `03_chapter_01_理论基础.md`
  - `03_chapter_02_参数设置.md`
  - `03_chapter_03_案例实战.md`
```
替换为：
```
  - `process/03_00_前言.md`
  - `process/03_chapter_01_理论基础.md`
  - `process/03_chapter_02_参数设置.md`
  - `process/03_chapter_03_案例实战.md`
```

**Step 4: 修改附录保存路径**

将：
```
- 保存为 `03_99_附录.md`
```
替换为：
```
- 保存为 `process/03_99_附录.md`
```

**Step 5: 修改初稿汇总路径（含完成标志）**

将：
```
- 保存为 `03_draft_full.md`（完整初稿）
```
替换为：
```
- 保存为 `process/03_draft_full.md`（完整初稿）
```

将：
```
**完成标志：** `03_draft_full.md` 包含完整初稿（前言+正文+附录）
```
替换为：
```
**完成标志：** `process/03_draft_full.md` 包含完整初稿（前言+正文+附录）
```

**Step 6: 验证**

```bash
grep -n "03_" CLAUDE.md | grep -v "process/"
```

预期：无输出（即所有 `03_` 开头的文件路径都已加 `process/` 前缀）。

---

## Task 5: 更新 Steps 4-6（04、05、06 审查文件，6处）

**Files:**
- Modify: `CLAUDE.md`（Steps 4-6 区域，约第 201-286 行）

**Step 1: 修改 04_review_content.md（2处）**

将：
```
- 创建 `04_review_content.md`
```
替换为：
```
- 创建 `process/04_review_content.md`
```

将：
```
**完成标志：** `04_review_content.md` 包含审查报告和修改后的稿件
```
替换为：
```
**完成标志：** `process/04_review_content.md` 包含审查报告和修改后的稿件
```

**Step 2: 修改 05_review_case.md（2处）**

将：
```
- 创建 `05_review_case.md`
```
替换为：
```
- 创建 `process/05_review_case.md`
```

将：
```
**完成标志：** `05_review_case.md` 包含审查报告和修改后的稿件
```
替换为：
```
**完成标志：** `process/05_review_case.md` 包含审查报告和修改后的稿件
```

**Step 3: 修改 06_review_style.md（2处）**

将：
```
- 创建 `06_review_style.md`
```
替换为：
```
- 创建 `process/06_review_style.md`
```

将：
```
**完成标志：** `06_review_style.md` 包含审查报告和修改后的稿件
```
替换为：
```
**完成标志：** `process/06_review_style.md` 包含审查报告和修改后的稿件
```

**Step 4: 验证**

```bash
grep -n "04_review\|05_review\|06_review" CLAUDE.md | grep -v "process/"
```

预期：无输出。

---

## Task 6: 更新 Step 7.1 和 7.1b（07_fix.md，3处）

**Files:**
- Modify: `CLAUDE.md`（Step 7 区域，约第 289-336 行）

**Step 1: 修改 7.1 整合步骤**

将：
```
- 整合 Step 4、5、6 的所有修改，生成`07_fix.md`
```
替换为：
```
- 整合 Step 4、5、6 的所有修改，生成`process/07_fix.md`
```

**Step 2: 修改 7.1b 的 orbital_validator 命令**

将：
```
python tools/orbital_validator.py 07_fix.md --fix
```
替换为：
```
python tools/orbital_validator.py process/07_fix.md --fix
```

**Step 3: 修改 7.1b 的「更新 07_fix.md」**

将：
```
- 有 NG 问题并自动修正时：确认修正内容合理，更新 07_fix.md
```
替换为：
```
- 有 NG 问题并自动修正时：确认修正内容合理，更新 process/07_fix.md
```

**Step 4: 验证 07_fix.md 全部加了前缀，且 07_Final_Tutorial 和 08_summary 未加**

```bash
grep -n "07_fix\|07_Final\|08_summary" CLAUDE.md
```

预期：
- `07_fix` 的所有引用均含 `process/`
- `07_Final_Tutorial` 和 `08_summary` 均**不含** `process/`

---

## Task 7: 最终验证与提交

**Step 1: 整体扫描，确认无遗漏的裸路径**

```bash
grep -n "`00_brief\|`01_research\|`02_outline\|`03_\|`04_review\|`05_review\|`06_review\|`07_fix" CLAUDE.md | grep -v "process/"
```

预期：无输出（所有过程文件路径均含 `process/` 前缀）。

**Step 2: 确认最终文件路径未被修改**

```bash
grep -n "07_Final_Tutorial\|08_summary" CLAUDE.md
```

预期：不含 `process/`，路径正确。

**Step 3: 提交**

```bash
git add CLAUDE.md
git commit -m "feat: store all process files in process/ subfolder"
```
