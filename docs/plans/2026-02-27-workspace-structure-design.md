# Workspace Structure Design

**日期：** 2026-02-27
**状态：** 已确认，待实现

---

## 问题背景

当前 `_workspace/YYYYMMDD_主题/` 文件夹内所有文件平铺在同一层，包括过程文件（00-07_fix）和最终文件（07_Final_Tutorial_*.md、08_summary.md），打开文件夹时难以快速定位最终产物。

## 目标结构

```
_workspace/20260203_105918_弹性常数计算/
├── 07_Final_Tutorial_弹性常数计算.md   ← 打开即见
├── 08_summary.md                        ← 打开即见
└── process/
    ├── 00_brief.md
    ├── 01_research.md
    ├── 02_outline.md
    ├── 03_00_前言.md
    ├── 03_chapter_01_xxx.md
    ├── 03_draft_full.md
    ├── 04_review_content.md
    ├── 05_review_case.md
    ├── 06_review_style.md
    └── 07_fix.md
```

## 设计决策

- **子目录名：** `process/`（非 `processed/`）
- **范围：** 只影响未来新建目录，现有 `_workspace/` 下旧文件夹不做整理
- **变更类型：** 纯路径前缀变更，无逻辑改动

## 变更列表（CLAUDE.md）

| 步骤 | 文件 | 变更 |
|------|------|------|
| Step 0 | 创建目录 | 新增同时创建 `process/` 子目录 |
| Step 0 | `00_brief.md` | → `process/00_brief.md` |
| Step 1（3处） | `01_research.md` | → `process/01_research.md` |
| Step 2（3处） | `01_research.md`、`02_outline.md` | 加 `process/` 前缀 |
| Step 3（6处） | `03_00_前言.md`、`03_chapter_XX_标题.md`、`03_99_附录.md`、`03_draft_full.md` | 加 `process/` 前缀 |
| Step 4（2处） | `04_review_content.md` | → `process/04_review_content.md` |
| Step 5（2处） | `05_review_case.md` | → `process/05_review_case.md` |
| Step 6（2处） | `06_review_style.md` | → `process/06_review_style.md` |
| Step 7.1（1处） | `07_fix.md` | → `process/07_fix.md` |
| Step 7.1b（2处） | `07_fix.md`（命令参数） | → `process/07_fix.md` |
| Step 7.3、7.4 | `07_Final_Tutorial_*.md`、`08_summary.md` | **不变**，留在任务文件夹根 |
