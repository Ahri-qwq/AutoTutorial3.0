# Workspace Consolidation Design

**日期：** 2026-02-27
**状态：** 已确认，待实现

---

## 问题背景

项目根目录同时存在两个顶层目录：
- `_workspace/YYYYMMDD_主题/` — CLAUDE.md 生成的教程文件
- `test_workspace/YYYYMMDD_tutorial/` — testCLAUDE.md 每次测试运行的目录

每个教程对应的生成文件和测试文件分散在两处，可读性差。

## 目标结构

```
_workspace/20260203_105918_弹性常数计算/
├── 07_Final_Tutorial_弹性常数计算.md   ← 打开即见
├── 08_summary.md                        ← 打开即见
├── test_report.md                       ← 测试报告（测试完成后写入）
├── process/                             ← 生成过程文件（现有）
│   ├── 00_brief.md
│   ├── 01_research.md
│   └── ...
└── test_20260226_144949/                ← 每次测试一个带时间戳子目录
    ├── 01_analysis.json
    ├── Si/
    │   ├── 01_relax/
    │   └── 02_elastic/
    └── job_logs/
        ├── job_22134839/
        └── job_22134840/
```

**不再有** `test_workspace/` 项目根目录（现有旧数据删除）。

## 设计决策

- **test 子目录命名：** `test_YYYYMMDD_HHMMSS/`（每次测试独立子目录，保留历史记录）
- **test_report.md 位置：** 写到教程文件夹根（`_workspace/XXX/test_report.md`），不在 test 子目录内
- **旧数据：** 现有 `test_workspace/` 直接删除，不迁移
- **范围：** 只影响未来新建测试，现有 `_workspace/` 文件夹不做整理

## 变更列表

### testCLAUDE.md（6处）

| 位置 | 旧内容 | 新内容 |
|------|--------|--------|
| Step 1.1（2行） | `mkdir test_workspace/${timestamp}_${tutorial_name}` + `cd ...` | 从教程路径提取父目录，`mkdir _workspace/XXX/test_${timestamp}`；从项目根运行 |
| Step 1.2, Step 5.3（各1行） | `python ../../tools/test_framework_integrated.py ...` | `python tools/test_framework_integrated.py ...` |
| Step 2.3b（1行） | `report_path = 'orbital_fix_report.json'` | `report_path = '$test_dir/orbital_fix_report.json'` |
| Step 4.4 示例（1行） | `test_workspace/job_logs/job_21984267/` | `_workspace/XXX/test_YYYYMMDD/job_logs/job_21984267/` |
| Step 6.4 示例模板（1行） | `**测试目录：** test_workspace/20260226_103000_elastic_tutorial/` | `**测试目录：** _workspace/XXX/test_20260226_103000/` |
| 工作总结（1行） | `test_workspace/.../test_report.md` | `_workspace/XXX/test_report.md` |

**不变：** Step 5.1 的 `../../job_logs`（相对于 `$test_dir/Si/01_relax/`，深度不变）、Step 3.2 的 `cd Si/01_relax`

### test_framework_integrated.py（3处）

| 位置 | 变更 |
|------|------|
| `__init__`（约5行） | 新增 `test_dir: str = None` 参数；若提供则直接使用，不再内部创建时间戳子目录 |
| `phase7_generate_report`（1行） | `self.test_dir / "08_test_report.md"` → `self.test_dir.parent / "test_report.md"` |
| `__main__` 块（约3行） | continue 模式改用 `test_dir=` 参数；更新默认 tutorial_path 示例 |

### 项目根目录

- 删除整个 `test_workspace/`（旧数据，不迁移）

### 不修改的文件

- `docs/` 内各文件（历史调试记录，非运行指令）
- `fix_stru.py`（help 文本路径示例，不影响功能）
