# Workspace Consolidation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 将测试输出从独立的 `test_workspace/` 目录合并到对应教程的 `_workspace/YYYYMMDD_主题/test_YYYYMMDD/` 子目录，删除项目根目录的 `test_workspace/`。

**Architecture:** 修改 `test_framework_integrated.py` 支持外部传入 `test_dir`（跳过内部时间戳子目录创建），同时将 `test_report.md` 写到 `test_dir.parent`（教程根目录）；修改 `testCLAUDE.md` 将 Step 1.1 的 mkdir/cd 指向 `_workspace/XXX/test_YYYYMMDD/` 并更新后续路径引用；最后删除旧 `test_workspace/`。

**Tech Stack:** Python 3.8+, pathlib, 纯文本 Markdown 替换。

---

## Task 1: 修改 test_framework_integrated.py（__init__ + __main__）

**Files:**
- Modify: `tools/test_framework_integrated.py:41-49`（`__init__`）
- Modify: `tools/test_framework_integrated.py:387-408`（`__main__`）

**背景：**
当前 `__init__` 接受 `work_dir` 参数，在其内部自动创建时间戳子目录作为 `test_dir`。
新设计中，testCLAUDE.md 在外部创建 `test_dir`，通过 `--test-dir` 参数传入，脚本直接使用。

**Step 1: 修改 `__init__`**

将（第 41-49 行）：
```python
    def __init__(self, tutorial_path: str, work_dir: str = "./test_workspace"):
        self.tutorial_path = Path(tutorial_path).resolve()
        self.work_dir = Path(work_dir).resolve()
        self.work_dir.mkdir(exist_ok=True)

        # 创建测试目录
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.test_dir = self.work_dir / f"{timestamp}_{self.tutorial_path.stem}"
        self.test_dir.mkdir(exist_ok=True)
```

替换为：
```python
    def __init__(self, tutorial_path: str, work_dir: str = None, test_dir: str = None):
        self.tutorial_path = Path(tutorial_path).resolve() if tutorial_path else Path(".")

        if test_dir:
            # 直接使用外部传入的目录（testCLAUDE.md 新流程）
            self.test_dir = Path(test_dir).resolve()
            self.test_dir.mkdir(exist_ok=True, parents=True)
        else:
            # 兼容旧流程：在 work_dir 下自动创建时间戳子目录
            wd = Path(work_dir or "./test_workspace").resolve()
            wd.mkdir(exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.test_dir = wd / f"{timestamp}_{self.tutorial_path.stem}"
            self.test_dir.mkdir(exist_ok=True)
```

**Step 2: 修改 `__main__` 块**

将（第 387-408 行）：
```python
if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "continue":
        # 继续模式
        test_dir = sys.argv[2] if len(sys.argv) > 2 else None
        if test_dir:
            executor = FullTestExecutor("", work_dir="./test_workspace")
            executor.test_dir = Path(test_dir)
            executor.continue_after_jobs_done()
        else:
            print("用法: python test_framework_integrated.py continue <test_dir>")
    else:
        # 正常模式
        if len(sys.argv) > 1:
            tutorial_path = sys.argv[1]
        else:
            tutorial_path = "test_workspace/elastic_tutorial.md"

        executor = FullTestExecutor(
            tutorial_path=tutorial_path,
            work_dir="./test_workspace"
        )
        executor.run_full_test()
```

替换为：
```python
if __name__ == "__main__":
    # 手动解析 --test-dir 参数（避免引入额外依赖）
    test_dir_arg = None
    remaining = []
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '--test-dir' and i + 1 < len(sys.argv):
            test_dir_arg = sys.argv[i + 1]
            i += 2
        else:
            remaining.append(sys.argv[i])
            i += 1

    if remaining and remaining[0] == "continue":
        # 继续模式：python test_framework_integrated.py continue <test_dir>
        td = remaining[1] if len(remaining) > 1 else test_dir_arg
        if td:
            executor = FullTestExecutor("", test_dir=td)
            executor.continue_after_jobs_done()
        else:
            print("用法: python tools/test_framework_integrated.py continue <test_dir>")
    else:
        # 正常模式：python test_framework_integrated.py <tutorial_path> [--test-dir <test_dir>]
        tutorial_path = remaining[0] if remaining else "_workspace/elastic_tutorial.md"
        executor = FullTestExecutor(
            tutorial_path=tutorial_path,
            test_dir=test_dir_arg
        )
        executor.run_full_test()
```

**Step 3: 验证语法无误**

```bash
python -c "import sys; sys.path.insert(0, 'tools'); from test_framework_integrated import FullTestExecutor; print('import OK')"
```

预期输出：`import OK`

**Step 4: 验证 `__init__` 使用 test_dir 时不创建额外子目录**

```bash
python -c "
import sys, tempfile, pathlib
sys.path.insert(0, 'tools')
from test_framework_integrated import FullTestExecutor
import tempfile, os
td = tempfile.mkdtemp()
e = FullTestExecutor('', test_dir=td)
assert str(e.test_dir) == str(pathlib.Path(td).resolve()), f'test_dir mismatch: {e.test_dir}'
print('test_dir OK:', e.test_dir)
"
```

预期：打印 `test_dir OK: <路径>`，无报错。

**Step 5: Commit**

```bash
git add tools/test_framework_integrated.py
git commit -m "feat: add test_dir param to FullTestExecutor, skip internal timestamp subdir"
```

---

## Task 2: 修改 test_framework_integrated.py（test_report.md 路径）

**Files:**
- Modify: `tools/test_framework_integrated.py`（`phase7_generate_report`，约第 321 行）

**背景：**
当前 `test_report.md` 写到 `self.test_dir` 内（`test_dir/08_test_report.md`）。
新设计：写到 `self.test_dir.parent`（即 `_workspace/XXX/test_report.md`）。

**Step 1: 找到并确认当前行**

```bash
grep -n "test_report\|08_test_report" tools/test_framework_integrated.py
```

预期：找到类似 `report_file = self.test_dir / "08_test_report.md"` 的行。

**Step 2: 修改报告路径**

将：
```python
        report_file = self.test_dir / "08_test_report.md"
```

替换为：
```python
        report_file = self.test_dir.parent / "test_report.md"
```

**Step 3: 验证**

```bash
grep -n "test_report\|08_test_report" tools/test_framework_integrated.py
```

预期：只剩 `self.test_dir.parent / "test_report.md"`，无 `08_test_report`。

**Step 4: Commit**

```bash
git add tools/test_framework_integrated.py
git commit -m "fix: write test_report.md to tutorial root (test_dir.parent) instead of test_dir"
```

---

## Task 3: 修改 testCLAUDE.md — Step 1.1（目录创建逻辑）

**Files:**
- Modify: `testCLAUDE.md`（第 258-266 行，Step 1.1）

**Step 1: 确认当前内容**

```bash
grep -n "test_workspace\|mkdir\|cd test" testCLAUDE.md | head -10
```

预期：找到第 264-265 行的 `mkdir test_workspace/...` 和 `cd test_workspace/...`。

**Step 2: 修改 Step 1.1**

将：
```markdown
**1.1 创建测试工作目录**

```bash
# 生成时间戳目录
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$tutorial_name = "elastic_tutorial"  # 从教程文件名提取
mkdir test_workspace/${timestamp}_${tutorial_name}
cd test_workspace/${timestamp}_${tutorial_name}
```
```

替换为：
```markdown
**1.1 创建测试工作目录**

```bash
# 从教程路径提取父目录，在其中创建带时间戳的 test 子目录
# 假设教程路径为 $tutorial_path（如 "_workspace/20260203_XXX/07_Final_Tutorial_*.md"）
$tutorial_dir = Split-Path $tutorial_path -Parent   # → _workspace/20260203_XXX
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$test_dir = "$tutorial_dir/test_$timestamp"          # → _workspace/20260203_XXX/test_20260226_144949
mkdir $test_dir
# 不再 cd，后续所有命令从项目根目录运行，路径均用 $test_dir 前缀
```
```

**Step 3: 验证**

```bash
grep -n "mkdir test_workspace\|cd test_workspace" testCLAUDE.md
```

预期：无输出（已替换）。

```bash
grep -n "test_dir\|tutorial_dir" testCLAUDE.md | head -5
```

预期：找到新添加的变量定义行。

**Step 4: Commit**

```bash
git add testCLAUDE.md
git commit -m "docs: update testCLAUDE.md Step 1.1 to create test_dir inside _workspace/XXX/"
```

---

## Task 4: 修改 testCLAUDE.md — Step 1.2 和 Step 5.3 工具路径

**Files:**
- Modify: `testCLAUDE.md`（Step 1.2 约第 271 行，Step 5.3 约第 755 行）

**背景：**
之前 `cd test_workspace/XXX/` 后，工具路径用 `../../tools/`（2层上跳到项目根）。
现在从项目根运行，直接用 `tools/`。

**Step 1: 找到所有 `../../tools/` 引用**

```bash
grep -n "../../tools/" testCLAUDE.md
```

预期：找到 Step 1.2 和 Step 5.3 共2处。

**Step 2: 修改 Step 1.2**

将（Step 1.2 的命令行）：
```bash
python ../../tools/test_framework_integrated.py --tutorial ../../_workspace/XXX/07_final.md --phase analyze
```

替换为：
```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir"
```

**Step 3: 修改 Step 5.3**

将（Step 5.3 的命令行）：
```bash
cd Si/02_elastic
python ../../../tools/test_framework_integrated.py --phase postprocess --type elastic
```

替换为：
```bash
python tools/test_framework_integrated.py continue "$test_dir"
```

**Step 4: 修改 Step 2.3（download_files 命令）**

找到 Step 2.3 的 download_files 命令：
```bash
python ../../tools/test_framework_integrated.py --phase download_files --case Si
```

替换为：
```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir"
```

（注：该命令在新流程中由 run_full_test 内部处理，不需要单独调用；如有独立出现则删除或注释。）

**Step 5: 验证无遗留 `../../tools/` 引用**

```bash
grep -n "../../tools/" testCLAUDE.md
```

预期：无输出。

**Step 6: Commit**

```bash
git add testCLAUDE.md
git commit -m "docs: update testCLAUDE.md tool paths from ../../tools/ to tools/"
```

---

## Task 5: 修改 testCLAUDE.md — 其余路径引用（4处）

**Files:**
- Modify: `testCLAUDE.md`（Steps 2.3b, 4.4, 6.4, 工作总结）

**Step 1: 修改 Step 2.3b orbital_fix_report 路径**

找到（Step 2.3b 内）：
```python
report_path = 'orbital_fix_report.json'
```

替换为：
```python
report_path = f'{test_dir}/orbital_fix_report.json'
```

**Step 2: 修改 Step 4.4 示例路径**

找到：
```
- 查看日志：test_workspace/job_logs/job_21984267/
```

替换为：
```
- 查看日志：_workspace/XXX/test_YYYYMMDD/job_logs/job_21984267/
```

**Step 3: 修改 Step 6.4 测试报告模板中的"测试目录"字段**

找到：
```markdown
**测试目录：** test_workspace/20260226_103000_elastic_tutorial/
```

替换为：
```markdown
**测试目录：** _workspace/20260203_105918_弹性常数计算/test_20260226_103000/
```

**Step 4: 修改工作总结中的测试报告路径**

找到：
```
  test_workspace/20260226_103000_elastic_tutorial/test_report.md
```

替换为：
```
  _workspace/20260203_105918_弹性常数计算/test_report.md
```

**Step 5: 验证无遗留 test_workspace 引用（docs/plans 除外）**

```bash
grep -n "test_workspace" testCLAUDE.md
```

预期：无输出。

**Step 6: Commit**

```bash
git add testCLAUDE.md
git commit -m "docs: update testCLAUDE.md remaining test_workspace path references"
```

---

## Task 6: 删除 test_workspace/ 并最终验证

**Files:**
- Delete: `test_workspace/`（项目根目录）

**Step 1: 确认 test_workspace 内容（查看后再删除）**

```bash
ls test_workspace/
```

预期：看到若干旧测试目录（`20260210_xxx`、`20260226_xxx`、`job_logs/`）。

**Step 2: 删除 test_workspace/**

```bash
rm -rf test_workspace/
```

**Step 3: 验证已删除**

```bash
ls test_workspace/ 2>/dev/null || echo "已删除"
```

预期：`已删除`

**Step 4: 整体扫描确认无遗留 test_workspace 引用（排除 docs/ 历史文档）**

```bash
grep -rn "test_workspace" --include="*.py" --include="*.md" . \
  | grep -v "^./docs/" \
  | grep -v "^./test_workspace" \
  | grep -v "^Binary"
```

预期：无输出（testCLAUDE.md 和所有 .py 文件中均无 test_workspace 引用）。

**Step 5: 验证 Python 导入仍正常**

```bash
python -c "import sys; sys.path.insert(0, 'tools'); from test_framework_integrated import FullTestExecutor; print('OK')"
```

预期：`OK`

**Step 6: 更新 .gitignore（如有 test_workspace 条目则删除）**

```bash
grep -n "test_workspace" .gitignore
```

如果有 `test_workspace/` 条目，删除该行（test 数据现在在 `_workspace/` 下，已被 `_workspace/` 的 gitignore 规则覆盖，或者需要在 gitignore 中添加 `_workspace/*/test_*/`）。

检查当前 _workspace 相关的 gitignore 规则：
```bash
grep -n "_workspace\|test_" .gitignore
```

根据输出决定是否需要添加 `_workspace/*/test_*/` 规则（避免测试中间文件入库）。

**Step 7: Commit**

```bash
git add -A
git commit -m "chore: delete test_workspace/, consolidate test outputs under _workspace/"
```
