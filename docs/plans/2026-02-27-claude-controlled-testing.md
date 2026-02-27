# Claude-Controlled Testing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 让 Python 框架只负责"文件准备"，Claude 完全控制任务提交决策；同时在生成和测试两端双重拦截 `nbands auto` 等兼容性问题，并在测试结束后将修复建议反向写回教程原文。

**Architecture:** 新增 `--phase prepare` 参数使 Python 框架仅运行 Phase 1+2（解析+准备文件），不提交不等待；扩展 `orbital_validator.py` 同时检测 INPUT 参数问题（如 `nbands auto`）；testCLAUDE.md 新增 Step 2.5（INPUT 参数检查）和 Step 7（修复反向写回教程），并修正 Step 1.2/2.3 的双重执行 Bug。

**Tech Stack:** Python 3.8+, pathlib, re, argparse，纯文本 Markdown 替换。

---

## Task 1: 在 test_framework_integrated.py 中新增 `run_prepare_only()` + `--phase prepare`

**Files:**
- Modify: `tools/test_framework_integrated.py`（`FullTestExecutor` 类 + `__main__` 块）

**背景：**
当前 `testCLAUDE.md` Step 1.2 调用的是 `run_full_test()`，它会自动提交任务（Phase 3）和监控（Phase 4），导致 Claude 的手工提交步骤（Step 3）重复提交。需要新增一个只跑 Phase 1+2 的方法，并通过 `--phase prepare` 参数触发。

---

**Step 1: 在 `FullTestExecutor` 类中新增 `run_prepare_only()` 方法**

定位文件 `tools/test_framework_integrated.py`，在 `run_full_test()` 方法（约第 77 行）**之后**、`continue_after_jobs_done()` 方法**之前**插入：

```python
    def run_prepare_only(self):
        """仅运行 Phase 1+2：解析教程 + 准备输入文件（不提交、不监控）"""
        self.start_time = datetime.now()

        print("\n" + "=" * 70)
        print("AutoTutorial 3.0 - 准备输入文件（--phase prepare）")
        print("=" * 70)

        try:
            # Phase 1: 分析文章
            self.phase1_analyze()

            # Phase 2: 准备输入文件
            self.phase2_prepare_inputs()

            print("\n" + "=" * 70)
            print("[OK] 准备完成，请按教程内容决定任务提交顺序")
            print("=" * 70)
            print(f"\n测试目录：{self.test_dir}")
            print(f"解析结果：{self.test_dir / '01_analysis.json'}")
            print("\n下一步：")
            print("  1. 读取教程内容，理解计算流程和任务依赖关系")
            print("  2. 按顺序用 bohr job submit 提交任务（串行/并行由你决定）")
            print(f"  3. 所有任务完成后运行：")
            print(f"     python tools/test_framework_integrated.py continue {self.test_dir}")

        except Exception as e:
            print(f"\n[ERROR] 准备失败: {e}")
            import traceback
            traceback.print_exc()
```

---

**Step 2: 更新 `__main__` 块，支持解析 `--phase` 参数**

定位 `__main__` 块（约第 392 行），将整个 `if __name__ == "__main__":` 块替换为：

```python
if __name__ == "__main__":
    # 手动解析 --test-dir 和 --phase 参数（避免引入额外依赖）
    test_dir_arg = None
    phase_arg = None
    remaining = []
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '--test-dir' and i + 1 < len(sys.argv):
            test_dir_arg = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--phase' and i + 1 < len(sys.argv):
            phase_arg = sys.argv[i + 1]
            i += 2
        else:
            remaining.append(sys.argv[i])
            i += 1

    if remaining and remaining[0] == "continue":
        # 继续模式：python tools/test_framework_integrated.py continue <test_dir>
        td = remaining[1] if len(remaining) > 1 else test_dir_arg
        if td:
            executor = FullTestExecutor("", test_dir=td)
            executor.continue_after_jobs_done()
        else:
            print("用法: python tools/test_framework_integrated.py continue <test_dir>")
    else:
        # 正常模式：python tools/test_framework_integrated.py <tutorial_path> [--test-dir <dir>] [--phase prepare]
        tutorial_path = remaining[0] if remaining else "_workspace/elastic_tutorial.md"
        executor = FullTestExecutor(
            tutorial_path=tutorial_path,
            test_dir=test_dir_arg
        )
        if phase_arg == "prepare":
            executor.run_prepare_only()
        else:
            executor.run_full_test()
```

---

**Step 3: 验证语法正确**

```bash
python -c "import sys; sys.path.insert(0, 'tools'); from test_framework_integrated import FullTestExecutor; print('import OK')"
```

预期输出：`import OK`

---

**Step 4: 验证 `run_prepare_only` 方法存在**

```bash
python -c "import sys; sys.path.insert(0, 'tools'); from test_framework_integrated import FullTestExecutor; e = FullTestExecutor(''); print(hasattr(e, 'run_prepare_only'))"
```

预期输出：`True`

---

**Step 5: Commit**

```bash
git add tools/test_framework_integrated.py
git commit -m "feat: add --phase prepare to run Phase 1+2 only, defer job submission to Claude"
```

---

## Task 2: 扩展 `orbital_validator.py` 检测 INPUT 参数兼容性问题

**Files:**
- Modify: `tools/orbital_validator.py`

**背景：**
当前 `orbital_validator.py` 只检测 `.orb` 轨道文件名。需要扩展：同时检测 Markdown 代码块中的 `nbands auto`（ABACUS v3.10.x 不支持），并在 `--fix` 模式下自动删除该行。

---

**Step 1: 在 `find_orbital_filenames()` 之后（约第 36 行）插入新函数**

```python
def find_input_param_issues(content: str):
    """
    从 Markdown 内容中检测 INPUT 参数兼容性问题。
    Returns: list of (param_desc, line_number, fix_action)
      fix_action: 'remove_line' = 删除该行
    """
    results = []
    for i, line in enumerate(content.splitlines(), start=1):
        # nbands auto - ABACUS v3.10.x 不支持，删除后 ABACUS 自动计算
        if re.search(r'\bnbands\s+auto\b', line, re.IGNORECASE):
            results.append(('nbands auto', i, 'remove_line'))
    return results
```

---

**Step 2: 在 `validate_tutorial()` 函数的 `orb_files` 检查代码块之后（约第 88 行，`return issues, fixed_count` 之前）插入 INPUT 参数检查逻辑**

找到 `validate_tutorial()` 中以下行：
```python
    return issues, fixed_count
```

在这一行**之前**插入：

```python
    # ── INPUT 参数兼容性检查 ──────────────────────────────────────────
    param_issues = find_input_param_issues(content)

    if param_issues:
        print(f"\n[扫描] 发现 {len(param_issues)} 处 INPUT 参数兼容性问题：\n")
        if fix:
            # 按行号倒序删除，避免行号偏移
            lines = fixed_content.splitlines(keepends=True)
            lines_to_remove = set()
            for param_desc, line_no, action in param_issues:
                print(f"  第{line_no}行  NG  {param_desc}")
                print(f"             -> 删除此行（ABACUS 自动计算，无需指定）")
                if action == 'remove_line':
                    lines_to_remove.add(line_no)
            if lines_to_remove:
                fixed_content = ''.join(
                    line for i, line in enumerate(lines, start=1)
                    if i not in lines_to_remove
                )
                fixed_count += len(lines_to_remove)
                # 若已有轨道修正，fixed_content 在前面已替换，需重新应用行删除
                # 注：当 fix=True 且两类问题同时存在时，行号基于原始 content 计算
                # 轨道替换不影响行号，因此这里直接删除即可
        else:
            for param_desc, line_no, action in param_issues:
                print(f"  第{line_no}行  NG  {param_desc}")
                print(f"             -> 删除此行（ABACUS 自动计算，无需指定）")
            print(f"\n[摘要] 发现 {len(param_issues)} 处 INPUT 参数问题，均可自动修正。")
            print(f"       使用 --fix 参数自动修正。")

    if fix and fixed_count > 0:
```

**注意：** 上面的插入需要同时将原来的 `if fix and fixed_count > 0:` 代码块（原 `return` 之前的最后几行）**去掉重复**，因为我们把它合并进来了。具体地，找到原来的：

```python
    if fix and fixed_count > 0:
        out_path = Path(output_path) if output_path else path
        out_path.write_text(fixed_content, encoding='utf-8')
        print(f"\n[修正] 已自动修正 {fixed_count} 处，保存到：{out_path}")
    elif issues:
        print(f"\n[摘要] 发现 {len(issues)} 处问题，{sum(1 for _, _, c, _ in issues if c)} 处可自动修正。")
        print(f"       使用 --fix 参数自动修正可修正的问题。")

    return issues, fixed_count
```

替换为（已包含 param_issues 之后的最终保存逻辑）：

```python
    # ── INPUT 参数兼容性检查 ──────────────────────────────────────────
    param_issues = find_input_param_issues(content)

    if param_issues:
        print(f"\n[扫描] 发现 {len(param_issues)} 处 INPUT 参数兼容性问题：\n")
        if fix:
            lines = fixed_content.splitlines(keepends=True)
            lines_to_remove = set()
            for param_desc, line_no, action in param_issues:
                print(f"  第{line_no}行  NG  {param_desc}")
                print(f"             -> 删除此行（ABACUS 自动计算，无需指定）")
                if action == 'remove_line':
                    lines_to_remove.add(line_no)
            if lines_to_remove:
                fixed_content = ''.join(
                    line for i, line in enumerate(lines, start=1)
                    if i not in lines_to_remove
                )
                fixed_count += len(lines_to_remove)
        else:
            for param_desc, line_no, action in param_issues:
                print(f"  第{line_no}行  NG  {param_desc}")
                print(f"             -> 删除此行（ABACUS 自动计算，无需指定）")
            print(f"\n[摘要] 发现 {len(param_issues)} 处 INPUT 参数问题，均可自动修正。")
            print(f"       使用 --fix 参数自动修正。")

    # ── 最终保存 ─────────────────────────────────────────────────────
    if fix and fixed_count > 0:
        out_path = Path(output_path) if output_path else path
        out_path.write_text(fixed_content, encoding='utf-8')
        print(f"\n[修正] 已自动修正 {fixed_count} 处，保存到：{out_path}")
    elif issues or param_issues:
        fixable = sum(1 for _, _, c, _ in issues if c) + len(param_issues)
        print(f"\n[摘要] 发现问题：轨道文件 {len(issues)} 处，INPUT参数 {len(param_issues)} 处。")
        print(f"       共 {fixable} 处可自动修正，使用 --fix 参数修正。")

    return issues, fixed_count
```

---

**Step 3: 验证：创建临时测试文件，检测 `nbands auto`**

```bash
python -c "
import sys, tempfile, pathlib
sys.path.insert(0, 'tools')
from orbital_validator import validate_tutorial

# 创建临时 markdown
content = '''
## 示例

\`\`\`
INPUT_PARAMETERS
nbands    auto
ecutwfc   100
\`\`\`
'''
tmp = pathlib.Path(tempfile.mktemp(suffix='.md'))
tmp.write_text(content, encoding='utf-8')

issues, fixed = validate_tutorial(str(tmp), fix=False)
print('param detected:', fixed == 0)  # fix=False 时不应修改文件

tmp.unlink()
print('DONE')
"
```

预期输出包含 `NG  nbands auto` 和 `DONE`。

---

**Step 4: 验证 `--fix` 模式删除 `nbands auto`**

```bash
python -c "
import sys, tempfile, pathlib
sys.path.insert(0, 'tools')
from orbital_validator import validate_tutorial

content = 'nbands    auto\necutwfc   100\n'
tmp = pathlib.Path(tempfile.mktemp(suffix='.md'))
tmp.write_text(content, encoding='utf-8')

issues, fixed = validate_tutorial(str(tmp), fix=True)
result = tmp.read_text(encoding='utf-8')
assert 'nbands' not in result, f'nbands still present: {result!r}'
assert 'ecutwfc' in result, f'ecutwfc was removed: {result!r}'
print('fix OK')
tmp.unlink()
"
```

预期输出：`fix OK`

---

**Step 5: Commit**

```bash
git add tools/orbital_validator.py
git commit -m "feat: extend orbital_validator to detect and remove 'nbands auto' (ABACUS v3.10.x incompatible)"
```

---

## Task 3: 更新 `testCLAUDE.md`（4处修改）

**Files:**
- Modify: `testCLAUDE.md`

**背景：**
- Step 1.2 需要加 `--phase prepare`（防止自动提交任务）
- Step 2.3 需要改成说明（去掉重复的 Python 命令）
- 新增 Step 2.5（INPUT 参数兼容性检查）
- 新增 Step 7（修复反向写回教程）

---

**Step 1: 修改 Step 1.2 命令，加 `--phase prepare`**

找到（约第 273 行）：
```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir"
```

替换为：
```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir" --phase prepare
```

---

**Step 2: 替换 Step 2.3（去掉重复命令，改成说明文字）**

找到（约第 362-395 行）：
```markdown
**2.3 下载赝势和轨道文件**

使用内置的文件管理器自动下载：
```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir"
```

**下载逻辑：**
1. 检查本地缓存（`tools/pseudopotentials/`、`tools/orbitals/`）
2. 如已缓存，复制到案例目录
3. 如未缓存，从ABACUS仓库下载：
   - 赝势：`http://abacus.deepmodeling.com/...[文件名].upf`
   - 轨道：`http://abacus.deepmodeling.com/...[文件名].orb`
4. 验证下载内容（检查HTML/404错误）
5. 应用文件映射表（如`Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`）

**输出示例：**
```
=== 下载依赖文件 ===
...
✅ 所有依赖文件已准备
```
```

替换为：
```markdown
**2.3 下载赝势和轨道文件**

Step 1.2 的 `--phase prepare` 命令已自动完成以下操作：
1. 从本地缓存（`tools/pseudopotentials/`、`tools/orbitals/`）复制已缓存文件
2. 若未缓存，从 ABACUS 仓库下载并验证（检查 HTML/404 错误）
3. 自动应用文件映射表（如 `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`）

确认 Step 1.2 输出中显示"✅ 所有依赖文件已准备"后继续。若出现 `OrbitalNotFoundError`，按 Step 2.3b 处理。
```

---

**Step 3: 在 Step 2.4（job.json）和 Step 2.5（STRU修复）之间，新增 Step 2.5（INPUT参数检查）**

找到（约第 459-503 行，Step 2.4 和 Step 2.5 之间）：

```markdown
**2.5 修复STRU文件格式（如需要）**
```

在这一行**之前**插入新的 **Step 2.5** 节（并将原 **2.5** 改为 **2.6**）：

```markdown
**2.5 检查 INPUT 参数兼容性**

扫描已准备的 INPUT 文件，查找已知不兼容参数：
```bash
grep -rn "nbands.*auto" "$test_dir/"
```

**若发现 `nbands auto`：**

1. 修改对应 INPUT 文件，删除 `nbands auto` 这一行（ABACUS 未指定时自动计算，效果相同）

2. 将修改记录追加到 `"$test_dir/param_fix_report.json"`：
```bash
python -c "
import json, os
from datetime import datetime
report_path = '$test_dir/param_fix_report.json'
report = json.loads(open(report_path).read()) if os.path.exists(report_path) else {'tutorial': '', 'fixes': []}
report['tutorial'] = '$tutorial_path'
report['fixes'].append({
    'type': 'input_param',
    'param': 'nbands auto',
    'action': 'removed',
    'reason': 'ABACUS v3.10.x不支持auto关键字，删除后ABACUS自动计算',
    'location': '教程对应代码块（已修改INPUT文件）',
    'timestamp': datetime.now().isoformat()
})
open(report_path, 'w', encoding='utf-8').write(json.dumps(report, ensure_ascii=False, indent=2))
print('[记录] param_fix_report.json 已更新')
"
```

**Think Aloud：** 说明发现了哪些参数问题，如何处理，是否影响后续计算

**若未发现任何问题：** 跳过此步骤。

---

**2.6 修复STRU文件格式（如需要）**
```

**注意：** 将原来的 `**2.5 修复STRU文件格式（如需要）**` 改为 `**2.6 修复STRU文件格式（如需要）**`，以及后续 Think Aloud 中的序号对应更新（Step 2 完成标志的序号不影响）。

---

**Step 4: 在"工作总结"章节（约第 1020 行）之前，插入 Step 7**

找到（约第 1018-1020 行）：
```markdown
---

## 工作总结
```

在这一段**之前**插入：

```markdown
---

### Step 7: 将测试发现的问题反向修正教程原文

**目标：** 将测试过程中发现的轨道文件名错误、INPUT 参数问题等，写回教程原文，避免下次测试时重复出现相同问题。

**前提条件：** 仅当 Step 6 测试结果为 **全部 PASS** 时执行。若测试失败，只记录不修改（避免用错误参数修改教程）。

**执行：**

**7.1 检查所有修正记录**

```bash
ls "$test_dir/orbital_fix_report.json" 2>/dev/null && echo "轨道修正记录：存在" || echo "轨道修正记录：无"
ls "$test_dir/param_fix_report.json"   2>/dev/null && echo "参数修正记录：存在" || echo "参数修正记录：无"
```

**7.2 汇总并 Think Aloud**

读取所有存在的 fix report，输出待修改清单：

```
[汇总] 本次测试发现以下需回写教程的修正：

  轨道文件名：
    - Si_gga_7au_100Ry_2s2p1d.orb → Si_gga_8au_100Ry_2s2p1d.orb（第385行）

  INPUT 参数：
    - nbands auto → 删除（ABACUS 自动计算）（第47行代码块）
```

若两个 report 均不存在，输出"本次测试无需回写教程"并跳过 7.3-7.4。

**7.3 修正教程原文**

```bash
python tools/orbital_validator.py "$tutorial_path" --fix
```

`orbital_validator.py` 已同时处理：
- 轨道文件名替换（来自 orbital_fix_report + orbital_db）
- `nbands auto` 删除

**7.4 Think Aloud：说明修改结果**

- 说明修改了哪些行
- 判断是否需要重走 CLAUDE.md 审查流程：
  - 只删除 `nbands auto`：**无需**重走审查（参数删除不影响文章内容）
  - 轨道文件名替换：**建议**重走 Step 5（案例审查），确保文件名一致性

**完成标志：**
- 所有 fix report 中的修正已应用到教程
- 输出："✅ Step 7完成，教程原文已更新"

---
```

---

**Step 5: 验证 Step 1.2 已包含 `--phase prepare`**

```bash
grep -n "phase prepare" testCLAUDE.md
```

预期：找到包含 `--phase prepare` 的行。

**Step 6: 验证 Step 2.3 无重复的 Python 命令**

```bash
grep -c "test_framework_integrated.py.*tutorial_path.*test-dir\"$" testCLAUDE.md
```

预期：输出 `1`（只剩 Step 1.2 一处，Step 2.3 已改为说明文字）。

**Step 7: 验证 Step 7 和 param_fix_report 已存在**

```bash
grep -n "param_fix_report\|Step 7" testCLAUDE.md | head -10
```

预期：找到 Step 7 标题行和 param_fix_report 相关内容。

---

**Step 8: Commit**

```bash
git add testCLAUDE.md
git commit -m "docs: refactor testCLAUDE.md - --phase prepare, Step 2.5 nbands check, Step 7 feedback loop"
```

---

## Task 4: 更新 `CLAUDE.md` Step 7.1b 的说明文字

**Files:**
- Modify: `CLAUDE.md`

**背景：**
`orbital_validator.py` 现在同时检测轨道文件名和 INPUT 参数问题（`nbands auto`），但 CLAUDE.md 的说明只提到轨道文件名。需要更新说明文字。

---

**Step 1: 找到 Step 7.1b 说明（约第 298-314 行）**

```bash
grep -n "7\.1b\|orbital_validator\|轨道文件名自动" CLAUDE.md
```

---

**Step 2: 更新 Step 7.1b 的标题和说明**

找到：
```markdown
**7.1b 轨道文件名自动验证（如教程含 LCAO 计算）**

如果教程涉及 LCAO 基组（即包含 `.orb` 轨道文件），执行：
```bash
python tools/orbital_validator.py process/07_fix.md --fix
```

- 有 NG 问题并自动修正时：确认修正内容合理，更新 process/07_fix.md
- 有 ?? 问题无法自动修正时：手动查阅 ABACUS GitHub 确认正确文件名
- 全部 OK 时：直接继续 7.2
```

替换为：
```markdown
**7.1b 自动验证：轨道文件名 + INPUT 参数兼容性（如教程含 LCAO 计算）**

如果教程涉及 LCAO 基组（即包含 `.orb` 轨道文件），执行：
```bash
python tools/orbital_validator.py process/07_fix.md --fix
```

此命令同时检测并修正：
1. **轨道文件名错误**（如 `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`）
2. **INPUT 参数兼容性问题**（如 `nbands auto` → 删除，ABACUS v3.10.x 不支持）

处理规则：
- 有 NG/自动修正时：确认修正内容合理，更新 `process/07_fix.md`
- 有 ?? 无法自动修正时：手动查阅 ABACUS GitHub 确认正确文件名
- 全部 OK 时：直接继续 7.2
```

---

**Step 3: 验证**

```bash
grep -A 10 "7\.1b" CLAUDE.md | head -12
```

预期：看到新的说明文字，包含"INPUT 参数兼容性问题"和"nbands auto"。

---

**Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md Step 7.1b - orbital_validator now also checks INPUT param compatibility"
```

---

## Task 5: 整体验证

**Step 1: Python 导入测试**

```bash
python -c "
import sys
sys.path.insert(0, 'tools')
from test_framework_integrated import FullTestExecutor
from orbital_validator import validate_tutorial, find_input_param_issues
print('all imports OK')
"
```

预期：`all imports OK`

**Step 2: 确认 testCLAUDE.md 中无遗留问题**

```bash
# 1. Step 1.2 含 --phase prepare
grep -n "phase prepare" testCLAUDE.md

# 2. Step 2.3 不再有重复的 run_full_test 命令
grep -n "test_framework_integrated.py.*tutorial_path.*--test-dir" testCLAUDE.md

# 3. Step 7 存在
grep -n "Step 7" testCLAUDE.md

# 4. param_fix_report 存在
grep -n "param_fix_report" testCLAUDE.md
```

预期：
- 第1条：找到至少1行
- 第2条：只找到1行（Step 1.2，Step 2.3 已改为说明文字）
- 第3条：找到 Step 7 标题
- 第4条：找到至少2行

**Step 3: 确认 CLAUDE.md Step 7.1b 已更新**

```bash
grep -n "INPUT 参数兼容性\|nbands" CLAUDE.md
```

预期：找到新添加的说明。

**Step 4: Final Commit（如 Task 3/4 已分开提交，此步跳过）**

若所有改动已分步提交，只需运行：
```bash
git log --oneline -5
```

确认看到 4 个新 commit。
