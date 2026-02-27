# Orbital Runtime Correction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 在测试阶段（testCLAUDE.md Step 2.3）遇到未知错误轨道文件名时，自动查询 GitHub 候选列表，由 Claude 推荐并让用户确认正确文件名，写入 `orbital_db.py`，计算通过后自动修正教程文章。

**Architecture:** `orbital_db.py` 新增 `query_github_candidates()` 和 `add_correction()` 两个函数，提供结构化能力给 Claude 调用。`test_framework_phase3_7_impl.py` 的 `get_file()` 在下载失败时抛出携带候选列表的 `OrbitalNotFoundError`，而不是崩溃。`testCLAUDE.md` 补充 Claude 的处理流程（展示候选、推荐、确认、保存修正报告、测试通过后改教程）。

**Tech Stack:** Python 3.8+, urllib.request (GitHub API), pathlib, json；不引入新依赖。

---

## 背景：当前 get_file() 行为

`tools/test_framework_phase3_7_impl.py` 第 163-233 行：
- 第 167-169 行：已知错误文件名 → 自动映射（`ORBITAL_MAPPING = ORBITAL_CORRECTIONS`）
- 第 197-232 行：尝试从 GitHub 两个 URL 下载，验证不是 HTML/404
- **第 233 行：所有下载源失败 → `raise Exception("下载失败...")`（程序崩溃）**

目标：将第 233 行的 generic Exception 替换为携带候选列表的 `OrbitalNotFoundError`。

---

## Task 1: 向 orbital_db.py 添加 query_github_candidates()

**Files:**
- Modify: `tools/orbital_db.py`（在文件末尾追加）

**Step 1: 在 orbital_db.py 末尾追加函数**

在 `validate_orbital()` 函数之后追加：

```python
def query_github_candidates(element: str) -> list:
    """
    查询 ABACUS GitHub，返回该元素的所有已知 .orb 文件名。

    Args:
        element: 元素符号，如 "Si"、"Ti"

    Returns:
        list of str: 文件名列表，网络失败时返回空列表
    """
    import urllib.request
    import json as _json

    url = "https://api.github.com/repos/deepmodeling/abacus-develop/contents/tests/PP_ORB"
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "AutoTutorial3.0-OrbitalChecker"}
        )
        with urllib.request.urlopen(req, timeout=10) as response:
            data = _json.loads(response.read().decode())
        return [
            item["name"]
            for item in data
            if item["name"].startswith(element + "_") and item["name"].endswith(".orb")
        ]
    except Exception:
        return []
```

**Step 2: 验证函数可导入**

```bash
python -c "from tools.orbital_db import query_github_candidates; print('import OK')"
```

预期输出：`import OK`

**Step 3: 验证网络查询（可选，需要网络）**

```bash
python -c "
from tools.orbital_db import query_github_candidates
result = query_github_candidates('Si')
print('Si candidates:', result[:3])
print('network OK' if result else 'network failed or no Si files')
"
```

预期：打印出包含 `Si_gga_8au_100Ry_2s2p1d.orb` 的列表，或网络失败时打印 `[]`（均为正常）。

**Step 4: Commit**

```bash
git add tools/orbital_db.py
git commit -m "feat: add query_github_candidates() to orbital_db"
```

---

## Task 2: 向 orbital_db.py 添加 add_correction()

**Files:**
- Modify: `tools/orbital_db.py`（在 `query_github_candidates` 之后追加）

**Step 1: 追加 add_correction() 函数**

```python
def add_correction(wrong: str, correct: str):
    """
    向 ORBITAL_CORRECTIONS 添加新映射，同时更新内存状态和源文件。

    Args:
        wrong: 错误的文件名（如 "Si_gga_7au_100Ry_2s2p1d.orb"）
        correct: 正确的文件名（如 "Si_gga_8au_100Ry_2s2p1d.orb"）
    """
    from pathlib import Path as _Path

    # 1. 更新内存状态
    ORBITAL_CORRECTIONS[wrong] = correct
    ALL_KNOWN_ORBITALS.add(correct)

    # 2. 读取源文件
    db_path = _Path(__file__).resolve()
    content = db_path.read_text(encoding='utf-8')

    # 3. 在 ORBITAL_CORRECTIONS 结束的 } 前插入新条目
    #    文件中的标记：}\n\n# 扁平化已知文件名集合
    element = wrong.split("_")[0] if "_" in wrong else "unknown"
    new_entry = f'    # {element}: 运行时自动发现（测试阶段用户确认）\n    "{wrong}": "{correct}",\n'
    marker = '}\n\n# 扁平化已知文件名集合'
    if marker not in content:
        raise RuntimeError("orbital_db.py 格式已变化，无法自动插入新映射，请手动添加。")
    content = content.replace(marker, new_entry + marker, 1)

    # 4. 写回文件
    db_path.write_text(content, encoding='utf-8')
    print(f"[orbital_db] 已写入新映射：{wrong} -> {correct}")
```

**Step 2: 验证函数逻辑（在临时副本上测试，不修改真实文件）**

```bash
python -c "
import shutil, sys
# 备份
shutil.copy('tools/orbital_db.py', 'tools/orbital_db.py.bak')
try:
    from tools.orbital_db import add_correction, ORBITAL_CORRECTIONS
    # 测试写入
    add_correction('Fe_gga_7au_100Ry_4s2p2d1f.orb', 'Fe_gga_9au_100Ry_4s2p2d1f.orb')
    # 验证内存更新
    assert 'Fe_gga_7au_100Ry_4s2p2d1f.orb' in ORBITAL_CORRECTIONS
    print('内存更新 OK')
    # 验证文件写入
    content = open('tools/orbital_db.py', encoding='utf-8').read()
    assert 'Fe_gga_7au_100Ry_4s2p2d1f.orb' in content
    print('文件写入 OK')
finally:
    # 还原
    shutil.copy('tools/orbital_db.py.bak', 'tools/orbital_db.py')
    import os; os.remove('tools/orbital_db.py.bak')
    print('已还原')
"
```

预期输出：
```
[orbital_db] 已写入新映射：Fe_gga_7au_100Ry_4s2p2d1f.orb -> Fe_gga_9au_100Ry_4s2p2d1f.orb
内存更新 OK
文件写入 OK
已还原
```

**Step 3: Commit**

```bash
git add tools/orbital_db.py
git commit -m "feat: add add_correction() to orbital_db for runtime db updates"
```

---

## Task 3: 添加 OrbitalNotFoundError 并修改 get_file()

**Files:**
- Modify: `tools/test_framework_phase3_7_impl.py`

**Step 1: 在文件顶部 import 区域之后（第 19 行之后）添加异常类**

在 `class BohriumJobManager:` 之前添加：

```python
class OrbitalNotFoundError(Exception):
    """
    轨道文件下载失败且无已知映射时抛出。
    携带 GitHub 候选列表，供 Claude 在 testCLAUDE.md 流程中处理。

    Attributes:
        filename: 请求的错误文件名
        candidates: GitHub 返回的该元素候选文件名列表（可能为空）
    """
    def __init__(self, filename: str, candidates: list):
        self.filename = filename
        self.candidates = candidates
        super().__init__(
            f"轨道文件不存在：{filename}\n"
            f"候选文件：{candidates if candidates else '（GitHub 查询失败，请手动提供）'}"
        )
```

**Step 2: 修改 get_file() 第 233 行**

将：
```python
        raise Exception(f"下载失败: {filename} - 所有下载源都失败")
```

替换为：
```python
        # 如果是轨道文件，查询 GitHub 候选列表并抛出结构化异常
        if file_type == "orbital":
            element = filename.split("_")[0] if "_" in filename else ""
            candidates = query_github_candidates(element) if element else []
            raise OrbitalNotFoundError(filename, candidates)
        raise Exception(f"下载失败: {filename} - 所有下载源都失败")
```

同时在文件顶部的 import 区域添加（`from orbital_db import ORBITAL_CORRECTIONS` 之后）：

```python
from orbital_db import ORBITAL_CORRECTIONS, query_github_candidates
```

**Step 3: 验证导入和类定义正确**

```bash
python -c "
from tools.test_framework_phase3_7_impl import OrbitalNotFoundError, PseudopotentialManager
e = OrbitalNotFoundError('Si_gga_7au_100Ry_2s2p1d.orb', ['Si_gga_8au_100Ry_2s2p1d.orb'])
assert e.filename == 'Si_gga_7au_100Ry_2s2p1d.orb'
assert e.candidates == ['Si_gga_8au_100Ry_2s2p1d.orb']
print('OrbitalNotFoundError OK')
print('message:', str(e))
"
```

预期输出：
```
OrbitalNotFoundError OK
message: 轨道文件不存在：Si_gga_7au_100Ry_2s2p1d.orb
候选文件：['Si_gga_8au_100Ry_2s2p1d.orb']
```

**Step 4: Commit**

```bash
git add tools/test_framework_phase3_7_impl.py
git commit -m "feat: raise OrbitalNotFoundError with GitHub candidates instead of crashing"
```

---

## Task 4: 更新 testCLAUDE.md

**Files:**
- Modify: `testCLAUDE.md`

**Step 1: 在 Step 2.3「下载逻辑」末尾添加异常处理子步骤**

找到 Step 2.3 中的「**输出示例：**」代码块（大约第 378-393 行），在其后、Step 2.4 之前插入：

````markdown
**2.3b 处理轨道文件不存在（OrbitalNotFoundError）**

如果下载过程中抛出 `OrbitalNotFoundError`，执行以下流程：

**情况1：GitHub 查询成功，有候选列表**

```
⚠️ 轨道文件不存在：{filename}

[查询] 已从 ABACUS GitHub 获取候选文件：
  1. {candidate_1}  ← 推荐（说明推荐理由，如"8au 替换 7au，命名规律一致"）
  2. {candidate_2}
  ...

请选择正确文件名（输入序号，或直接输入文件名）：
```

**Think Aloud：** 说明为什么推荐某个候选（命名模式相似度、截断能是否一致、轨道数是否匹配）

用户确认后：
```bash
# 写入 orbital_db.py
python -c "
import sys; sys.path.insert(0, 'tools')
from orbital_db import add_correction
add_correction('{wrong}', '{correct}')
"
```

然后重新调用 get_file() 下载正确文件，继续测试流程。

**情况2：GitHub 查询失败（网络问题/API 限速），候选列表为空**

```
⚠️ 轨道文件不存在：{filename}
   GitHub 查询失败，请手动提供正确文件名。

参考：https://github.com/deepmodeling/abacus-develop/tree/develop/tests/PP_ORB
请输入正确的文件名：
```

用户提供文件名后，同样调用 `add_correction()` 写入，然后重试下载。

**每次发生修正后，记录到修正报告：**

```bash
python -c "
import json, os
from datetime import datetime
report_path = 'orbital_fix_report.json'
report = json.loads(open(report_path).read()) if os.path.exists(report_path) else {'tutorial': '', 'corrections': []}
report['tutorial'] = '[当前教程路径]'
report['corrections'].append({
    'wrong': '{wrong}',
    'correct': '{correct}',
    'source': 'user_confirmed',
    'timestamp': datetime.now().isoformat()
})
open(report_path, 'w').write(json.dumps(report, ensure_ascii=False, indent=2))
print('[记录] 修正报告已更新')
"
```
````

**Step 2: 在 Step 6.4「生成测试报告」末尾添加修正报告处理**

找到 Step 6.4 结尾（大约第 814 行之后），在「**Think Aloud：**」之前插入：

````markdown
**6.4b 处理轨道文件修正记录**

检查当前测试目录下是否存在 `orbital_fix_report.json`：

```bash
ls orbital_fix_report.json 2>/dev/null && echo "存在" || echo "不存在"
```

**如果存在且本次测试通过：**

1. 将修正记录追加到 `test_report.md`：

```markdown
## 轨道文件修正记录

本次测试发现并修正了以下轨道文件名错误（已验证计算通过）：

| 错误文件名 | 正确文件名 | 来源 |
|-----------|-----------|------|
| {wrong} | {correct} | 用户确认 |

**已自动修正教程文章** `{tutorial_path}`（见下方修正详情）
```

2. 用 `orbital_validator.py --fix` 修正教程文章：

```bash
python tools/orbital_validator.py "{tutorial_path}" --fix
```

3. 说明修正了哪些行（orbital_validator 的输出）

**如果存在但本次测试不通过：**

将修正尝试记录到报告，但**不修改教程**：

```markdown
## 轨道文件修正记录

本次测试尝试了以下轨道文件名修正，但计算结果不通过，教程未修改：

| 错误文件名 | 尝试修正为 |
|-----------|-----------|
| {wrong} | {correct} |

请检查计算参数或联系教程作者。
```

**如果不存在：** 跳过此步骤。

**Think Aloud：** 说明是否发生了轨道文件修正，如何处理
````

**Step 3: 验证修改位置正确**

```bash
grep -n "OrbitalNotFoundError\|2.3b\|6.4b\|orbital_fix_report" testCLAUDE.md | head -20
```

预期：能看到 `2.3b`、`6.4b`、`OrbitalNotFoundError`、`orbital_fix_report` 四个关键词所在行。

**Step 4: Commit**

```bash
git add testCLAUDE.md
git commit -m "docs: add orbital runtime correction flow to testCLAUDE.md (Step 2.3b, 6.4b)"
```

---

## Task 5: 端到端验证

**Step 1: 验证所有新接口可导入**

```bash
python -c "
import sys; sys.path.insert(0, 'tools')
from orbital_db import query_github_candidates, add_correction, validate_orbital
from test_framework_phase3_7_impl import OrbitalNotFoundError, PseudopotentialManager

# 验证 OrbitalNotFoundError 结构
e = OrbitalNotFoundError('test.orb', ['candidate.orb'])
assert e.filename == 'test.orb'
assert e.candidates == ['candidate.orb']

# 验证 ORBITAL_MAPPING 仍然指向共享数据库
from orbital_db import ORBITAL_CORRECTIONS
assert PseudopotentialManager.ORBITAL_MAPPING is ORBITAL_CORRECTIONS

print('所有接口验证通过')
"
```

预期输出：`所有接口验证通过`

**Step 2: 验证现有功能不受影响**

```bash
python -c "
from tools.orbital_db import validate_orbital
assert validate_orbital('Si_gga_8au_100Ry_2s2p1d.orb') == (True, 'Si_gga_8au_100Ry_2s2p1d.orb', 'OK')
assert validate_orbital('Si_gga_7au_100Ry_2s2p1d.orb')[0] == False
print('现有功能未受影响')
"
```

```bash
python tools/orbital_validator.py "_workspace/20260203_105918_弹性常数计算/07_final.md"
```

预期：与之前相同的输出（2 处 NG，1 处 OK），不受本次修改影响。

**Step 3: 确认 git status 无遗漏**

```bash
git status
```

预期：仅显示预先存在的未跟踪文件（_workspace、data、docs 等），无计划范围内的未提交修改。

---

## 变更范围汇总

| 文件 | 操作 | 内容 |
|------|------|------|
| `tools/orbital_db.py` | **修改** | 新增 `query_github_candidates()`、`add_correction()` |
| `tools/test_framework_phase3_7_impl.py` | **修改** | 新增 `OrbitalNotFoundError` 类；`get_file()` 失败时抛出结构化异常 |
| `testCLAUDE.md` | **修改** | 新增 Step 2.3b（异常处理流程）、Step 6.4b（修正报告与教程修改） |
| `CLAUDE.md` | **不变** | Step 7.1b 保留 |
| `orbital_validator.py` | **不变** | 功能不变，Step 6.4b 中复用 |
