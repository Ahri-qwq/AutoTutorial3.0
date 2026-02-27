# Orbital Runtime Correction Design

**日期：** 2026-02-27
**状态：** 已确认，待实现

---

## 问题背景

当前 `orbital_db.py` 的 `ORBITAL_CORRECTIONS` 只能处理**已知**错误文件名（5条映射）。当 RAG 检索到一个**未知**错误文件名时，`PseudopotentialManager.get_file()` 下载失败后直接报错退出，测试卡死，无法自动恢复。

**目标：** 在 testCLAUDE.md 测试流程中，遇到未知错误轨道文件名时：
1. 自动查 GitHub 获取候选列表
2. 由 Claude 推荐 + 用户确认正确文件名
3. 写入 `orbital_db.py`，本次及未来测试自动处理
4. 继续计算，测试通过后自动修正教程文章

---

## 设计原则

**代码做结构化的事**（确定的、机械的）：
- 检测下载失败（404/HTML 响应）
- 查 GitHub API，返回候选文件名列表
- 提供 `add_correction()` 接口写入数据库
- 生成结构化修正报告

**Claude 做判断性的事**（模糊的、需要上下文的）：
- 从候选列表里推荐最合理的选项（结合教程主题、元素、计算类型）
- 向用户展示选项并说明推荐理由
- 获取确认后调用接口
- 保存修正报告，决定是否继续

---

## 架构概览

```
testCLAUDE.md Step 2.3 下载文件时
        ↓
  已知错误文件名？──YES──→ ORBITAL_CORRECTIONS 自动映射（现有逻辑，不变）
        ↓ NO
    下载失败？──NO──→ 正常继续
        ↓ YES
  Python: 查 GitHub API → 返回候选列表
        ↓
  Claude: 展示候选 + 推荐 → 等用户确认
        ↓
  Python: add_correction() 写入 orbital_db.py
  Python: 保存修正报告 orbital_fix_report.json
        ↓
  用正确文件名重试下载 → 继续计算
        ↓
  计算完成，Step 6 结果对比
        ↓
  通过？──YES──→ 用 orbital_fix_report.json 修改教程文章
       └─NO──→ 报告不通过（修正报告仍保存供参考，不修改教程）
```

**不变的部分：**
- `CLAUDE.md Step 7.1b`：保留，静态扫描已知错误，测试前第一道防线
- `orbital_validator.py`：保留，功能不变
- `ORBITAL_CORRECTIONS` 里已有的映射：继续自动处理，不打扰用户

---

## Python 代码变更

### `tools/orbital_db.py` 新增两个函数

```python
def query_github_candidates(element: str) -> list[str]:
    """
    查询 ABACUS GitHub，返回该元素的所有已知轨道文件名。

    调用：
        GET https://api.github.com/repos/deepmodeling/abacus-develop/contents/tests/PP_ORB

    过滤：文件名以 element + "_" 开头的 .orb 文件。
    失败时（网络不通、API 限速）返回空列表。
    不引入新依赖，使用 urllib.request。
    """

def add_correction(wrong: str, correct: str):
    """
    向 orbital_db.py 的 ORBITAL_CORRECTIONS 写入新映射。

    操作：
    1. 更新内存中的 ORBITAL_CORRECTIONS 和 ALL_KNOWN_ORBITALS
    2. 读取 orbital_db.py 源文件
    3. 在 ORBITAL_CORRECTIONS 字典末尾插入新条目（含注释说明来源）
    4. 写回文件
    """
```

### `tools/test_framework_phase3_7_impl.py` 修改 `get_file()`

当前下载失败行为：直接报错退出。

改为：
```
下载内容 → 检测是否为 HTML（404 页面）
  → 如果是：
      candidates = query_github_candidates(element)
      raise OrbitalNotFoundError(filename, candidates)
  → 如果不是：正常返回
```

`OrbitalNotFoundError` 携带：
- `filename`：请求的错误文件名
- `candidates`：GitHub 返回的候选列表（可能为空）

Claude 捕获这个异常并接管后续处理，程序不崩溃。

### 新增文件 `orbital_fix_report.json`（测试运行时生成）

保存在测试工作目录下：

```json
{
  "tutorial": "_workspace/XXX/07_final.md",
  "corrections": [
    {
      "wrong": "Si_gga_7au_100Ry_2s2p1d.orb",
      "correct": "Si_gga_8au_100Ry_2s2p1d.orb",
      "source": "user_confirmed",
      "timestamp": "2026-02-27T10:30:00"
    }
  ]
}
```

`source` 字段取值：`"user_confirmed"`（用户选择）或 `"user_manual"`（用户手动输入）。

---

## testCLAUDE.md 流程变更

### Step 2.3 新增：遇到 `OrbitalNotFoundError` 时

```
⚠️ 轨道文件不存在：Si_gga_7au_100Ry_2s2p1d.orb

[查询] 正在查询 ABACUS GitHub...
  Si 元素的已知轨道文件：
    1. Si_gga_8au_100Ry_2s2p1d.orb  ← 推荐（8au 替换 7au，命名规律一致）
    2. Si_gga_8au_60Ry_2s2p1d.orb
    3. （其他候选）

请选择正确文件名（输入序号或直接输入文件名）：
```

用户确认后：
```
[写入] 新映射已保存到 tools/orbital_db.py
[记录] 保存到 orbital_fix_report.json
[重试] 下载 Si_gga_8au_100Ry_2s2p1d.orb ... ✅ 继续
```

若 GitHub 查询失败（无网络/API 限速）：
```
⚠️ GitHub 查询失败，请手动提供正确文件名：
```

**Think Aloud：** 说明为什么推荐某个候选（命名模式、元素匹配、截断能一致性等）

### Step 6.4 生成报告时新增

若本次测试中发生过轨道文件修正（`orbital_fix_report.json` 非空）：

**测试通过时：**
```markdown
## 轨道文件修正记录

本次测试发现并修正了以下轨道文件名错误：

| 错误文件名 | 正确文件名 | 来源 |
|-----------|-----------|------|
| Si_gga_7au_100Ry_2s2p1d.orb | Si_gga_8au_100Ry_2s2p1d.orb | 用户确认 |

**已自动修正教程文章** `07_final.md`：
- 第 385 行：Si_gga_7au → Si_gga_8au
```

**测试不通过时：**
```markdown
## 轨道文件修正记录

本次测试尝试了以下轨道文件名修正，但计算结果不通过，教程未修改：

| 错误文件名 | 尝试修正为 |
|-----------|-----------|
| Si_gga_7au_100Ry_2s2p1d.orb | Si_gga_8au_100Ry_2s2p1d.orb |

请检查计算参数或联系教程作者。
```

---

## 变更范围汇总

| 文件 | 操作 | 内容 |
|------|------|------|
| `tools/orbital_db.py` | **修改** | 新增 `query_github_candidates()`、`add_correction()` |
| `tools/test_framework_phase3_7_impl.py` | **修改** | `get_file()` 失败时抛出 `OrbitalNotFoundError` |
| `testCLAUDE.md` | **修改** | Step 2.3 加处理流程、Step 6.4 加修正报告逻辑 |
| `CLAUDE.md` | **不变** | Step 7.1b 保留 |
| `orbital_validator.py` | **不变** | 功能不变 |

**不引入新依赖**（GitHub API 调用使用标准库 `urllib.request`）。

---

## 后续可扩展方向（不在本计划范围）

- **定期同步 GitHub 列表到 `KNOWN_ORBITALS`**：自动保持数据库更新
- **为赝势文件实现同样的 runtime correction**：当前只覆盖轨道文件
- **批量处理多个未知文件**：当前设计为逐个处理
