# Claude-Controlled Testing Design

**日期：** 2026-02-27
**状态：** 已确认，待实现

---

## 问题背景

当前 testCLAUDE.md 与 Python 框架存在三个核心冲突：

1. **双重执行**：Step 1.2 和 Step 2.3 都调用 `run_full_test()`，导致整个流程跑两遍、任务重复提交
2. **控制权模糊**：Python 框架（elastic_plugin 4阶段串行、band_plugin 同步提交）与 testCLAUDE.md 手工 bohr 命令冲突
3. **nbands auto 未拦截**：历史上导致 `stoi` 崩溃的参数，在生成和测试两端均无自动检测

根本设计原则确认：**Claude 读教程理解物理逻辑，决定所有任务提交策略；Python 框架只做文件准备工具。**

---

## 目标架构

```
Claude-controlled testing flow:

  Step 1.2: python --phase prepare    ← Python 准备文件，输出 analysis.json
      ↓
  Step 2.5: Claude 检查 INPUT 参数    ← 发现 nbands auto 等问题 → param_fix_report.json
      ↓
  Step 3: Claude 读教程 → Think Aloud → 手工 bohr job submit（串行/并行由Claude决定）
      ↓
  Step 4: Claude 监控（bohr job status）
      ↓
  Step 5: Claude 下载 + python continue（post-process）
      ↓
  Step 6: 对比结果，生成 test_report.md
      ↓
  Step 7 (新): 读取所有 fix report → 反向修正教程原文（仅测试通过时）
```

---

## 变更列表

### 1. `tools/test_framework_integrated.py`

新增 `--phase prepare` 参数：
- 只运行 Phase 1（解析 → analysis.json）+ Phase 2（下载赝势/轨道，写 INPUT/STRU/KPT）
- 退出，不提交，不监控
- `continue` 模式不变（Phase 5-7）
- 默认无参数行为不变（兼容旧用法）

`analysis.json` 输出格式：
```json
{
  "tutorial_path": "...",
  "case_name": "Si",
  "calc_types": ["relax", "elastic"],
  "prepared_dirs": {
    "relax": "...Si/01_relax/",
    "elastic": "...Si/02_elastic/"
  },
  "orbitals_used": ["Si_gga_8au_100Ry_2s2p1d.orb"],
  "pseudopotentials_used": ["Si_ONCV_PBE-1.0.upf"],
  "expected_results": {}
}
```

### 2. `tools/orbital_validator.py`

扩展检测范围，同时处理：
- 轨道文件名错误（已有）
- INPUT 参数兼容性问题（新增）：
  - `nbands\s+auto` → 删除该行（ABACUS 自动计算）
  - 后续可扩展其他不兼容参数

`--fix` 模式下自动修正，并在 stdout 输出修改摘要。

### 3. `CLAUDE.md` Step 7.1b

说明更新：orbital_validator 现在同时检测轨道文件名 **和** INPUT 参数兼容性问题（无需改命令，改说明）。

### 4. `testCLAUDE.md`（主要改动）

| 位置 | 变化 |
|------|------|
| Step 1.2 | 命令改为 `--phase prepare` |
| Step 2.3 | **删除** |
| Step 2.5（新） | INPUT 参数兼容性检查（grep + param_fix_report.json） |
| Step 3 | 重写：Claude Think Aloud 分析依赖，手工 bohr 提交 |
| Step 7（新） | 读取所有 fix report → 反向修正教程原文 |

### Step 7 详细设计

```
Step 7: 将测试发现的问题反向修正教程

条件：仅当本次测试通过（Step 6 结果为 PASS）时才修改教程

7.1 检查所有修正记录
    ls "$test_dir/orbital_fix_report.json"
    ls "$test_dir/param_fix_report.json"

7.2 汇总待修改项（Think Aloud 说明每一处修改的理由）

7.3 修正教程原文
    python tools/orbital_validator.py "$tutorial_path" --fix

7.4 说明是否需要重走 CLAUDE.md 审查流程
    - 只删除 nbands auto：无需重走审查
    - 轨道文件名替换：建议重走 Step 5（案例审查）
```

### param_fix_report.json 格式

```json
{
  "tutorial": "_workspace/XXX/07_final.md",
  "fixes": [
    {
      "type": "input_param",
      "param": "nbands auto",
      "action": "removed",
      "reason": "ABACUS v3.10.x不支持auto关键字，删除后ABACUS自动计算",
      "location": "教程第XX行代码块",
      "timestamp": "2026-02-27T10:30:00"
    }
  ]
}
```

---

## 不修改的文件

- `tools/test_plugins/elastic_plugin.py`：内部 submit+wait 逻辑在 `--phase prepare` 下不被调用，保留做兜底
- `tools/test_plugins/band_plugin.py`：同上
- `docs/` 内历史文档
- `tools/fix_stru.py`

---

## 历史问题解决状态

| 历史问题 | 解决方案 | 状态 |
|---------|---------|------|
| 轨道文件名错误 | orbital_validator（生成）+ orbital_db（测试） | ✅ 已解决 |
| STRU 文件格式 | fix_stru.py，elastic_plugin 自动调用 | ✅ 已解决 |
| Bohrium 环境变量 | job.json 固定命令 | ✅ 已解决 |
| nbands auto | orbital_validator 扩展（生成）+ Step 2.5 grep（测试） | 🔧 本次解决 |
| 双重任务提交 | --phase prepare 替换 run_full_test | 🔧 本次解决 |
| 串行/并行决策僵化 | Claude 读教程 Think Aloud，完全灵活 | 🔧 本次解决 |
| 测试修复反向写回教程 | Step 7（新） | 🔧 本次解决 |
