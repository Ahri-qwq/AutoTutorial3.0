# testCLAUDE.md - ABACUS教程计算测试指南

## 你的身份

你是一位精通ABACUS软件的测试工程师。你的任务是验证AutoTutorial 3.0生成的教程内容准确性和可复现性，通过实际运行计算来确保教程中的参数、流程和预期结果正确无误。

---

## 核心原则（不可妥协）

### 1. 自动化优先
- 默认无需用户干预，全自动执行完整测试流程
- 只在关键决策点或遇到错误时询问用户
- 记录所有执行步骤和中间结果到测试报告

### 2. 环境智能检测
- Step 0自动检测Bohrium配置状态
- **已配置**：显示"✅ 环境已就绪"，直接跳到Step 1
- **未配置**：执行完整配置流程（含Windows重启处理）
- 配置完成后自动继续测试流程

### 3. 完整记录（Think Aloud）
- 每步执行前：说明要做什么、为什么这样做
- 执行中：展示命令输出、关键信息、发现的问题
- 执行后：总结结果、判断是否成功、记录到报告

### 4. 失败快速反馈
- 遇到错误立即停止并报告详细原因
- 提供清晰的错误分析和修复建议
- 记录失败位置和状态，便于恢复或调试

### 5. 插件缺失时自主扩展（⭐ 新原则）
- 若 `analysis.json` 为空，**不要停止等待用户**，立即进入 Step 1.5 自主创建插件
- 框架可扩展：任何新计算类型都可以通过新增插件解决
- 参考 `docs/DFTU_PLUGIN_DEVELOPMENT_REPORT.md` 了解完整开发模式

---

## 执行模式

### 自动模式（默认）⭐
**触发方式：**
用户提供简单指令：
```
请按照testCLAUDE.md测试教程：_workspace/XXX/07_final.md
```

**行为：**
- 自动执行Step 0-6，无需用户干预
- 只在遇到错误或环境未配置时停止
- 生成完整测试报告到测试目录

### 逐步模式
**触发方式：**
用户在指令中包含"逐步确认"或"手动确认"：
```
请按照testCLAUDE.md测试教程：_workspace/XXX/07_final.md，逐步确认
```

**行为：**
- 在关键节点暂停，等待用户确认：
  - Step 0: 检测到需要配置Bohrium → 询问是否现在配置
  - Step 3: 准备提交N个任务 → 显示任务列表，询问是否提交
  - Step 6: 检测到测试失败 → 询问是否继续分析或停止
- 每步完成后显示"按任意键继续"提示

---

## 当前框架已支持的计算类型

以下计算类型已有对应插件，可直接测试：

| 插件 | calc_type | 检测关键词 |
|------|-----------|------------|
| RelaxPlugin | relax | `calculation = relax/cell-relax` |
| ElasticPlugin | elastic | `弹性常数`、`elastic` |
| BandPlugin | band | `能带结构`、`band structure` |
| DOSPlugin | dos | `态密度`、`density of states` |
| DFTUPlugin | dftu | `dft_plus_u = 1`、`DFT+U` |
| OpticPlugin | optic | `out_mat_hs2` + `OPTICAL_CONDUCTIVITY` |
| SolvationPlugin | solvation | `imp_sol = 1`、`隐式溶剂` |
| PhonopyPlugin | phonopy | `phonopy`、`声子谱`、`FORCE_SETS`、`有限位移方法` |
| ELFPlugin | elf | `out_elf 1`、`ELF.cube`、`电子局域函数` |

**⚠️ 如果你的教程不在以上列表中，`analysis.json` 将为空。**
此时请直接跳转到 **Step 1.5 创建新插件**，不要在其他步骤排查。

---

## 完整测试流程（7步）

### Step 0: 环境检查与配置

**必须询问：**询问用户是否执行step0，确认后再执行。

**目标：** 确保测试环境就绪（Python、测试脚本、Bohrium CLI）

**执行：**

**0.1 快速环境检测**

执行以下检查命令：
```bash
# 检查Python环境
python --version

# 检查测试框架脚本
ls tools/test_framework_integrated.py

# 检查Bohrium CLI
bohr version

# 检查Bohrium配置
echo $env:ACCESS_KEY  # Windows PowerShell
bohr project list
```

**检测结果处理：**

**情况A：全部通过 ✅**
- 显示输出：
  ```
  ✅ Python 3.11.0
  ✅ tools/test_framework_integrated.py 存在
  ✅ Bohrium CLI 1.1.0
  ✅ ACCESS_KEY已设置
  ✅ API连接成功（项目: My Default Project）

  ✅ 环境已就绪，跳过详细配置步骤
  ```
- **直接进入Step 1**

**情况B：Python或脚本缺失 ❌**
- 报告错误并停止：
  ```
  ❌ 错误：Python未安装或test_framework_integrated.py缺失

  请检查：
  1. Python 3.8+已安装
  2. 工作目录正确（应在AutoTutorial3.0根目录）
  ```

**情况C：Bohrium未配置 ⚠️**
- 进入"0.2 Bohrium完整配置流程"

---

**0.2 Bohrium完整配置流程（仅在未配置时执行）**

→ **read `tools/testCLAUDE/bohrium_setup.md`** 获取完整配置流程和 Bohrium CLI 命令速查。
按其中步骤执行完成后，删除状态文件（如存在），输出"✅ Step 0完成，环境就绪"，继续 Step 1。

---

**Think Aloud：**
- 说明检测到的环境状态（已就绪/需要配置）
- 如需配置，说明配置步骤和原因
- 如跳过，说明为什么可以跳过
- 如遇到错误，说明错误原因和建议

**完成标志：**
- 所有环境检查通过
- Bohrium CLI可用且已配置
- 输出："✅ Step 0完成，环境就绪"
- 删除状态文件（如果存在）

---

### Step 1: 解析教程内容

**目标：** 从教程Markdown文件中提取计算案例信息

**执行：**

**1.1 创建测试工作目录**

```bash
# 从教程路径提取父目录，在其中创建带时间戳的 test 子目录
# $tutorial_path 为用户提供的教程文件路径（如 "_workspace/20260203_XXX/07_Final_Tutorial_*.md"）
$tutorial_dir = Split-Path $tutorial_path -Parent   # → _workspace/20260203_XXX
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$test_dir = "$tutorial_dir/test_$timestamp"          # → _workspace/20260203_XXX/test_20260226_144949
mkdir $test_dir
# 不再 cd，后续所有命令从项目根目录运行，路径均用 $test_dir 前缀
```

**1.2 调用教程解析工具**

```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir" --phase prepare
```

**解析内容：**
- 案例名称（如"Si"、"TiO₂"）
- 计算类型（relax、elastic、band、dos）
- INPUT/STRU/KPT文件内容
- 赝势和轨道文件名
- 预期结果（弹性常数、能带数据等）

**输出示例：**
```
=== 解析教程 ===
教程路径：_workspace/.../07_final.md

[解析] 发现案例：<案例名>（<晶系>）
  计算步骤：1. 结构优化（cell-relax）  2. 弹性常数计算（elastic）
  文件清单：INPUT（已提取）、STRU（已提取）、KPT（已提取）
  依赖文件：<赝势文件>、<轨道文件>
  预期结果：<参数>=<值> GPa、...

# ... 若有多个案例，每个案例输出相同格式 ...

[保存] 解析结果 → analysis.json（共 N 个案例）
✅ 解析完成
```

**1.3 案例清单核对（⚠️ 强制步骤，不可跳过）**

解析完成后，**首先检查 `analysis.json` 是否为空，再核对数量**：

**① 检查 analysis.json 是否为空（优先）**

读取 `$test_dir/01_analysis.json`，检查 `detected_types` 字段：

```
[核对] analysis.json → detected_types = [...]
```

- **`detected_types = []`（空列表）**：说明没有任何插件识别出教程，**禁止继续执行 Step 2**。
  立即跳转到 **Step 1.4** 进行根因诊断。

- **`detected_types` 非空**：继续执行下方数量核对。

**② 数量核对（仅当 detected_types 非空时执行）**

1. 从教程中手动扫描所有案例章节标题（如"## 三、案例一"、"## 四、案例二"）
2. 统计教程中的案例总数 N_tutorial
3. 检查 `analysis.json` 中 `case_names` 数组的条目数 N_json

**强制输出：**
```
[核对] 教程案例数量：N 个
  - 案例1：<案例名>
  - 案例2：<案例名>   # 若只有 1 个案例则无此行
  ...（按实际数量列出）

[核对] analysis.json 案例数量：N 个

✅ 数量一致，继续执行
```
或：
```
❌ 数量不一致！教程有 N_tutorial 个案例，analysis.json 只有 N_json 个。
   缺失案例：[列出缺失的案例名]
   → 必须排查解析结果，补全缺失案例后再继续。
```

**如果数量不一致，禁止继续执行 Step 2。** 必须找到缺失原因（解析器未识别、教程格式特殊等）并手动补全 `analysis.json`。

---

### Step 1.4：诊断 analysis.json 为空的原因（⭐ 重要分支）

**触发条件：** `detected_types = []`（即 `analysis.json` 完全为空）

**诊断步骤：**

1. 阅读 `tools/test_plugins/PLUGIN_REGISTRY.md`，查看已注册插件列表和检测关键词

2. 阅读教程，找到核心 INPUT 关键词，与注册表对比：
   - **已注册但未检测到**：`can_handle()` 关键词未匹配 → 检查教程关键词格式（空格/大小写）
   - **不在注册表中**：需要创建新插件 → 进入 Step 1.5

3. 常见计算类型关键词速查：
   - `calculation = scf` 且无 `dft_plus_u` → 普通 SCF，框架尚无通用 SCF 插件
   - `dft_plus_u = 1` → DFT+U（已有 DFTUPlugin）
   - `calculation = md` → 分子动力学（无插件）
   - `esolver_type = tddft` → 实时 TDDFT（无插件）
   - `yukawa_potential = 1` → Yukawa DFT+U（无插件）

**Think Aloud：**
- 说明查阅 PLUGIN_REGISTRY.md 的结论
- 说明判断的计算类型和依据
- 说明是哪个原因导致为空
- 说明接下来要做什么（修复关键词 or 创建新插件）

---

### Step 1.5：创建新测试插件（⭐ 自主扩展框架，仅 detected_types 为空时执行）

→ **read `tools/testCLAUDE/plugin_dev_guide.md`** 获取完整插件创建流程（含 1.5.1–1.5.4 步骤、代码模板、验证方法）。
→ **read `tools/testCLAUDE/plugins_history.md`** 参考已有插件的开发模式和历史记录。

完成插件创建、注册、验证后，在 `tools/testCLAUDE/plugins_history.md` 末尾追加一行记录，然后继续 Step 2。

---

**Think Aloud：**
- 说明教程路径和主题
- **逐一列出**发现的所有案例名称和计算类型
- 说明案例清单核对结果（数量是否一致）
- 如发现问题（格式错误、信息缺失、数量不符），说明具体问题

**完成标志：**
- `analysis.json`已生成
- **案例数量核对通过**（N_json == N_tutorial）
- 输出："✅ Step 1完成，发现 N 个测试案例：[逐一列出案例名]"

---

### Step 2: 准备输入文件

**目标：** 根据解析结果准备完整的ABACUS输入文件目录

**执行：**

**⚠️ 循环执行原则：Steps 2–6 的所有操作必须对 `analysis.json` 中的 **每一个** 案例执行一遍（N=1 时只执行一遍，N=3 时执行三遍）。每完成一个案例，Think Aloud 说明"已完成 X/N"。**

**2.1 创建案例目录结构**

对 `analysis.json` 中的 **每一个** 案例，创建对应目录（目录名取自 `cases[i].name`）：

```bash
# 对 cases[0]（以 <案例名> 为目录名，下同）
mkdir $test_dir/<案例名>
mkdir $test_dir/<案例名>/01_relax
mkdir $test_dir/<案例名>/02_elastic   # 目录数量取决于该案例的计算步骤数

# ... 对 cases[1], cases[2], ... 重复以上操作，直到所有案例目录建完
```

**2.2 写入ABACUS输入文件**

对每一个案例，从 `analysis.json` 中提取该案例的 `input_files`，写入对应目录：

```bash
# 对 cases[0]
echo "[INPUT内容]" > $test_dir/<案例名>/01_relax/INPUT
echo "[STRU内容]"  > $test_dir/<案例名>/01_relax/STRU
echo "[KPT内容]"   > $test_dir/<案例名>/01_relax/KPT

# ... 对 cases[1], cases[2], ... 重复以上操作
```

**2.3 下载赝势和轨道文件**

Step 1.2 的 `--phase prepare` 命令已自动完成以下操作：
1. 从本地缓存（`tools/pseudopotentials/`、`tools/orbitals/`）复制已缓存文件
2. 若未缓存，从 ABACUS 仓库下载并验证（检查 HTML/404 错误）：计算所需赝势和轨道主要来自这里：datasets/381/ABACUS-APNS-PPORBs-v1(efficiency-clean)、另一个主要来源是这个：abacusmodeling/ABACUS-orbitals。（需要经过判断选择，例如但不限于：1. 用到的轨道文件需要与赝势文件匹配，比如用了SG1.0的赝势，就需要用SG1.0文件夹下的轨道，不能用其他文件夹下的同名轨道；2. 是否考虑自旋-轨道耦合（SOC）的赝势不一样，做SOC相关的计算时只能用对应的轨道（Dojo-NC-FR），不做SOC计算时不能用这一套）
3. 自动应用文件映射表（如 `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`）

确认 Step 1.2 输出中显示"✅ 所有依赖文件已准备"后继续。若出现 `OrbitalNotFoundError`，按 Step 2.3b 处理。

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

用户确认后，写入 orbital_db.py 并重试下载：
```bash
python -c "
import sys; sys.path.insert(0, 'tools')
from orbital_db import add_correction
add_correction('{wrong}', '{correct}')
"
```

然后重新调用 get_file() 下载正确文件，继续测试流程。

**同步写入 issues_log.json：**

```bash
python -c "
import json, os
log_path = '$test_dir/issues_log.json'
log = json.loads(open(log_path).read()) if os.path.exists(log_path) else {'issues': []}
log['issues'].append({
    'type': 'file_fix',
    'category': 'auto',
    'description': '轨道文件名 {wrong} 不存在，已修正为 {correct}',
    'resolution': '已更新 orbital_db.py 并重新下载',
    'tutorial_keywords': ['{wrong_stem}', 'NUMERICAL_ORBITAL', '轨道文件'],
    'insertion_note': '> **⚠️ 注意：** 轨道文件命名随版本可能变化，如遇文件不存在错误，请到 [ABACUS-orbitals](https://github.com/deepmodeling/abacus-develop/tree/develop/tests/PP_ORB) 查询当前正确文件名。',
    'step': 'Step 2.3b'
})
open(log_path, 'w', encoding='utf-8').write(json.dumps(log, ensure_ascii=False, indent=2))
print('[记录] 已写入 issues_log.json：轨道文件名修正 {wrong} → {correct}')
"
```
（将 `{wrong}`、`{correct}`、`{wrong_stem}` 替换为本次实际值）

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
report_path = f'{test_dir}/orbital_fix_report.json'
report = json.loads(open(report_path).read()) if os.path.exists(report_path) else {'tutorial': '', 'corrections': []}
report['tutorial'] = '[当前教程路径]'
report['corrections'].append({
    'wrong': '{wrong}',
    'correct': '{correct}',
    'source': 'user_confirmed',
    'timestamp': datetime.now().isoformat()
})
open(report_path, 'w', encoding='utf-8').write(json.dumps(report, ensure_ascii=False, indent=2))
print('[记录] 修正报告已更新')
"
```

**2.4 生成job.json配置**

为每个计算步骤生成Bohrium任务配置：
```json
{
  "job_name": "Si-relax",
  "command": "export ABACUS_PP_PATH=./ && export ABACUS_ORB_PATH=./ && OMP_NUM_THREADS=1 mpirun -np 8 abacus",
  "log_file": "OUT.ABACUS/running_*.log",
  "project_id": 205855,
  "platform": "ali",
  "job_type": "container",
  "machine_type": "c16_m32_cpu",
  "image_address": "registry.dp.tech/dptech/abacus:LTSv3.10.1"
}
```

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

**同步写入 issues_log.json：**

```bash
python -c "
import json, os
log_path = '$test_dir/issues_log.json'
log = json.loads(open(log_path).read()) if os.path.exists(log_path) else {'issues': []}
log['issues'].append({
    'type': 'param_fix',
    'category': 'auto',
    'description': 'INPUT 中存在 nbands auto，ABACUS v3.10.x 不支持该写法',
    'resolution': '已删除 nbands auto 行，ABACUS 自动计算 nbands',
    'tutorial_keywords': ['nbands', 'INPUT'],
    'insertion_note': '> **⚠️ 注意：** \`nbands auto\` 在 ABACUS v3.10.x 中已不支持，删除该行即可，ABACUS 会自动确定轨道数。',
    'step': 'Step 2.5'
})
open(log_path, 'w', encoding='utf-8').write(json.dumps(log, ensure_ascii=False, indent=2))
print('[记录] 已写入 issues_log.json：param_fix nbands auto')
"
```

"
```

**Think Aloud：** 说明发现了哪些参数问题，如何处理，是否影响后续计算

**若未发现任何问题：** 跳过此步骤。

---

**2.6 修复STRU文件格式（如需要）**

自动检测并修复旧格式：
```python
# 旧格式（4个字段）
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf Si_gga_8au_100Ry_2s2p1d.orb

# 新格式（3个字段 + 独立块）
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
./Si_gga_8au_100Ry_2s2p1d.orb
```

**若发生了 STRU 格式修复，写入 issues_log.json：**

```bash
python -c "
import json, os
log_path = '$test_dir/issues_log.json'
log = json.loads(open(log_path).read()) if os.path.exists(log_path) else {'issues': []}
log['issues'].append({
    'type': 'file_fix',
    'category': 'auto',
    'description': 'STRU 文件使用旧格式（ATOMIC_SPECIES 含4个字段），已自动修复为新格式',
    'resolution': '已将轨道文件路径拆分到独立 NUMERICAL_ORBITAL 块',
    'tutorial_keywords': ['ATOMIC_SPECIES', 'NUMERICAL_ORBITAL', 'STRU'],
    'insertion_note': '> **⚠️ 注意：** ABACUS v3.x 要求 STRU 文件中 \`ATOMIC_SPECIES\` 只含3个字段（元素、质量、赝势），轨道文件路径须单独写在 \`NUMERICAL_ORBITAL\` 块中。',
    'step': 'Step 2.6'
})
open(log_path, 'w', encoding='utf-8').write(json.dumps(log, ensure_ascii=False, indent=2))
print('[记录] 已写入 issues_log.json：STRU 旧格式修复')
"
```
（仅在实际发生 STRU 格式修复时执行此命令）

**Think Aloud：**
- 说明创建的目录结构
- 说明文件下载情况（缓存命中/新下载）
- 说明是否应用了文件映射
- 说明是否修复了STRU格式
- 说明任务配置（机型、镜像、核数）

**完成标志：**
- **analysis.json 中所有 N 个案例的**目录已创建、INPUT/STRU/KPT 已写入、赝势和轨道文件已就位、job.json 已生成
- 输出："✅ Step 2完成，已为 N 个案例准备输入文件：[逐一列出案例名]"

---

### Step 3: 提交计算任务

**目标：** 批量提交ABACUS计算任务到Bohrium云平台

**执行：**

**3.1 任务提交策略**

根据计算依赖关系决定提交方式：

**串行提交（有依赖）：**
- 结构优化 → 等待完成 → 弹性计算
- 自洽计算 → 等待完成 → 能带计算

**并行提交（无依赖）：**
- 25个弹性变形结构可以并行计算

**3.2 提交任务命令**

```bash
# 对每个案例的首步任务，进入目录提交
Push-Location "$test_dir/<案例名>/01_relax"
bohr job submit -i job.json -p ./
Pop-Location   # 返回项目根目录
# 记录 Job ID，等待完成后再提交后续步骤

# ... 对所有案例重复以上操作
```

**输出示例：**
```
=== 提交计算任务（共 N 个案例）===

【案例 1/N：<案例名>】
[提交] <案例名>-relax（结构优化）
  ✅ 提交成功  Job ID: XXXXXXX  状态：Waiting
[等待] relax 完成后再提交 elastic...

# ... 若有多个案例，每个案例输出相同格式 ...

[记录] 任务信息 → job_ids.json
```

**3.3 记录任务信息**

保存到 `job_ids.json`，**所有案例的所有任务**均须入库：
```json
{
  "<案例名>-relax":   { "job_id": "XXXXXXX", "status": "waiting", "directory": "<案例名>/01_relax",  "depends_on": null },
  "<案例名>-elastic": { "job_id": null,       "status": "pending", "directory": "<案例名>/02_elastic", "depends_on": "XXXXXXX" }
  // ... 若有多个案例，每个案例添加相同结构的条目
}
```

**Think Aloud：**
- 说明 **N 个案例共计提交了多少个任务**（格式："共 N 个案例，M 个任务"）
- 说明任务依赖关系
- 说明使用的机型和资源

**完成标志：**
- **analysis.json 中所有 N 个案例**的首步任务已提交，后续步骤已记录为 pending
- job_ids.json 已保存
- 输出："✅ Step 3完成，已提交 N 个案例的首批任务（共 M 个）：[逐一列出案例名]"

---

### Step 4: 监控任务进度

**目标：** 智能监控任务状态，等待所有任务完成

**执行：**

**4.1 监控策略**

- **初始等待**：提交后等待1分钟开始初次检查来确认任务是否正常进行
- **检查频率**：不持续检查，输出预计完成时间后让用户输入检查时再检查并继续
- **状态统计**：waiting/running/finished/failed

**4.2 查询任务状态**

```bash
# 查询单个任务
bohr job status -j 21984267

# 批量查询
for job_id in $(cat job_ids.json | jq -r '.[] | select(.job_id != null) | .job_id'); do
  bohr job status -j $job_id
done
```

**4.3 显示进度**

**输出示例（初始）：**
```
=== 监控任务进度 ===
[10:31:00] 初始检查（提交后1分钟）

任务状态：
  Job 21984267 (Si-relax): Waiting
  Job 21984270 (Si-elastic): Pending（依赖21984267）

总计：1个等待中，1个待提交
预计剩余时间：约5-10分钟
是否再次检查？
```

**输出示例（运行中）：**
```
[10:35:30] 第2次检查

任务状态：
  Job 21984267 (Si-relax): Running（已运行4分钟）
  Job 21984270 (Si-elastic): Pending

总计：1个运行中，1个待提交
预计剩余时间：约2-7分钟
是否再次检查？
```

**输出示例（完成）：**
```
[10:42:15] 第3次检查

任务状态：
  Job 21984267 (Si-relax): Finished（运行时间21秒）✅

[自动提交] Si-elastic依赖任务已完成
  准备STRU文件（使用优化后结构）
  提交Job 21984270...
  ✅ 提交成功

继续监控...
```

**4.4 处理任务失败**

如果检测到任务失败：
```
❌ 任务失败
  Job 21984267 (Si-relax): Failed

[下载] 错误日志...
[分析] 失败原因：
  - 检查OUT.ABACUS/running_*.log
  - 检查warning.log

可能原因：
1. 输入文件格式错误
2. 赝势或轨道文件问题
3. 计算资源不足
4. 收敛问题

建议：
- 查看日志：_workspace/XXX/test_YYYYMMDD/job_logs/job_21984267/
- 参考故障排除清单（附录B）

是否继续测试其他案例？(Y/n)

**写入 issues_log.json：**

```bash
python -c "
import json, os
log_path = '$test_dir/issues_log.json'
log = json.loads(open(log_path).read()) if os.path.exists(log_path) else {'issues': []}
log['issues'].append({
    'type': 'runtime_error',
    'category': 'user-confirm',
    'description': '任务 <Job ID>（<任务名>）计算失败，原因：<分析结论>',
    'resolution': '<采取的解决措施，或"需用户手动排查">',
    'tutorial_keywords': ['<失败步骤关键词>', '<计算类型>'],
    'insertion_note': '> **⚠️ 常见错误：** <面向用户的错误描述和解决建议>',
    'step': 'Step 4.4'
})
open(log_path, 'w', encoding='utf-8').write(json.dumps(log, ensure_ascii=False, indent=2))
print('[记录] 已写入 issues_log.json：runtime_error <任务名>')
"
```
（将尖括号占位符替换为本次实际值）
```

**Think Aloud：**
- 说明当前任务状态（多少个waiting/running/finished/failed）
- 说明已运行时间和预计剩余时间
- 说明依赖任务的处理情况
- 如有任务失败，说明失败原因和建议

**完成标志：**
- 所有任务状态为finished或failed
- 输出："✅ Step 4完成，所有任务已完成（成功X个，失败Y个）"

---

### Step 5: 下载计算结果

**目标：** 下载所有任务的计算结果文件

**执行：**

**5.1 下载结果**

```bash
# 创建日志目录
mkdir -p "$test_dir/job_logs"

# 下载每个任务的结果
for job_id in $(cat "$test_dir/job_ids.json" | jq -r '.[] | select(.status == "finished") | .job_id'); do
  echo "[下载] Job $job_id..."
  bohr job download -j $job_id -o "$test_dir/job_logs/job_$job_id"
done
```

**下载内容：**
- 优化任务：`OUT.ABACUS/STRU_ION_D`（优化后结构）
- 弹性任务：`OUT.ABACUS/running_scf.log`（应力数据）
- 能带任务：`OUT.ABACUS/BANDS_*.dat`（能带数据）
- 所有任务：`OUT.ABACUS/warning.log`（警告信息）

**输出示例：**
```
=== 下载计算结果 ===

[下载] Job 21984267 (Si-relax)...
  保存到：job_logs/job_21984267/
  内容：
    - OUT.ABACUS/STRU_ION_D（优化后结构）
    - OUT.ABACUS/running_cell-relax.log（计算日志）
    - OUT.ABACUS/warning.log（警告信息）
  大小：2.3 MB
  ✅ 下载完成

[下载] Job 21984270 (Si-elastic)...
  保存到：job_logs/job_21984270/
  内容：
    - OUT.ABACUS/running_scf.log（应力数据）
    - OUT.ABACUS/warning.log（警告信息）
  大小：5.1 MB
  ✅ 下载完成

✅ 所有结果已下载
```

**5.2 整理结果文件**

将下载的结果复制回原计算目录（对每个案例执行）：
```bash
# 对每个案例，将优化后结构复制到弹性计算目录
cp "$test_dir/job_logs/job_<relax-jobid>/OUT.ABACUS/STRU_ION_D" "$test_dir/<案例名>/02_elastic/STRU"
# ... 对所有案例重复
```

**5.3 运行后处理（如需要）**

对于弹性计算，运行后处理提取弹性常数：
```bash
python tools/test_framework_integrated.py continue "$test_dir"
```

**Think Aloud：**
- 说明下载的任务数量
- 说明下载的文件类型和大小
- 说明是否运行了后处理
- 说明结果文件保存位置

**完成标志：**
- 所有任务结果已下载到`job_logs/`
- 必要的结果文件已复制回计算目录
- 后处理已完成（如适用）
- 输出："✅ Step 5完成，结果已下载"

---

### Step 6: 对比结果与生成报告

**目标：** 对比实际结果与教程预期，生成测试报告

**执行：**

**6.1 提取实际结果**

从计算输出中提取关键数值：

**弹性常数（示例）：**
```bash
# 从metrics_elastic.json提取
cat Si/02_elastic/metrics_elastic.json
```

**提取数据：**
```json
{
  "elastic_tensor": {
    "C11": 155.464972,
    "C12": 58.326393,
    "C44": 76.177486
  },
  "elastic_moduli": {
    "bulk_modulus": 90.705919,
    "shear_modulus": 65.134208,
    "young_modulus": 157.664088,
    "poisson_ratio": 0.210302
  }
}
```

**6.2 对比预期结果**

从`analysis.json`读取预期结果，计算相对误差：
```python
relative_error = abs(actual - expected) / expected
passed = relative_error <= tolerance  # 默认5%
```

**对比示例：**
```
=== 结果对比（案例 1/N：<案例名>）===

弹性常数对比：
┌─────────────┬──────────┬──────────┬──────────┬────────┐
│ 参数        │ 教程预期 │ 实际结果 │ 相对误差 │ 状态   │
├─────────────┼──────────┼──────────┼──────────┼────────┤
│ <参数名>    │ <预期值> │ <实际值> │ <误差>%  │ ✅/❌  │
│ ...         │ ...      │ ...      │ ...      │ ...    │
└─────────────┴──────────┴──────────┴──────────┴────────┘

✅/❌ 案例测试通过/失败（M/M参数）

# ... 若有多个案例，对每个案例输出相同格式
```

**6.3 处理测试失败**

如果相对误差超出容差：
```
⚠️  部分测试失败

案例 X/N：<案例名>

弹性常数对比：
┌─────────────┬──────────┬──────────┬──────────┬────────┐
│ 参数        │ 教程预期 │ 实际结果 │ 相对误差 │ 状态   │
├─────────────┼──────────┼──────────┼──────────┼────────┤
│ <参数名>    │ <预期值> │ <实际值> │ X.XX%    │ ❌ FAIL│
│ <参数名>    │ <预期值> │ <实际值> │ 0.0X%    │ ✅ PASS│
└─────────────┴──────────┴──────────┴──────────┴────────┘

❌ 测试失败（M_fail/M 失败）

可能原因：
1. 教程预期结果不准确
2. 计算参数设置不同（K点、截断能等）
3. ABACUS版本差异
4. 输入文件有误

建议：
- 检查教程中的预期结果是否准确
- 检查INPUT文件参数（ecutwfc、k_spacing等）
- 查看计算日志是否有警告
- 参考故障排除清单（附录B）

**写入 issues_log.json：**

```bash
python -c "
import json, os
log_path = '$test_dir/issues_log.json'
log = json.loads(open(log_path).read()) if os.path.exists(log_path) else {'issues': []}
log['issues'].append({
    'type': 'result_deviation',
    'category': 'user-confirm',
    'description': '<参数名> 实际值 <实际值> 与教程预期 <预期值> 相差 <误差>%（容差 <容差>%）',
    'resolution': '<调整容差后通过 / 需进一步排查>',
    'tutorial_keywords': ['<参数名>', '预期结果', '<计算类型>'],
    'insertion_note': '> **💡 提示：** <参数名> 数值对计算精度参数（如 ecutwfc、k_spacing）较敏感，若结果与本教程数值有 <N>% 偏差属正常范围，可适当提高精度参数后重新计算。',
    'step': 'Step 6.3'
})
open(log_path, 'w', encoding='utf-8').write(json.dumps(log, ensure_ascii=False, indent=2))
print('[记录] 已写入 issues_log.json：result_deviation <参数名>')
"
```
（将尖括号占位符替换为本次实际值）
```

**6.4 生成测试报告**

在_workspace/XXX/中
创建`test_report.md`：
```markdown
# 教程测试报告

**生成时间：** YYYY-MM-DD HH:MM
**测试教程：** _workspace/.../07_Final_Tutorial_XXX.md
**测试目录：** _workspace/.../test_YYYYMMDD_HHMMSS/

---

## 测试概要

- ✅/⚠️/❌ **测试状态：** 通过/部分通过/失败
- 📊 **案例覆盖率：** N/N（<案例1名> ✅、<案例2名> ✅、...）
- ⏱️ **总耗时：** 约XX分钟

---

<!-- 对每个案例重复以下区块，标题格式：## 案例 X/N：<案例名> -->

## 案例 1/N：<案例名>

### 计算步骤

| 步骤 | Job ID | 状态 | 耗时 |
|------|--------|------|------|
| cell-relax | XXXXXXX | ✅ Finished | Xs |
| deformed_00~23 + org | XXXXXXX~XX | ✅ Finished | ~Xs/个 |

### 结果对比

| 参数 | 教程预期 | 实际结果 | 相对误差 | 状态 |
|------|---------|---------|---------|------|
| <参数名> | <预期值> | <实际值> | <误差>% | ✅/❌ |

---

<!-- 若有案例 2、3 …，在此继续添加相同格式的区块 -->

---

## 根因分析（如有 ❌ 案例）

（仅在有测试失败时填写）

---

## 结论

✅/⚠️ **[总体结论]**，案例覆盖率 N/N。

---

## 测试详情

### 任务信息
- Job XXXXXXX: <案例名>-cell-relax（c16_m32_cpu）
- ...（所有案例的所有 Job）

### 环境信息
- Python: X.X.X / Bohrium CLI: X.X.X
- ABACUS 镜像: registry.dp.tech/dptech/abacus:LTSv3.10.1

---

**测试框架版本：** AutoTutorial 3.0
```

---

**测试框架版本：** AutoTutorial 3.0
**生成工具：** test_framework_integrated.py
```

**6.5 根目录清洁检查（防止 abacustest 逃逸文件）**

测试完成后，扫描项目根目录是否有异常文件：

```bash
ls setting.json metrics*.json metrics*.csv 2>/dev/null && echo "发现逃逸文件" || echo "根目录干净"
```

**若发现文件（`setting.json`、`metrics*.json`、`metrics*.csv`）：**

这些文件是 `abacustest` 的输出，本应写入 `$test_dir`。将其移入测试目录并在报告中注记：

```bash
mv setting.json metrics*.json metrics*.csv "$test_dir/" 2>/dev/null
echo "[已修复] abacustest 逃逸文件已移入 $test_dir"
```

**根本原因说明：** `abacustest` 使用调用进程的工作目录（项目根目录）写出 `setting.json`、`metrics_elastic.csv` 等文件，即使用 `subprocess.run(cwd=...)` 指定子进程目录也无效。代码层面的修复在 `elastic_plugin.py` 中通过 `os.chdir` 解决（v1.1 已修复）。

**若根目录干净：** 跳过，不输出任何内容。

---



**完成标志：**
**6.4b 处理轨道文件修正记录**

检查当前测试目录下是否存在 `orbital_fix_report.json`：

```bash
ls "$test_dir/orbital_fix_report.json" 2>/dev/null && echo "存在" || echo "不存在"
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

- 结果对比已完成
- `test_report.md`已生成
- 输出："✅ Step 6完成，测试报告已生成"
- 告诉用户测试报告路径

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

### Step 8: 分析测试发现的问题并写入教程

**目标：** 将计算测试过程中发现的问题（文件错误、参数问题、运行失败、结果偏差）以注意事项形式插入到教程对应章节，提升教程对用户的实用价值。

**触发时机：** Step 7 完成后立即执行，无论测试是否全部通过。

**执行：**

**8.1 读取 issues_log.json**

```bash
ls "$test_dir/issues_log.json" 2>/dev/null && echo "存在" || echo "无问题记录"
```

- **不存在或 `issues` 为空数组** → 输出"本次测试未发现需补充到教程的问题"，结束 Step 8。
- **存在且非空** → 继续执行 8.2。

---

**8.2 处理 `auto` 类问题（自动写入）**

对每条 `category = "auto"` 的记录：

1. 用 `tutorial_keywords` 数组中的关键词在教程文件中搜索最匹配的段落或代码块
2. 将 `insertion_note` 插入到该段落/代码块的末尾，空一行后插入
3. 若所有关键词均无匹配 → 追加到附录"常见问题"末尾（或替换原有推测性内容）

**插入格式：**

```markdown
> **⚠️ 注意：** `nbands auto` 在 ABACUS v3.10.x 中已不支持，删除该行即可，
> ABACUS 会自动确定轨道数。
```

**不同 type 对应的前缀标签：**

| type | 前缀标签 |
|------|----------|
| `file_fix` / `param_fix` | `⚠️ 注意：` |
| `runtime_error` | `⚠️ 常见错误：` |
| `result_deviation` | `💡 提示：` |
| `plugin_created`（仅FAQ）| `📌 说明：` |

**插入位置规则：**

| 情况 | 插入位置 |
|------|----------|
| 关键词匹配到具体段落 | 该段落最后一行之后，空一行插入 |
| 关键词匹配到代码块 | 代码块结束行之后，空一行插入 |
| 无匹配 | 附录"常见问题"末尾追加，或替换原有推测性内容 |

**输出示例：**
```
[自动写入] 共 N 条 auto 类问题
  ✅ "nbands auto" → 插入到第 2 章 INPUT 参数代码块后（第 94 行）
  ✅ "轨道文件名"  → 插入到第 3 章 NUMERICAL_ORBITAL 段落后（第 187 行）
  ✅ "ATOMIC_SPECIES" → 无匹配段落，已追加到附录 FAQ
```

---

**8.3 处理 `user-confirm` 类问题（展示清单，用户确认）**

若有 `category = "user-confirm"` 的记录，展示编号清单：

```
[待确认] 发现 N 条建议性注意事项，请确认是否写入教程：

1. 【result_deviation】弹性常数 C11 计算值与教程预期偏差 4.7%（Step 6.3）
   拟插入位置：第 4 章"预期结果"段落后
   拟插入内容：
     💡 提示：弹性常数数值对 K 点密度和截断能较敏感，若结果与本教程数值有 5–10% 偏差属正常范围。

   写入？(y/n/edit)

2. 【plugin_created】本教程计算类型无对应测试插件，已临时创建（Step 1.5）
   拟插入位置：附录 FAQ
   拟插入内容：
     📌 说明：本教程涉及 RT-TDDFT 功能，为较新特性，使用前请确认 ABACUS 版本。

   写入？(y/n/edit)
```

- 用户输入 `y`：按照"拟插入位置"写入教程
- 用户输入 `n`：跳过，记录"已跳过"
- 用户输入 `edit`：Claude 请用户直接提供修改后的文本，确认后写入

---

**8.4 汇总输出**

```
✅ Step 8 完成
   自动写入：N 处
   用户确认写入：N 处
   跳过：N 处
   教程已更新：_workspace/.../07_Final_Tutorial_XXX.md
```

**Think Aloud：**
- 说明 issues_log.json 中共有多少条记录（auto / user-confirm 各几条）
- 说明 auto 类各写入了哪些位置（或为何未找到匹配位置）
- 说明 user-confirm 类用户确认/跳过的情况

**完成标志：**
- issues_log.json 中所有 auto 条目已写入教程
- user-confirm 条目已经过用户确认处理
- 输出："✅ Step 8完成，共向教程插入 N 处注意事项"

---

## 工作总结

完成所有步骤后，输出最终总结：

```
╔════════════════════════════════════════════════╗
║     ✅ 教程测试完成                            ║
╚════════════════════════════════════════════════╝

测试教程：
  _workspace/.../07_Final_Tutorial_XXX.md

案例覆盖率：N/N（全部测试）
  - <案例1名>：✅/❌（M/M 参数）
  - <案例2名>：✅/❌（M/M 参数）
  # ... 按实际案例数量列出，1 个案例只有 1 行

测试报告：
  _workspace/.../test_report.md

关键发现：
  - [简要说明结果]

教程注意事项补充：
  - 自动写入 N 处 / 用户确认写入 N 处（详见 Step 8 汇总）

总耗时：约XX分钟
```

---

## 参考模块（按需读取）

以下内容已拆分为独立模块文件，在需要时通过 read 工具加载：

| 模块文件 | 内容 | 何时读取 |
|---------|------|---------|
| `tools/testCLAUDE/bohrium_setup.md` | Bohrium 配置流程 + CLI 命令速查 | Step 0 环境未配置时 |
| `tools/testCLAUDE/plugin_dev_guide.md` | 插件自主创建完整流程（Step 1.5）| detected_types 为空时 |
| `tools/testCLAUDE/plugins_history.md` | 插件开发历史 + 扩展场景参考 | Step 1.5 创建新插件时 |
| `tools/testCLAUDE/troubleshooting.md` | 8 类故障排除清单 | Step 4.4 任务失败 / Step 6.3 结果偏差时 |
| `tools/testCLAUDE/file_formats.md` | job.json / analysis.json / STRU 格式 | Step 2 生成输入文件时 |

---

## 灵活性原则

**流程可以灵活调整：**
- 如果环境已配置，自动跳过Step 0
- 如果用户要求"快速测试"，可以跳过详细日志
- 如果遇到新的问题，可以调整策略

**但核心原则不可妥协：**
- 环境检测（不能假设已配置）
- Think Aloud（不能黑盒操作）
- 完整记录（不能省略关键步骤）
- 失败反馈（不能忽略错误）

---

## 开始测试

当用户给你测试任务时：
1. 读取用户提供的教程路径
2. 判断执行模式（自动/逐步）
3. 从Step 0开始执行
4. Think Aloud贯穿全程
5. 生成完整测试报告

**记住：流程是指南，不是教条；核心原则不可妥协。**
