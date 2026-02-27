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

## 完整测试流程（7步）

### Step 0: 环境检查与配置

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

**0.2.1 检查Bohrium CLI安装**

```bash
bohr version
```

**预期输出：** `1.1.0` 或其他版本号

**如果未安装：**
```
❌ Bohrium CLI未安装

请在PowerShell中运行以下命令安装：
curl -o install_bohr_windows_curl.bat https://dp-public.oss-cn-beijing.aliyuncs.com/bohrctl/1.0.0/install_bohr_windows_curl.bat && .\install_bohr_windows_curl.bat

安装完成后，重启VSCode，告诉我"我已重启VSCode，继续测试"
```
- 创建状态文件`BOHRIUM_CONFIG_STATUS.md`（记录进度）
- **停止执行，等待用户重启**

---

**0.2.2 配置AccessKey**

```bash
# 检查AccessKey
echo $env:ACCESS_KEY
```

**如果为空：**
```
❌ ACCESS_KEY未设置

请执行以下步骤：

1. 获取AccessKey
   - 访问 https://bohrium.dp.tech/
   - 登录后在"个人中心"→"AccessKey"获取

2. 在PowerShell中配置（请替换your_access_key_here）：
   setx ACCESS_KEY "your_access_key_here"

3. 重启VSCode使环境变量生效

完成后，告诉我"我已重启VSCode，继续测试"
```

- 创建状态文件：
  ```markdown
  # Bohrium配置状态

  ## 当前任务
  测试教程：[教程路径]

  ## 当前进度
  - [x] 检查Python环境
  - [x] 检查测试框架脚本
  - [x] 检查Bohrium CLI安装
  - [x] 配置AccessKey
  - [ ] 验证配置（等待重启）

  ## 下一步
  重启VSCode后，告诉Claude：
  "我已重启VSCode，继续测试"

  Claude将从"验证配置"步骤继续。
  ```
- **停止执行，等待用户重启**

---

**0.2.3 验证配置（用户重启后）**

当用户说"我已重启VSCode"或类似表述时：

1. 读取状态文件（如果存在）
2. 执行验证命令：
```bash
# 验证AccessKey
echo $env:ACCESS_KEY

# 测试API连接
bohr project list
```

**预期输出：**
- AccessKey显示正确的值
- 项目列表显示用户的项目

**如果成功：**
```
✅ Bohrium CLI配置成功！

可用项目：
- My Default Project (ID: 205855)

环境配置完成！
```

**如果失败：**
- 检查AccessKey是否正确
- 确认用户已完全重启VSCode（不是重新加载窗口）
- 检查网络连接
- 提供故障排除建议（见附录B）

---

**0.2.4 配置项目ID（可选）**

如果用户有多个项目，设置默认项目：
```bash
bohr config set project_id 205855
```

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
# 生成时间戳目录
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$tutorial_name = "elastic_tutorial"  # 从教程文件名提取
mkdir test_workspace/${timestamp}_${tutorial_name}
cd test_workspace/${timestamp}_${tutorial_name}
```

**1.2 调用教程解析工具**

```bash
python ../../tools/test_framework_integrated.py --tutorial ../../_workspace/XXX/07_final.md --phase analyze
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
教程标题：弹性常数计算
教程路径：_workspace/20260203_105918_弹性常数计算/07_final.md

[解析] 发现案例：Si（立方晶系）
  计算步骤：
    1. 结构优化（cell-relax）
    2. 弹性常数计算（elastic）

  文件清单：
    - INPUT（已提取）
    - STRU（已提取）
    - KPT（已提取）

  依赖文件：
    - Si_ONCV_PBE-1.0.upf（赝势）
    - Si_gga_7au_100Ry_2s2p1d.orb（轨道）

  预期结果：
    - C₁₁ = 155.5 GPa
    - C₁₂ = 58.3 GPa
    - C₄₄ = 76.2 GPa
    - 体模量 = 90.7 GPa

[保存] 解析结果 → analysis.json
✅ 解析完成
```

**1.3 验证解析结果**

检查`analysis.json`内容：
- 案例信息完整
- 文件内容有效（非空）
- 依赖文件列表清晰

**Think Aloud：**
- 说明教程路径和主题
- 说明发现的案例数量和类型
- 说明提取的文件和参数
- 说明预期结果
- 如发现问题（格式错误、信息缺失），说明具体问题

**完成标志：**
- `analysis.json`已生成
- 输出："✅ Step 1完成，发现N个测试案例"

---

### Step 2: 准备输入文件

**目标：** 根据解析结果准备完整的ABACUS输入文件目录

**执行：**

**2.1 创建案例目录结构**

```bash
# 示例：Si案例
mkdir Si
mkdir Si/01_relax
mkdir Si/02_elastic
```

**2.2 写入ABACUS输入文件**

从`analysis.json`提取内容，写入文件：
```bash
# 写入INPUT
echo "[INPUT内容]" > Si/01_relax/INPUT

# 写入STRU
echo "[STRU内容]" > Si/01_relax/STRU

# 写入KPT
echo "[KPT内容]" > Si/01_relax/KPT
```

**2.3 下载赝势和轨道文件**

使用内置的文件管理器自动下载：
```bash
python ../../tools/test_framework_integrated.py --phase download_files --case Si
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

[检查] Si_ONCV_PBE-1.0.upf
  ✅ 已缓存（tools/pseudopotentials/Si_ONCV_PBE-1.0.upf）
  [复制] → Si/01_relax/Si_ONCV_PBE-1.0.upf

[检查] Si_gga_7au_100Ry_2s2p1d.orb
  ⚠️  文件不存在，应用映射表
  [映射] Si_gga_7au_100Ry_2s2p1d.orb → Si_gga_8au_100Ry_2s2p1d.orb
  [下载] https://abacus.deepmodeling.com/.../Si_gga_8au_100Ry_2s2p1d.orb
  ✅ 下载完成（78 KB）
  [缓存] → tools/orbitals/Si_gga_8au_100Ry_2s2p1d.orb
  [复制] → Si/01_relax/Si_gga_8au_100Ry_2s2p1d.orb

✅ 所有依赖文件已准备
```

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

**2.5 修复STRU文件格式（如需要）**

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

**Think Aloud：**
- 说明创建的目录结构
- 说明文件下载情况（缓存命中/新下载）
- 说明是否应用了文件映射
- 说明是否修复了STRU格式
- 说明任务配置（机型、镜像、核数）

**完成标志：**
- 所有案例目录已创建
- INPUT/STRU/KPT文件已写入
- 赝势和轨道文件已就位
- job.json已生成
- 输出："✅ Step 2完成，输入文件已准备"

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
# 进入案例目录
cd Si/01_relax

# 提交任务
bohr job submit -i job.json -p ./

# 记录Job ID
```

**输出示例：**
```
=== 提交计算任务 ===

[提交] Si-relax（结构优化）
  目录：Si/01_relax
  配置：c16_m32_cpu, 8核
  命令：export ABACUS_PP_PATH=./ && export ABACUS_ORB_PATH=./ && OMP_NUM_THREADS=1 mpirun -np 8 abacus

  ✅ 提交成功
  Job ID: 21984267
  状态：Waiting

[等待] Si-relax完成后再提交Si-elastic...

[记录] 任务信息 → job_ids.json
```

**3.3 记录任务信息**

保存到`job_ids.json`：
```json
{
  "Si-relax": {
    "job_id": "21984267",
    "status": "waiting",
    "submit_time": "2026-02-26 10:30:00",
    "directory": "Si/01_relax",
    "depends_on": null
  },
  "Si-elastic": {
    "job_id": null,
    "status": "pending",
    "submit_time": null,
    "directory": "Si/02_elastic",
    "depends_on": "21984267"
  }
}
```

**Think Aloud：**
- 说明提交的任务数量和类型
- 说明任务依赖关系
- 说明使用的机型和资源
- 说明预计运行时间
- 说明Job ID记录位置

**完成标志：**
- 所有任务已提交或等待依赖完成
- job_ids.json已保存
- 输出："✅ Step 3完成，已提交N个任务"

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
- 查看日志：test_workspace/job_logs/job_21984267/
- 参考故障排除清单（附录B）

是否继续测试其他案例？(Y/n)
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
mkdir -p ../../job_logs

# 下载每个任务的结果
for job_id in $(cat job_ids.json | jq -r '.[] | select(.status == "finished") | .job_id'); do
  echo "[下载] Job $job_id..."
  bohr job download -j $job_id -o ../../job_logs/job_$job_id
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

将下载的结果复制回原计算目录：
```bash
# 复制优化后的结构到弹性计算目录
cp ../../job_logs/job_21984267/OUT.ABACUS/STRU_ION_D Si/02_elastic/STRU
```

**5.3 运行后处理（如需要）**

对于弹性计算，运行后处理提取弹性常数：
```bash
cd Si/02_elastic
python ../../../tools/test_framework_integrated.py --phase postprocess --type elastic
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
=== 结果对比 ===

案例：Si（立方晶系）

弹性常数对比：
┌─────────────┬──────────┬──────────┬──────────┬────────┐
│ 参数        │ 教程预期 │ 实际结果 │ 相对误差 │ 状态   │
├─────────────┼──────────┼──────────┼──────────┼────────┤
│ C₁₁ (GPa)   │ 155.5    │ 155.46   │ 0.03%    │ ✅ PASS│
│ C₁₂ (GPa)   │ 58.3     │ 58.33    │ 0.05%    │ ✅ PASS│
│ C₄₄ (GPa)   │ 76.2     │ 76.18    │ 0.03%    │ ✅ PASS│
└─────────────┴──────────┴──────────┴──────────┴────────┘

弹性模量对比：
┌─────────────┬──────────┬──────────┬──────────┬────────┐
│ 参数        │ 教程预期 │ 实际结果 │ 相对误差 │ 状态   │
├─────────────┼──────────┼──────────┼──────────┼────────┤
│ 体模量(GPa) │ 90.7     │ 90.71    │ 0.01%    │ ✅ PASS│
│剪切模量(GPa)│ 65.1     │ 65.13    │ 0.05%    │ ✅ PASS│
│杨氏模量(GPa)│ 157.7    │ 157.66   │ 0.03%    │ ✅ PASS│
│ 泊松比      │ 0.21     │ 0.210    │ 0.00%    │ ✅ PASS│
└─────────────┴──────────┴──────────┴──────────┴────────┘

✅ 测试通过（7/7）
最大相对误差：0.05%
```

**6.3 处理测试失败**

如果相对误差超出容差：
```
⚠️  部分测试失败

案例：Si

弹性常数对比：
┌─────────────┬──────────┬──────────┬──────────┬────────┐
│ 参数        │ 教程预期 │ 实际结果 │ 相对误差 │ 状态   │
├─────────────┼──────────┼──────────┼──────────┼────────┤
│ C₁₁ (GPa)   │ 155.5    │ 142.3    │ 8.49%    │ ❌ FAIL│
│ C₁₂ (GPa)   │ 58.3     │ 58.33    │ 0.05%    │ ✅ PASS│
└─────────────┴──────────┴──────────┴──────────┴────────┘

❌ 测试失败（1/2 失败）

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
```

**6.4 生成测试报告**

在_workspace/XXX/中
创建`test_report.md`：
```markdown
# 教程测试报告

**生成时间：** 2026-02-26 11:15:32
**测试教程：** _workspace/20260203_105918_弹性常数计算/07_final.md
**测试目录：** test_workspace/20260226_103000_elastic_tutorial/

---

## 测试概要

- ✅ **测试状态：** 通过
- 📊 **测试案例：** 1个（Si）
- ⏱️ **总耗时：** 约15分钟

---

## 案例：Si（立方晶系）

### 计算步骤

1. ✅ 结构优化（Job 21984267）- 21秒
2. ✅ 弹性常数计算（Job 21984270）- 8分钟

### 弹性常数对比

| 参数 | 教程预期 | 实际结果 | 相对误差 | 状态 |
|------|---------|---------|---------|------|
| C₁₁ (GPa) | 155.5 | 155.46 | 0.03% | ✅ PASS |
| C₁₂ (GPa) | 58.3 | 58.33 | 0.05% | ✅ PASS |
| C₄₄ (GPa) | 76.2 | 76.18 | 0.03% | ✅ PASS |

### 弹性模量对比

| 参数 | 教程预期 | 实际结果 | 相对误差 | 状态 |
|------|---------|---------|---------|------|
| 体模量 (GPa) | 90.7 | 90.71 | 0.01% | ✅ PASS |
| 剪切模量 (GPa) | 65.1 | 65.13 | 0.05% | ✅ PASS |
| 杨氏模量 (GPa) | 157.7 | 157.66 | 0.03% | ✅ PASS |
| 泊松比 | 0.21 | 0.210 | 0.00% | ✅ PASS |

### 总结

✅ 所有测试通过（7/7）
- 最大相对误差：0.05%
- 所有参数在5%容差范围内

---

## 结论

✅ **教程内容准确**，计算流程可复现，结果与预期一致。

---

## 测试详情

### 任务信息
- Job 21984267: Si-relax（c16_m32_cpu, 8核）
- Job 21984270: Si-elastic（c16_m32_cpu, 8核）

### 文件清单
- 输入文件：`Si/01_relax/`, `Si/02_elastic/`
- 计算结果：`job_logs/job_21984267/`, `job_logs/job_21984270/`
- 解析数据：`analysis.json`

### 环境信息
- Python: 3.11.0
- Bohrium CLI: 1.1.0
- ABACUS镜像: registry.dp.tech/dptech/abacus:LTSv3.10.1

---

**测试框架版本：** AutoTutorial 3.0
**生成工具：** test_framework_integrated.py
```

**Think Aloud：**
- 说明对比的参数数量和通过率
- 说明最大相对误差
- 如有失败，说明失败原因和建议
- 说明测试报告保存位置

**完成标志：**
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

- 结果对比已完成
- `test_report.md`已生成
- 输出："✅ Step 6完成，测试报告已生成"
- 告诉用户测试报告路径

---

## 工作总结

完成所有步骤后，输出最终总结：

```
╔════════════════════════════════════════════════╗
║     ✅ 教程测试完成                            ║
╚════════════════════════════════════════════════╝

测试教程：
  _workspace/20260203_105918_弹性常数计算/07_final.md

测试结果：
  ✅ 通过（7/7参数）

测试报告：
  test_workspace/20260226_103000_elastic_tutorial/test_report.md

关键发现：
  - 所有弹性常数在5%容差范围内
  - 最大相对误差：0.05%
  - 教程内容准确可复现

总耗时：约15分钟

下一步建议：
  - 查看完整报告：test_report.md
  - 如需调试，查看日志：job_logs/
  - 如需重新测试，保留当前目录，重新运行即可
```

---

## 附录A: Bohrium CLI命令速查

### 项目管理
```bash
# 列出项目
bohr project list

# 设置默认项目
bohr config set project_id 205855
```

### 任务管理
```bash
# 提交任务
bohr job submit -i job.json -p ./

# 查询任务状态
bohr job status -j <job_id>

# 列出所有任务
bohr job list

# 下载任务结果
bohr job download -j <job_id> -o <output_dir>

# 终止任务
bohr job kill -j <job_id>
```

### 资源查询
```bash
# 列出可用机型
bohr machine list

# 列出可用镜像
bohr image list

# 搜索ABACUS镜像
bohr image list | grep abacus
```

### 配置管理
```bash
# 查看当前配置
bohr config list

# 设置AccessKey
bohr config set access_key <your_key>

# 设置项目ID
bohr config set project_id <project_id>
```

---

## 附录B: 故障排除清单

### 问题1：Bohrium CLI未找到

**症状：**
```
'bohr' 不是内部或外部命令
```

**解决方案：**
1. 确认已安装Bohrium CLI
2. 重启VSCode（PATH更新）
3. 在PowerShell中运行`refreshenv`

---

### 问题2：ACCESS_KEY未设置

**症状：**
```
echo $env:ACCESS_KEY
# 输出为空
```

**解决方案：**
1. 在PowerShell中运行：`setx ACCESS_KEY "your_key"`
2. 完全关闭并重新打开VSCode（不是重新加载窗口）
3. 验证：`echo $env:ACCESS_KEY`

---

### 问题3：API连接失败

**症状：**
```
bohr project list
Error: Failed to connect to API
```

**解决方案：**
1. 检查AccessKey是否正确
2. 检查网络连接
3. 尝试访问 https://bohrium.dp.tech/
4. 检查防火墙设置

---

### 问题4：任务提交失败

**症状：**
```
Error: Job submission failed
```

**可能原因：**
1. 项目ID未设置或错误
2. 项目配额不足
3. 机型或镜像不可用
4. job.json格式错误

**解决方案：**
1. 检查项目ID：`bohr config list`
2. 检查项目配额：`bohr project list`
3. 检查机型可用性：`bohr machine list`
4. 验证job.json格式

---

### 问题5：轨道文件下载为HTML

**症状：**
下载的`.orb`文件内容是HTML（404页面）

**解决方案：**
- 自动应用文件映射表（已内置）
- 常见映射：
  - `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`
  - `Ti_gga_8au_100Ry_2s2p2d1f.orb` → `Ti_gga_8au_100Ry_4s2p2d1f.orb`
  - `O_gga_6au_100Ry_2s2p1d.orb` → `O_gga_7au_100Ry_2s2p1d.orb`

---

### 问题6：STRU文件格式错误

**症状：**
```
ERROR: autotest2006 has detected a ABACUS STRU file format error
```

**解决方案：**
- 自动修复工具已内置（`fix_stru.py`）
- 新格式要求：
  - `ATOMIC_SPECIES`只包含3个字段
  - 轨道文件单独放在`NUMERICAL_ORBITAL`块

---

### 问题7：计算结果与预期差异大

**症状：**
相对误差超过5%

**可能原因：**
1. 教程预期结果不准确
2. 计算参数不同（K点、截断能）
3. ABACUS版本差异
4. 赝势或轨道文件不同

**解决方案：**
1. 检查教程中的预期结果来源
2. 对比INPUT文件参数
3. 查看ABACUS版本（教程生成时 vs 当前测试）
4. 检查赝势和轨道文件版本

---

### 问题8：任务长时间Waiting

**症状：**
任务提交后一直处于Waiting状态

**可能原因：**
1. 平台资源紧张
2. 机型不可用
3. 项目优先级低

**解决方案：**
1. 等待（通常5-15分钟）
2. 更换机型（如c8_m16_cpu）
3. 联系平台管理员

---

## 附录C: 文件格式参考

### job.json格式

```json
{
  "job_name": "任务名称",
  "command": "执行命令",
  "log_file": "日志文件路径（支持通配符）",
  "project_id": 205855,
  "platform": "ali",
  "job_type": "container",
  "machine_type": "c16_m32_cpu",
  "image_address": "registry.dp.tech/dptech/abacus:LTSv3.10.1"
}
```

**字段说明：**
- `job_name`: 任务名称（用于识别）
- `command`: 容器内执行的命令
- `log_file`: 日志文件路径（用于判断任务完成）
- `project_id`: Bohrium项目ID
- `platform`: 平台类型（ali/aws）
- `job_type`: 任务类型（container）
- `machine_type`: 机型（c8_m16_cpu, c16_m32_cpu等）
- `image_address`: Docker镜像地址

**常用机型：**
- `c8_m16_cpu`: 8核16GB（适合小计算）
- `c16_m32_cpu`: 16核32GB（推荐）
- `c32_m64_cpu`: 32核64GB（大计算）

**ABACUS镜像：**
- `registry.dp.tech/dptech/abacus:3.10.1`: 稳定版
- `registry.dp.tech/dptech/abacus:LTSv3.10.1`: 长期支持版（推荐）
- `registry.dp.tech/dptech/abacus:latest`: 最新版

---

### analysis.json格式

```json
{
  "tutorial_title": "教程标题",
  "tutorial_path": "教程路径",
  "cases": [
    {
      "name": "Si",
      "type": "elastic",
      "steps": [
        {
          "step": "relax",
          "input_files": {
            "INPUT": "...",
            "STRU": "...",
            "KPT": "..."
          },
          "pseudopotentials": ["Si_ONCV_PBE-1.0.upf"],
          "orbitals": ["Si_gga_8au_100Ry_2s2p1d.orb"]
        }
      ],
      "expected_results": {
        "C11": 155.5,
        "C12": 58.3,
        "C44": 76.2,
        "bulk_modulus": 90.7
      },
      "tolerance": 0.05
    }
  ]
}
```

---

### STRU文件新格式

```
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
./Si_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
1 0 0
0 1 0
0 0 1

ATOMIC_POSITIONS
Cartesian
Si
0.0
2
0.0 0.0 0.0 1 1 1
0.25 0.25 0.25 1 1 1
```

**关键要点：**
- `ATOMIC_SPECIES`只包含3个字段：元素、质量、赝势
- `NUMERICAL_ORBITAL`独立块，每行一个轨道文件
- 轨道文件路径需要`./`前缀

---

## 参考文档

### AutoTutorial 3.0相关
- 教程生成指南：`CLAUDE.md`
- 测试功能总结：`docs/CALCULATION_TESTING_SUMMARY.md`
- 成功报告：`docs/2SUCCESS_REPORT.md`
- 综合解决方案：`docs/1COMPREHENSIVE_SOLUTION.md`

### Bohrium平台
- 官网：https://bohrium.dp.tech/
- CLI文档：https://docs.bohrium.com/docs/bohrctl/about
- ABACUS软件案例：https://bohrium-doc.dp.tech/docs/software/ABACUS/

### ABACUS相关
- 官网：http://abacus.deepmodeling.com/
- 赝势库：http://abacus.deepmodeling.com/upf/
- 轨道文件：http://abacus.deepmodeling.com/orbital/

---

**文档版本：** 1.0
**创建日期：** 2026-02-26
**适用范围：** AutoTutorial 3.0 教程测试
**维护者：** AutoTutorial开发团队

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
