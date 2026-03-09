# testCLAUDE.md 模块化拆分实施计划

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 将 testCLAUDE.md 从 1934 行拆分为主文件（~1280 行）+ 5 个按需读取模块，减轻测试执行时的上下文压力。

**Architecture:** 从 testCLAUDE.md 中剪切 Step 0.2（Bohrium 配置）、Step 1.5（插件创建）及四个附录，分别存为独立模块文件。主文件中的原有位置替换为 3–5 行"桩代码"，明确告诉 Claude 在什么条件下 read 哪个模块。同时修改 CLAUDE.md Step 8，增加 compact 暂停点。

**Tech Stack:** Markdown 文件编辑（Read/Edit/Write 工具），无 Python 改动。

---

## 重要前置约定

1. **操作顺序**：先创建模块文件，再编辑 testCLAUDE.md 删除对应内容，最后提交。
2. **Edit 操作**：每次 Edit 前必须先 Read 确认 old_string 与文件完全一致。本计划中的行号仅供定位参考，以实际文件内容为准。
3. **验证方法**：每个任务完成后，通过行数变化（`wc -l` 或 Bash `(Get-Content file).Count`）验证编辑是否生效。

---

### Task 1: 创建目录和 bohrium_setup.md

**Files:**
- Create: `tools/testCLAUDE/bohrium_setup.md`（内容来自 testCLAUDE.md 两处）

**Step 1: 确认内容边界**

Read `testCLAUDE.md` 行 146–264（Step 0.2 完整配置流程）：
```bash
# 用 Bash 验证行号
(Get-Content testCLAUDE.md)[145..263] | Select-Object -First 5
(Get-Content testCLAUDE.md)[145..263] | Select-Object -Last 5
```
预期首行：`**0.2 Bohrium完整配置流程（仅在未配置时执行）**`
预期末行：`bohr config set project_id 205855`（0.2.4 结尾）

Read `testCLAUDE.md` 行 1563–1614（附录A）：
预期首行：`## 附录A: Bohrium CLI命令速查`
预期末行：`bohr config set project_id <project_id>`

**Step 2: 创建 tools/testCLAUDE/ 目录**

```bash
mkdir tools/testCLAUDE
```

**Step 3: 写入 bohrium_setup.md**

Write `tools/testCLAUDE/bohrium_setup.md`，内容如下（两部分拼接）：

```markdown
# Bohrium 配置流程与 CLI 命令速查

> 本文件由 testCLAUDE.md 按需 read，仅在 Step 0 检测到 Bohrium 未配置时加载。

---

## Step 0.2 Bohrium完整配置流程

**0.2.1 检查Bohrium CLI安装**

\`\`\`bash
bohr version
\`\`\`

**预期输出：** `1.1.0` 或其他版本号

**如果未安装：**
\`\`\`
❌ Bohrium CLI未安装

请在PowerShell中运行以下命令安装：
curl -o install_bohr_windows_curl.bat https://dp-public.oss-cn-beijing.aliyuncs.com/bohrctl/1.0.0/install_bohr_windows_curl.bat && .\install_bohr_windows_curl.bat

安装完成后，重启VSCode，告诉我"我已重启VSCode，继续测试"
\`\`\`
- 创建状态文件`BOHRIUM_CONFIG_STATUS.md`（记录进度）
- **停止执行，等待用户重启**

---

**0.2.2 配置AccessKey**

\`\`\`bash
# 检查AccessKey
echo $env:ACCESS_KEY
\`\`\`

**如果为空：**
\`\`\`
❌ ACCESS_KEY未设置

请执行以下步骤：

1. 获取AccessKey
   - 访问 https://bohrium.dp.tech/
   - 登录后在"个人中心"→"AccessKey"获取

2. 在PowerShell中配置（请替换your_access_key_here）：
   setx ACCESS_KEY "your_access_key_here"

3. 重启VSCode使环境变量生效

完成后，告诉我"我已重启VSCode，继续测试"
\`\`\`

- 创建状态文件：
  \`\`\`markdown
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
  \`\`\`
- **停止执行，等待用户重启**

---

**0.2.3 验证配置（用户重启后）**

当用户说"我已重启VSCode"或类似表述时：

1. 读取状态文件（如果存在）
2. 执行验证命令：
\`\`\`bash
# 验证AccessKey
echo $env:ACCESS_KEY

# 测试API连接
bohr project list
\`\`\`

**预期输出：**
- AccessKey显示正确的值
- 项目列表显示用户的项目

**如果成功：**
\`\`\`
✅ Bohrium CLI配置成功！

可用项目：
- My Default Project (ID: 205855)

环境配置完成！
\`\`\`

**如果失败：**
- 检查AccessKey是否正确
- 确认用户已完全重启VSCode（不是重新加载窗口）
- 检查网络连接
- 提供故障排除建议（read `tools/testCLAUDE/troubleshooting.md`）

---

**0.2.4 配置项目ID（可选）**

如果用户有多个项目，设置默认项目：
\`\`\`bash
bohr config set project_id 205855
\`\`\`

---

完成后，删除状态文件（如果存在），输出"✅ Step 0完成，环境就绪"，返回主流程继续 Step 1。

---

## Bohrium CLI 命令速查

### 项目管理
\`\`\`bash
# 列出项目
bohr project list

# 设置默认项目
bohr config set project_id 205855
\`\`\`

### 任务管理
\`\`\`bash
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
\`\`\`

### 资源查询
\`\`\`bash
# 列出可用机型
bohr machine list

# 列出可用镜像
bohr image list

# 搜索ABACUS镜像
bohr image list | grep abacus
\`\`\`

### 配置管理
\`\`\`bash
# 查看当前配置
bohr config list

# 设置AccessKey
bohr config set access_key <your_key>

# 设置项目ID
bohr config set project_id <project_id>
\`\`\`
```

**Step 4: 验证文件创建**

```bash
(Get-Content tools/testCLAUDE/bohrium_setup.md).Count
```
预期：约 110 行

**Step 5: Commit**

```bash
git add tools/testCLAUDE/bohrium_setup.md
git commit -m "feat(testCLAUDE): 新增 bohrium_setup.md 模块（Step 0.2 配置流程 + CLI 命令速查）"
```

---

### Task 2: 编辑 testCLAUDE.md — 替换 Step 0.2 为桩代码

**Files:**
- Modify: `testCLAUDE.md`（约 lines 143–264）

**Step 1: Read 确认 old_string**

Read `testCLAUDE.md` lines 143–270，找到以下唯一锚点：

```
**情况C：Bohrium未配置 ⚠️**
- 进入"0.2 Bohrium完整配置流程"

---

**0.2 Bohrium完整配置流程（仅在未配置时执行）**
```
以及结尾：
```
bohr config set project_id 205855
```
（0.2.4 的最后一行），后面紧跟 `---`

**Step 2: Edit 替换**

用 Edit 工具，old_string = 从 `**情况C：Bohrium未配置 ⚠️**` 到 `bohr config set project_id 205855\n\`\`\`` 这整块（含中间所有行）。

new_string:
```markdown
**情况C：Bohrium未配置 ⚠️**
- 进入"0.2 Bohrium完整配置流程"

---

**0.2 Bohrium完整配置流程（仅在未配置时执行）**

→ **read `tools/testCLAUDE/bohrium_setup.md`** 获取完整配置流程和 Bohrium CLI 命令速查。
按其中步骤执行完成后，删除状态文件（如存在），输出"✅ Step 0完成，环境就绪"，继续 Step 1。
```

**Step 3: 验证行数减少**

```bash
(Get-Content testCLAUDE.md).Count
```
预期：约 1814 行（原 1934 - 约 120 行）

**Step 4: Commit**

```bash
git add testCLAUDE.md
git commit -m "refactor(testCLAUDE): Step 0.2 替换为按需读取桩（→ bohrium_setup.md）"
```

---

### Task 3: 创建 plugin_dev_guide.md

**Files:**
- Create: `tools/testCLAUDE/plugin_dev_guide.md`（内容来自 testCLAUDE.md Step 1.5）

**Step 1: Read 确认内容范围**

Read `testCLAUDE.md` 当前文件，搜索：
```
### Step 1.5：创建新测试插件（⭐ 自主扩展框架）
```
到：
```
（将 `<新插件名>`、`<calc_type关键词>`、`<计算类型>`、`<最低版本>` 替换为本次实际值）
```

记录实际行号（因 Task 2 已修改行数，行号已偏移）。

**Step 2: 写入 plugin_dev_guide.md**

Write `tools/testCLAUDE/plugin_dev_guide.md`，内容 = Step 1.5 全部内容（从 `### Step 1.5：` 标题到 issues_log 代码块最后的说明行），并在顶部加说明头：

```markdown
# 新测试插件自主创建指南

> 本文件由 testCLAUDE.md Step 1.4 按需 read，仅在 `detected_types = []` 时加载。
> 完成后返回主流程继续 Step 2。

---

```

然后追加从 testCLAUDE.md 读取的 Step 1.5 全部内容（`### Step 1.5：创建新测试插件...` 直到 issues_log 说明行结尾）。

**Step 3: 验证**

```bash
(Get-Content tools/testCLAUDE/plugin_dev_guide.md).Count
```
预期：约 158 行

**Step 4: Commit**

```bash
git add tools/testCLAUDE/plugin_dev_guide.md
git commit -m "feat(testCLAUDE): 新增 plugin_dev_guide.md 模块（Step 1.5 插件自主创建流程）"
```

---

### Task 4: 编辑 testCLAUDE.md — 替换 Step 1.5 为桩代码

**Files:**
- Modify: `testCLAUDE.md`

**Step 1: Read 确认 old_string**

Grep `testCLAUDE.md` 找 `### Step 1.5：创建新测试插件`，读取该行到 `（将 \`<新插件名>\`` 所在行（含）。

**Step 2: Edit 替换**

old_string = `### Step 1.5：创建新测试插件（⭐ 自主扩展框架）` 开头，到 `（将 \`<新插件名>\`、\`<calc_type关键词>\`、\`<计算类型>\`、\`<最低版本>\` 替换为本次实际值）` 结尾（含所有中间行）。

new_string:
```markdown
### Step 1.5：创建新测试插件（⭐ 自主扩展框架，仅 detected_types 为空时执行）

→ **read `tools/testCLAUDE/plugin_dev_guide.md`** 获取完整插件创建流程（含 1.5.1–1.5.4 步骤、代码模板、验证方法）。
→ **read `tools/testCLAUDE/plugins_history.md`** 参考已有插件的开发模式和历史记录。

完成插件创建、注册、验证后，在 `tools/testCLAUDE/plugins_history.md` 末尾追加一行记录，然后继续 Step 2。
```

**Step 3: 同步更新 Step 1.5 完成标志中的引用**

在 testCLAUDE.md 中，Step 1.5 原"完成标志"里有：
```
- 在 `testCLAUDE.md` **附录D** 的"插件开发历史"表中追加一行（含插件文件名、计算类型、添加日期、首次测试教程）
```

替换为：
```
- 在 `tools/testCLAUDE/plugins_history.md` 末尾追加一行（含插件文件名、计算类型、添加日期、首次测试教程）
```

注意：经过 Task 4 Step 2 替换后，此行实际已在 `plugin_dev_guide.md` 中（而非主文件），所以需要同时修改 `tools/testCLAUDE/plugin_dev_guide.md` 中的这行。

**Step 4: 验证行数**

```bash
(Get-Content testCLAUDE.md).Count
```
预期：约 1660 行（再减约 152 行）

**Step 5: Commit**

```bash
git add testCLAUDE.md tools/testCLAUDE/plugin_dev_guide.md
git commit -m "refactor(testCLAUDE): Step 1.5 替换为按需读取桩（→ plugin_dev_guide.md）"
```

---

### Task 5: 创建三个参考模块文件

**Files:**
- Create: `tools/testCLAUDE/plugins_history.md`
- Create: `tools/testCLAUDE/troubleshooting.md`
- Create: `tools/testCLAUDE/file_formats.md`

**Step 1: Read 内容范围（testCLAUDE.md 当前状态）**

Grep 以下唯一标题确认当前行号：
- `## 附录B: 故障排除清单`
- `## 附录C: 文件格式参考`
- `## 附录D: 未来计算类型扩展参考`

**Step 2: 写入 plugins_history.md**

内容 = 附录D全部（从 `## 附录D: 未来计算类型扩展参考` 到文档版本信息行），顶部加说明头：

```markdown
# 插件开发历史与扩展参考

> 本文件由 testCLAUDE.md 按需 read。
> **新增插件时，只需在本文件"插件开发历史"表末尾追加一行，无需修改主文件。**

---

```

然后追加附录D全部内容（保持原有结构不变）。

**Step 3: 写入 troubleshooting.md**

内容 = 附录B全部（从 `## 附录B: 故障排除清单` 到最后一个问题的解决方案），顶部加说明头：

```markdown
# 故障排除清单

> 本文件由 testCLAUDE.md 按需 read，在 Step 4.4（任务失败）或 Step 6.3（结果偏差）时加载。

---

```

然后追加附录B全部内容。

**Step 4: 写入 file_formats.md**

内容 = 附录C全部（从 `## 附录C: 文件格式参考` 到 STRU 格式说明末尾），顶部加说明头：

```markdown
# 文件格式参考

> 本文件由 testCLAUDE.md 按需 read，在 Step 2 生成输入文件或调试格式问题时加载。

---

```

然后追加附录C全部内容。

**Step 5: 验证三个文件**

```bash
(Get-Content tools/testCLAUDE/plugins_history.md).Count   # 预期 ~62 行
(Get-Content tools/testCLAUDE/troubleshooting.md).Count   # 预期 ~138 行
(Get-Content tools/testCLAUDE/file_formats.md).Count      # 预期 ~111 行
```

**Step 6: Commit**

```bash
git add tools/testCLAUDE/plugins_history.md tools/testCLAUDE/troubleshooting.md tools/testCLAUDE/file_formats.md
git commit -m "feat(testCLAUDE): 新增三个参考模块（plugins_history / troubleshooting / file_formats）"
```

---

### Task 6: 编辑 testCLAUDE.md — 替换四个附录为紧凑参考块

**Files:**
- Modify: `testCLAUDE.md`

**Step 1: Read 确认当前附录区域**

Read `testCLAUDE.md`，找到 `## 附录A: Bohrium CLI命令速查` 所在行到 `## 灵活性原则` 所在行之间的全部内容。

**Step 2: Edit 替换整个附录区**

用 Edit 工具，将从 `## 附录A: Bohrium CLI命令速查` 到附录D末尾（`**维护者：** AutoTutorial开发团队` 所在行）的全部内容，替换为以下紧凑参考块：

new_string:
```markdown
## 参考模块（按需读取）

以下内容已拆分为独立模块文件，在需要时通过 read 工具加载：

| 模块文件 | 内容 | 何时读取 |
|---------|------|---------|
| `tools/testCLAUDE/bohrium_setup.md` | Bohrium 配置流程 + CLI 命令速查 | Step 0 环境未配置时 |
| `tools/testCLAUDE/plugin_dev_guide.md` | 插件自主创建完整流程（Step 1.5）| detected_types 为空时 |
| `tools/testCLAUDE/plugins_history.md` | 插件开发历史 + 扩展场景参考 | Step 1.5 创建新插件时 |
| `tools/testCLAUDE/troubleshooting.md` | 8 类故障排除清单 | Step 4.4 任务失败 / Step 6.3 结果偏差时 |
| `tools/testCLAUDE/file_formats.md` | job.json / analysis.json / STRU 格式 | Step 2 生成输入文件时 |
```

**Step 3: 验证行数**

```bash
(Get-Content testCLAUDE.md).Count
```
预期：约 1290 行（再减约 370 行）

**Step 4: Commit**

```bash
git add testCLAUDE.md
git commit -m "refactor(testCLAUDE): 四个附录替换为紧凑参考表（→ tools/testCLAUDE/）"
```

---

### Task 7: 更新 testCLAUDE.md 内部的附录交叉引用

**Files:**
- Modify: `testCLAUDE.md`

**Step 1: 搜索所有附录引用**

```bash
grep -n "附录B\|附录C\|附录D\|appendix" testCLAUDE.md
```

预期找到以下位置（确认是否还有遗留引用）：
- Step 4.4 区域："参考故障排除清单（附录B）"
- Step 6.3 区域："参考故障排除清单（附录B）"
- 可能还有其他位置

**Step 2: Edit — Step 4.4 引用**

找到 Step 4.4 中包含附录B引用的那行，用 Edit 替换：

old_string（包含足够上下文确保唯一性）：
```
- 参考故障排除清单（附录B）
```
（如有多处，使用 Step 4.4 的上下文行确保唯一，如前后各带 1–2 行上下文）

new_string:
```
- 参考故障排除清单（read `tools/testCLAUDE/troubleshooting.md`）
```

**Step 3: Edit — Step 6.3 引用**

同上，替换 Step 6.3 中的同一表达。

**Step 4: 再次验证无遗留引用**

```bash
grep -n "附录[ABCD]" testCLAUDE.md
```
预期：0 匹配（所有附录引用已移除）

**Step 5: Commit**

```bash
git add testCLAUDE.md
git commit -m "refactor(testCLAUDE): 更新内部附录引用 → tools/testCLAUDE/ 路径"
```

---

### Task 8: 修改 CLAUDE.md Step 8 — 增加 compact 暂停点

**Files:**
- Modify: `CLAUDE.md`（Step 8 的"若用户选择'是'"部分）

**Step 1: Read 确认 old_string**

Read `CLAUDE.md`，找到 Step 8 中：
```
**若用户选择"是"：** 按以下方式启动测试流程：

1. 设置路径变量（根据当前工作目录自动填写）：
```tutorial_path = "_workspace/<当前任务目录>/07_Final_Tutorial_<标题>.md"
test_dir      = "_workspace/<当前任务目录>/test_<YYYYMMDD_HHMMSS>/"
```

2. 切换到 `testCLAUDE.md` 定义的完整测试流程，从 **Step 1** 开始执行，`tutorial_path` 和 `test_dir` 已知，无需再次询问。
```

**Step 2: Edit 替换**

old_string = 上述"若用户选择'是'"开始到"无需再次询问。"结尾的完整段落。

new_string:
```markdown
**若用户选择"是"：** 按以下方式启动测试流程：

1. 设置路径变量（根据当前工作目录自动填写）：
```
tutorial_path = "_workspace/<当前任务目录>/07_Final_Tutorial_<标题>.md"
test_dir      = "_workspace/<当前任务目录>/test_<YYYYMMDD_HHMMSS>/"
```

2. 输出以下信息，暂停等待用户响应：
```
即将开始计算测试。关键信息已确认：
  tutorial_path = _workspace/<当前目录>/07_Final_Tutorial_<标题>.md
  test_dir      = _workspace/<当前目录>/test_<时间戳>/

建议：运行 /compact 压缩当前上下文（可释放 40–60K tokens），有助于测试流程顺利运行。
压缩后告诉我"继续测试"即可。若不压缩，直接说"继续测试"也可。
```

3. 收到"继续测试"后，切换到 `testCLAUDE.md` 定义的完整测试流程，从 **Step 1** 开始执行，`tutorial_path` 和 `test_dir` 已知，无需再次询问。
```

**Step 3: 验证修改**

Read `CLAUDE.md` 确认 Step 8 新内容已正确写入，且"compact"字样出现在合适位置。

**Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "feat(CLAUDE): Step 8 增加 compact 暂停点，输出路径变量后等待用户确认"
```

---

### Task 9: 最终验证

**Step 1: 行数检查**

```bash
(Get-Content testCLAUDE.md).Count          # 预期 1250–1300 行
(Get-Content tools/testCLAUDE/bohrium_setup.md).Count
(Get-Content tools/testCLAUDE/plugin_dev_guide.md).Count
(Get-Content tools/testCLAUDE/plugins_history.md).Count
(Get-Content tools/testCLAUDE/troubleshooting.md).Count
(Get-Content tools/testCLAUDE/file_formats.md).Count
```

**Step 2: 关键内容检查**

```bash
# 确认主文件中桩代码存在
grep -n "bohrium_setup.md\|plugin_dev_guide.md\|plugins_history.md\|troubleshooting.md\|file_formats.md" testCLAUDE.md

# 确认无遗留附录标题
grep -n "## 附录[ABCD]" testCLAUDE.md

# 确认模块文件可读
head -5 tools/testCLAUDE/bohrium_setup.md
head -5 tools/testCLAUDE/plugins_history.md
```

**Step 3: 检查 CLAUDE.md**

```bash
grep -n "compact\|继续测试" CLAUDE.md
```
预期：在 Step 8 区域找到这两个关键词。

**Step 4: 最终 Commit（如有未提交内容）**

```bash
git status
git add -A
git commit -m "chore: testCLAUDE 模块化拆分完成验证"
```
