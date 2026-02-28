# Bohrium CLI 配置指南（Windows）

## 📋 概述

本指南帮助Claude Code在Windows系统上快速配置Bohrium CLI，用于提交ABACUS计算任务。

**关键点**：
- ✅ 环境变量设置后**必须重启VSCode**才能生效
- ✅ 重启后从"验证配置"步骤继续
- ✅ 整个流程约10分钟

---

## 🚀 快速配置流程

### Step 1: 检查安装

```bash
bohr version
```

**预期输出**：`1.1.0` 或其他版本号

**如果未安装**，请用户在PowerShell中运行：
```powershell
curl -o install_bohr_windows_curl.bat https://dp-public.oss-cn-beijing.aliyuncs.com/bohrctl/1.0.0/install_bohr_windows_curl.bat && .\install_bohr_windows_curl.bat
```

---

### Step 2: 配置AccessKey

**检查环境变量**：
```bash
echo $ACCESS_KEY
```

**如果为空**，请用户在PowerShell中运行：
```powershell
setx ACCESS_KEY "your_access_key_here"
```

**⚠️ 关键步骤**：
1. 告诉用户：**必须重启VSCode**才能使环境变量生效
2. 在用户重启前，记录当前进度到工作目录
3. 创建状态文件，标记"等待重启"

---

### Step 3: 请求用户重启

**创建状态文件**（在工作目录）：
```markdown
# Bohrium配置状态

## 当前进度
- [x] 检查安装
- [x] 配置AccessKey
- [ ] 验证配置（等待重启）

## 下一步
重启VSCode后，告诉Claude：
"我已重启VSCode，继续Bohrium配置"

Claude将从"验证配置"步骤继续。
```

**告诉用户**：
```
✅ AccessKey已配置！

⚠️ 重要：请重启VSCode使环境变量生效

重启后，告诉我"我已重启VSCode，继续Bohrium配置"，
我将继续验证配置并开始测试。
```

**停止执行**，等待用户重启并返回。

---

### Step 4: 验证配置（重启后）

当用户返回后，运行：

```bash
# 检查环境变量
echo $ACCESS_KEY

# 测试API连接
bohr project list
```

**预期输出**：
- 环境变量显示AccessKey值
- 项目列表显示用户的项目

**如果成功**：
```
✅ Bohrium CLI配置成功！

可用功能：
- bohr project list    # 列出项目
- bohr machine list    # 列出机型
- bohr image list      # 列出镜像
- bohr job submit      # 提交任务
```

**如果失败**：
- 检查AccessKey是否正确
- 确认用户已重启VSCode
- 检查网络连接

---

## 🔄 重启后如何继续

### 识别重启场景

用户可能说：
- "我已重启VSCode"
- "重启完成"
- "继续Bohrium配置"
- "继续之前的任务"

### 恢复进度

1. **读取状态文件**（如果存在）
2. **检查环境变量**：`echo $ACCESS_KEY`
3. **从验证步骤继续**：`bohr project list`

### 状态文件位置

- 工作目录：`_workspace/YYYYMMDD_HHMMSS_主题/bohrium_status.md`
- 或根目录：`BOHRIUM_CONFIG_STATUS.md`

---

## 📝 完整示例对话

### 首次配置

**Claude**:
```
让我检查Bohrium CLI安装状态...
[运行 bohr version]
✅ 已安装版本 1.1.0

现在检查AccessKey配置...
[运行 echo $ACCESS_KEY]
❌ 未配置

请在PowerShell中运行以下命令：
setx ACCESS_KEY "your_access_key_here"

然后重启VSCode，告诉我"我已重启VSCode，继续Bohrium配置"
```

**用户**: [重启VSCode]

**用户**: "我已重启VSCode，继续Bohrium配置"

**Claude**:
```
欢迎回来！让我验证配置...
[运行 echo $ACCESS_KEY]
✅ AccessKey已设置

[运行 bohr project list]
✅ API连接成功！

找到项目：
- My Default Project (2 items)

配置完成！现在可以提交ABACUS任务了。
```

---

## 🎯 提交ABACUS任务

配置完成后，可以提交任务：

```bash
# 列出可用机型
bohr machine list

# 列出ABACUS镜像
bohr image list | grep abacus

# 提交任务（示例）
bohr job submit \
  --job-name "Si-elastic" \
  --command "OMP_NUM_THREADS=1 mpirun -np 8 abacus" \
  --machine c16_m32_cpu \
  --image registry.dp.tech/dptech/abacus:3.0.0 \
  --input-data ./Si_elastic
```

---

## ⚠️ 常见问题

### Q1: 重启后仍然找不到AccessKey
**A**:
- 确认在PowerShell中运行了`setx`命令
- 确认完全关闭并重新打开VSCode（不是重新加载窗口）
- 检查Windows环境变量设置

### Q2: bohr命令找不到
**A**:
- 确认已安装bohrctl
- 重启VSCode
- 或手动刷新PATH

### Q3: API连接失败
**A**:
- 检查AccessKey是否正确
- 检查网络连接
- 尝试访问 https://bohrium.dp.tech/

---

## 📊 配置检查清单

- [ ] `bohr version` 显示版本号
- [ ] `echo $ACCESS_KEY` 显示AccessKey
- [ ] `bohr project list` 显示项目列表
- [ ] `bohr machine list` 显示机型列表
- [ ] `bohr image list` 显示镜像列表

全部通过后，配置完成！

---

## 🔧 Claude Code使用指南

### 何时请求重启

当执行以下操作后：
1. 使用`setx`设置环境变量
2. 安装新的CLI工具
3. 修改系统PATH

### 如何请求重启

```markdown
⚠️ 重要：请重启VSCode使配置生效

重启后，告诉我"我已重启VSCode，继续XXX任务"
我将从[具体步骤]继续。

[可选：创建状态文件，记录当前进度]
```

### 重启后如何恢复

1. 用户说"已重启"或类似表述
2. 读取状态文件（如果有）
3. 验证配置是否生效
4. 从中断点继续执行

---

## 📚 参考链接

### 快速入门
- Status：https://bohrium-doc.dp.tech/docs/quickstart/Status
- KillJob：https://bohrium-doc.dp.tech/docs/quickstart/KillJob
- Result：https://bohrium-doc.dp.tech/docs/quickstart/Result

### 命令行工具
- Bohrium CLI 简介：https://docs.bohrium.com/docs/bohrctl/about
- 安装Bohrium CLI：https://docs.bohrium.com/docs/bohrctl/install
- 任务：https://docs.bohrium.com/docs/bohrctl/job
- 任务组：https://docs.bohrium.com/docs/bohrctl/job_group
- 节点：https://docs.bohrium.com/docs/bohrctl/node
- 机型：https://docs.bohrium.com/docs/bohrctl/machine
- 数据集：https://docs.bohrium.com/docs/bohrctl/dataset
- 镜像：https://docs.bohrium.com/docs/bohrctl/image
- 项目：https://docs.bohrium.com/docs/bohrctl/project

### 软件案例
- ABACUS：https://bohrium-doc.dp.tech/docs/software/ABACUS/

### 常见问题
- CLI：https://bohrium-doc.dp.tech/docs/faq/UtilityFaq
- 节点：https://bohrium-doc.dp.tech/docs/faq/Machine
- 任务：https://bohrium-doc.dp.tech/docs/faq/Job

---

**创建日期**: 2026-02-06
**适用系统**: Windows 10/11
**Bohrium CLI版本**: 1.1.0+
**预计配置时间**: 10分钟（含重启）
