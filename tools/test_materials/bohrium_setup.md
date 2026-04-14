# Bohrium 配置流程与 CLI 命令速查

> 本文件由 testCLAUDE.md 按需 read，仅在 Step 0 检测到 Bohrium 未配置时加载。

---

## Step 0.2 Bohrium完整配置流程

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
- 提供故障排除建议（read `tools/testCLAUDE/troubleshooting.md`）

---

**0.2.4 配置项目ID（可选）**

如果用户有多个项目，设置默认项目：
```bash
bohr config set project_id 205855
```

---

完成后，删除状态文件（如果存在），输出"✅ Step 0完成，环境就绪"，返回主流程继续 Step 1。

---

## Bohrium CLI 命令速查

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
