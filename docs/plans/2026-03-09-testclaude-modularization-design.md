# testCLAUDE.md 模块化拆分设计文档

**日期：** 2026-03-09
**目标：** 将 testCLAUDE.md（1934 行）拆分为主文件 + 5 个按需读取模块，减轻测试执行时的上下文压力，同时提升文件可维护性。

---

## 背景与问题

- `testCLAUDE.md` 当前共 1934 行，单次 Read 调用全部进入上下文，约占 **40K tokens**
- 测试流程执行时，上下文还包含：CLAUDE.md（~15K tokens）+ 教程生成的对话历史（~30–60K tokens）
- 实际运行中已出现自动 compact，说明 200K 窗口已接近上限
- compact 发生在关键步骤（Step 1.5 写插件、Step 3–4 监控任务），会导致状态丢失

## 设计决策

### 模块拆分原则

1. **按需读取**：模块只在需要时由 Claude 显式 read，不在会话开始时加载
2. **生命周期一致**：合并生命周期相同的内容（如 Bohrium CLI 命令随 Step 0 读入）
3. **功能命名**：文件名直接反映功能，不使用 appendix_A/B/C 等编号
4. **集中存放**：所有模块放入 `tools/testCLAUDE/` 目录

### 文件结构

```
testCLAUDE.md                         # 主文件，~1280 行（原 1934 行）
tools/testCLAUDE/
  bohrium_setup.md                    # Step 0 完整配置流程 + Bohrium CLI 命令速查
  plugin_dev_guide.md                 # Step 1.5 插件自主创建流程
  plugins_history.md                  # 插件开发历史表（每新增插件追加一行）
  troubleshooting.md                  # 故障排除清单（8 类场景）
  file_formats.md                     # 文件格式参考（job.json / analysis.json / STRU）
```

### 各模块内容来源

| 模块文件 | 来源（原 testCLAUDE.md）| 行数 |
|---------|----------------------|------|
| `bohrium_setup.md` | Step 0 完整配置流程（Lines 146–261）+ Appendix A（Lines 1561–1613）| ~238 行 |
| `plugin_dev_guide.md` | Step 1.5 插件自主创建流程（Lines 398–549）| ~152 行 |
| `plugins_history.md` | Appendix D 插件开发历史（Lines 1855–1909）| ~55 行 |
| `troubleshooting.md` | Appendix B 故障排除清单（Lines 1616–1746）| ~131 行 |
| `file_formats.md` | Appendix C 文件格式参考（Lines 1749–1852）| ~104 行 |

### 触发时机

| 模块 | 触发条件 | 典型频率 |
|-----|---------|--------|
| `bohrium_setup.md` | Step 0 检测到环境未配置 | 同一台机器仅首次 |
| `plugin_dev_guide.md` | Step 1.4 发现 `detected_types = []` | 新计算类型时 |
| `plugins_history.md` | Step 1.4 创建新插件，需参考历史模式 | 同上 |
| `troubleshooting.md` | Step 4.4 任务失败 / Step 6.3 结果偏差 | 出错时 |
| `file_formats.md` | Step 2 生成输入文件，需格式参考 | 偶尔 |

正常流程（环境已配置、插件已存在、测试通过）只读主文件 ~1280 行。

### 主文件中的桩代码模式

原有大块内容替换为统一格式的"桩"：

```markdown
#### Step 0.2 Bohrium 配置（仅环境未配置时执行）

检测到环境未配置，read `tools/testCLAUDE/bohrium_setup.md` 获取完整配置流程，
按其中的步骤执行后返回继续 Step 1。
```

```markdown
#### Step 1.5 自主创建新插件（仅 detected_types 为空时执行）

read `tools/testCLAUDE/plugin_dev_guide.md` 获取完整插件创建流程。
同时 read `tools/testCLAUDE/plugins_history.md` 参考已有插件的开发模式。
```

```markdown
#### 附录参考

- 故障排除：需要时 read `tools/testCLAUDE/troubleshooting.md`
- 文件格式：需要时 read `tools/testCLAUDE/file_formats.md`
- 插件历史：需要时 read `tools/testCLAUDE/plugins_history.md`
```

---

## CLAUDE.md 修改（Step 8 暂停点）

在用户同意执行测试后、读取 testCLAUDE.md 之前，增加显式暂停点：

**Claude 输出内容：**
```
即将开始计算测试。关键信息已确认：
  tutorial_path = _workspace/<当前目录>/07_Final_Tutorial_<标题>.md
  test_dir      = _workspace/<当前目录>/test_<时间戳>/

建议：运行 /compact 压缩当前上下文（可释放 40–60K tokens），
压缩后告诉我"继续测试"即可。若不压缩，直接说"继续测试"也可。
```

**原理：**
- Claude 无法主动触发 `/compact`（无对应工具接口）
- 通过明确输出路径变量，确保 compact 后摘要中保留关键信息
- compact 是可选的，不阻塞正常流程

---

## 预期效果

| 指标 | 改造前 | 改造后（正常流程）|
|------|--------|-----------------|
| testCLAUDE.md 行数 | 1934 行 | ~1280 行 |
| 读入上下文的 token 数 | ~40K | ~26K |
| compact 后关键状态保留 | 依赖摘要质量 | 路径变量在 compact 前明确输出 |
| 新增插件时需修改的文件 | testCLAUDE.md（主文件）| 仅 plugins_history.md |

---

## 不在此次改造范围内

- 不拆分 Step 2–8 的核心流程（这些是顺序执行的，拆分收益低于管理成本）
- 不修改任何 Python 工具文件
- 不改变测试逻辑、验证容差、Bohrium 提交参数
