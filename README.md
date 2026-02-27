# AutoTutorial 3.0

**AutoTutorial 3.0** 是一个专为 **ABACUS** 第一性原理计算软件设计的智能教程生成与验证系统。通过 **Claude Code** 协调 RAG 检索、案例解析和计算测试工具，自动生成高质量教程，并可将教程提交到 Bohrium 云平台进行实际计算验证。

---

## 核心特性

- **Claude Code 协调**：使用 Claude Code 作为主控，通过自然语言对话完成教程生成全流程
- **RAG 知识检索**：从本地知识库检索相关文档，确保内容准确性
- **案例驱动生成**：支持上传案例文件（.docx/.md），自动提取参数和文件结构
- **多轮审查机制**：内容审查、案例审查、风格审查，确保教程质量
- **自动文件验证**：`orbital_validator.py` 在生成阶段自动修正错误轨道文件名、删除不兼容参数（如 `nbands auto`）
- **计算测试验证**：教程生成完成后，可选择提交到 Bohrium 云平台实际运行，验证参数和流程的可复现性
- **反馈写回机制**：测试发现的问题（文件名错误、参数不兼容）自动写回教程原文

---

## 快速开始

### 1. 环境配置

```bash
pip install -r requirements.txt
```

### 2. 配置文件

在项目根目录创建 `config.yaml`：

```yaml
llm:
  api_key: "sk-xxxxxxxxxxxx"  # 阿里云 DashScope API Key（用于 Embedding）
```

### 3. 初始化知识库

将 ABACUS 相关资料（PDF/Markdown/Txt）放入 `data/knowledge_source/`，然后运行：

```bash
python src/knowledge_manager.py
```

### 4. 添加参考文章

将优质的 ABACUS 教程文章放入 `data/reference_materials/style_references/`，用于学习写作风格。

### 5. （可选）配置 Bohrium 测试

如需使用计算测试功能，配置 Bohrium CLI：

```bash
pip install bohrium-cli
bohr config set access_key <your-access-key>
bohr config set project_id <your-project-id>
```

### 6. 开始使用

在 Claude Code 中输入：

```
请按照CLAUDE.md生成一篇关于"ABACUS SCF自洽计算"的教程
```

---

## 工作流程

系统分为两个阶段：**教程生成**（CLAUDE.md）和**计算测试**（testCLAUDE.md）。

### 阶段一：教程生成（CLAUDE.md）

| 步骤 | 内容 |
|------|------|
| Step 0 | 任务初始化，创建工作目录 |
| Step 1 | RAG 知识检索、案例解析、风格参考学习 |
| Step 2 | 提供 3 个大纲方案，等待用户选择 |
| Step 3 | 逐章撰写初稿（每章双查询 RAG） |
| Step 4 | 内容审查（逻辑、事实、结构） |
| Step 5 | 案例审查（完整性、准确性，如有案例） |
| Step 6 | 风格审查（去 AI 腔、语言打磨） |
| Step 7 | 最终输出（含轨道文件名自动验证） |
| **Step 8** | **询问是否执行计算测试，确认后切换到 testCLAUDE.md** |

### 阶段二：计算测试（testCLAUDE.md）

| 步骤 | 内容 |
|------|------|
| Step 0 | 初始化 Bohrium 配置 |
| Step 1 | 调用 `--phase prepare` 解析教程、准备输入文件 |
| Step 2 | 验证输入文件（轨道文件、INPUT 参数、STRU 格式） |
| Step 3 | Claude 读教程，Think Aloud 分析任务依赖，手工提交 Bohrium 任务 |
| Step 4 | 监控任务状态 |
| Step 5 | 下载计算结果 |
| Step 6 | 对比结果，生成测试报告 |
| Step 7 | 将测试发现的问题反向写回教程原文 |

---

## 目录结构

```
AutoTutorial3.0/
├── CLAUDE.md                        # 教程生成指南（核心）
├── testCLAUDE.md                    # 计算测试指南
├── README.md                        # 本文件
│
├── tools/
│   ├── retriever.py                 # RAG 检索工具
│   ├── case_parser.py               # 案例解析工具
│   ├── orbital_validator.py         # 轨道文件名 + INPUT 参数验证工具
│   ├── orbital_db.py                # 轨道文件映射数据库
│   ├── fix_stru.py                  # STRU 文件格式修复工具
│   ├── test_framework_integrated.py # 测试框架主入口
│   ├── test_framework_phase1_analyzer.py
│   ├── test_framework_phase3_7_impl.py
│   ├── test_plugins/                # 计算类型插件
│   │   ├── relax_plugin.py          # 结构优化
│   │   ├── elastic_plugin.py        # 弹性常数
│   │   ├── band_plugin.py           # 能带计算
│   │   └── dos_plugin.py            # 态密度
│   ├── orbitals/                    # 轨道文件本地缓存
│   └── pseudopotentials/            # 赝势文件本地缓存
│
├── data/
│   ├── knowledge_source/            # 知识库源文件（PDF/MD/TXT）
│   ├── chroma_db/                   # 向量数据库
│   ├── reference_materials/
│   │   └── style_references/        # 风格参考文章
│   └── input/                       # 输入案例文件
│
├── _workspace/                      # 工作区（自动创建）
│   └── YYYYMMDD_HHMMSS_主题/
│       ├── process/
│       │   ├── 00_brief.md
│       │   ├── 01_research.md
│       │   ├── 02_outline.md
│       │   ├── 03_chapter_*.md
│       │   ├── 04_review_content.md
│       │   ├── 05_review_case.md
│       │   ├── 06_review_style.md
│       │   └── 07_fix.md
│       ├── 07_Final_Tutorial_<标题>.md   # 最终教程
│       └── test_<时间戳>/               # 测试工作目录（如执行测试）
│           ├── 01_analysis.json
│           ├── <案例名>/               # 计算输入文件
│           ├── orbital_fix_report.json
│           ├── param_fix_report.json
│           └── test_report.md
│
└── config.yaml                      # 配置文件
```

---

## 工具说明

### orbital_validator.py — 轨道文件名 + INPUT 参数验证

在教程生成的 Step 7 自动调用，检测并修正两类问题：

```bash
python tools/orbital_validator.py process/07_fix.md --fix
```

- **轨道文件名错误**：如 `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`
- **INPUT 参数不兼容**：如 `nbands auto` → 删除（ABACUS v3.10.x 不支持）

### test_framework_integrated.py — 计算测试框架

```bash
# 仅准备输入文件（Claude 控制任务提交）
python tools/test_framework_integrated.py tutorial.md --test-dir ./test_dir --phase prepare

# 任务完成后继续后处理
python tools/test_framework_integrated.py continue ./test_dir
```

### retriever.py — RAG 检索

```bash
python tools/retriever.py --query "ABACUS SCF计算 参数设置" --top_k 10
```

### case_parser.py — 案例解析

```bash
python tools/case_parser.py --input "data/input/case.docx"
```

### fix_stru.py — STRU 格式修复

自动将旧版 STRU 格式（4字段）转为 ABACUS v3.x 新格式（3字段 + NUMERICAL_ORBITAL 块）：

```bash
python tools/fix_stru.py --stru path/to/STRU
```

---

## 使用示例

### 生成教程（无案例）

```
请按照CLAUDE.md生成一篇关于"ABACUS能带计算"的教程
```

### 生成教程（基于案例）

```
请按照CLAUDE.md生成一篇关于"ABACUS弹性常数计算"的教程，
案例文件是 data/input/elastic_case.docx，请围绕案例展开
```

生成完成后，Claude 会询问是否执行计算测试。选择"是"后自动切换到 testCLAUDE.md 流程。

---

## 核心原则

1. **真实性优先**：所有参数、数值必须来自检索结果或案例，严禁编造
2. **案例完整性**：案例中的所有数据必须完整、准确地出现在文中
3. **Think Aloud**：每一步都说明思考过程，用户可见每一步决策
4. **先讨论大纲**：提供 3 个大纲方案供用户选择，避免方向错误
5. **多轮审查**：内容审查、案例审查、风格审查，确保质量

---

## 依赖说明

- **chromadb**：向量数据库
- **dashscope**：阿里云 DashScope SDK（用于 Qwen Embedding）
- **python-docx**：解析 Word 文档
- **PyYAML**：配置文件解析
- **pypdf**：PDF 文件解析（可选）
- **bohrium-cli**：Bohrium 云平台 CLI（计算测试功能可选）

---

## 版本历史

### v3.0 (2026.02)
- 从自动化 pipeline 模式升级为 Claude Code 协调模式
- 新增完整的计算测试验证系统（testCLAUDE.md + 插件化测试框架）
- 新增 `orbital_validator.py`：自动修正轨道文件名错误和 INPUT 参数兼容性问题
- 新增 `fix_stru.py`：自动修复 STRU 文件格式
- 新增 `--phase prepare` 模式：Python 仅负责文件准备，Claude 控制所有任务提交决策
- 新增测试反馈写回机制（Step 7）：测试通过后自动将修正写回教程原文

### v2.0 (2026.01)
- 生成引擎从 Qwen 迁移到 Gemini
- 新增案例驱动生成功能

---

## 许可证

MIT License
