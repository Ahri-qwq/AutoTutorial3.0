# AutoTutorial 3.0

**AutoTutorial 3.0** 是一个专为 **ABACUS** 第一性原理计算软件设计的智能教程生成系统。通过 **Claude Code** 协调 RAG 检索和案例解析工具，自动生成高质量、包含实战参数配置的物理教程。

> **2026年2月 重大更新**：
> 系统从自动化pipeline模式升级为 **Claude Code 协调模式**，通过智能对话和工具调用，实现更灵活、更可控的教程生成流程。

---

## 🌟 核心特性

- **🤖 Claude Code 协调**：使用 Claude Code 作为主控，通过自然语言对话完成教程生成
- **🔍 RAG 知识检索**：从本地知识库检索相关文档，确保内容准确性
- **📂 案例驱动生成**：支持上传案例文件（.docx/.md/.txt），自动提取参数和结构
- **📝 多轮审查机制**：内容审查、案例审查、风格审查，确保教程质量
- **🎯 Think Aloud**：全程透明化思考过程，用户可见每一步决策

---

## 🛠️ 快速开始

### 1. 环境配置

确保已安装 Python 3.8+，并安装依赖：

```bash
pip install -r requirements.txt
```

### 2. 配置文件

在项目根目录创建 `config.yaml`：

```yaml
llm:
  api_key: "sk-xxxxxxxxxxxx"  # 阿里云 DashScope API Key（用于Embedding）
```

### 3. 初始化知识库

将 ABACUS 相关资料（PDF/Markdown/Txt）放入 `data/knowledge_source/` 目录，然后运行：

```bash
python src/knowledge_manager.py
```

### 4. 添加参考文章

将优质的 ABACUS 教程文章放入 `data/reference_materials/style_references/` 目录，用于学习写作风格。

### 5. 开始使用

在 Claude Code 中输入：

```
请按照CLAUDE.md生成一篇关于"ABACUS SCF自洽计算"的教程
```

Claude Code 会：
1. 询问任务信息（是否有案例等）
2. 创建工作目录
3. 执行知识调研（调用 retriever.py）
4. 提供 3 个大纲方案供选择
5. 撰写初稿
6. 三轮审查（内容、案例、风格）
7. 生成最终版本

---

## 📚 使用方式

### 方式1：标准教程生成

适合没有现成案例，需要从零开始生成教程的场景。

```
你：请按照CLAUDE.md生成一篇关于"ABACUS能带计算"的教程

Claude Code会询问：
- 是否有案例文件？
- 是否有特殊要求？

然后执行完整的7步流程
```

### 方式2：基于案例的教程

适合有现成案例，希望在教程中引用案例的场景。

```
你：
请按照CLAUDE.md生成一篇关于"ABACUS弹性常数计算"的教程，
案例文件是 data/input/elastic_case.docx

Claude Code会：
1. 解析案例（调用 case_parser.py）
2. 在教程中引用案例数据
3. 确保案例数据完整准确
```

### 方式3：案例驱动教程

适合整篇教程围绕案例展开的场景。

```
你：
请按照CLAUDE.md生成一篇关于"ABACUS弹性常数计算"的教程，
案例文件是 data/input/elastic_case.docx，
请围绕案例展开

Claude Code会：
1. 以案例为核心设计大纲
2. 所有参数说明都基于案例
3. 确保案例数据完整性
```

---

## 📁 目录结构

```
AutoTutorial3.0/
├── CLAUDE.md                    # Claude Code 指南（核心）
├── 快速开始.md                  # 快速开始指南
├── README.md                    # 本文件
│
├── tools/                       # 工具
│   ├── retriever.py            # RAG 检索工具
│   ├── case_parser.py          # 案例解析工具
│   └── README.md               # 工具使用说明
│
├── data/
│   ├── knowledge_source/       # 知识源文件（PDF/MD/TXT）
│   ├── chroma_db/              # 向量数据库
│   ├── reference_materials/    # 参考文章
│   │   └── style_references/   # 风格参考文章
│   └── input/                  # 输入案例文件
│
├── _workspace/                 # 工作区（自动创建）
│   └── YYYYMMDD_HHMMSS_主题/  # 每次生成的工作目录
│       ├── 00_brief.md        # 任务简报
│       ├── 01_research.md     # 知识调研结果
│       ├── 02_outline.md      # 大纲方案
│       ├── 03_chapter_*.md    # 各章节初稿
│       ├── 04_review_content.md  # 内容审查
│       ├── 05_review_case.md     # 案例审查
│       ├── 06_review_style.md    # 风格审查
│       └── 07_final.md           # 最终版本
│
└── config.yaml                 # 配置文件
```

---

## 🔧 工具说明

### retriever.py - RAG 检索工具

从知识库检索相关文档。

```bash
python tools/retriever.py --query "ABACUS SCF计算 参数设置" --top_k 10
```

### case_parser.py - 案例解析工具

解析案例文件，提取文件结构、参数、计算流程。

```bash
python tools/case_parser.py --input "data/input/case.docx"
```

详细使用说明见 `tools/README.md`。

---

## 📝 核心原则

1. **真实性优先**：所有参数、数值必须来自检索结果或案例，严禁编造
2. **案例完整性**：案例中的所有数据必须完整、准确地出现在文中
3. **Think Aloud**：每一步都说明思考过程，让用户看到决策过程
4. **先讨论大纲**：提供 3 个大纲方案供用户选择，避免方向错误
5. **多轮审查**：内容审查、案例审查、风格审查，确保质量

---

## 🎯 工作流程

### Step 0: 任务初始化
- 询问任务信息
- 判断任务类型（新建/基于案例/案例驱动/修改）
- 创建工作目录

### Step 1: 知识调研与案例解析
- RAG 知识检索
- 案例解析（如有案例）
- 风格参考学习

### Step 2: 大纲讨论
- 提供 3 个大纲方案
- 等待用户选择
- 根据反馈调整

### Step 3: 撰写初稿
- 逐章撰写
- 使用双查询 RAG 检索
- 保存独立章节文件

### Step 4: 内容审查
- 逻辑审查
- 事实审查
- 结构审查

### Step 5: 案例审查（如有案例）
- 案例完整性检查
- 案例数据准确性检查

### Step 6: 风格审查
- AI 腔检测与修改
- 语言风格对齐
- 细节打磨

### Step 7: 最终输出
- 整合所有修改
- 添加前言和附录
- 生成元数据

---

## ❓ 常见问题

### Q1: 我没有参考文章怎么办？

可以：
1. 从 `data/runs/` 中找之前生成的好的教程
2. 从网上找 ABACUS 相关教程
3. 或者先跳过，直接测试（但风格可能不够好）

### Q2: 知识库如何更新？

将新的资料放入 `data/knowledge_source/`，然后重新运行：

```bash
python src/knowledge_manager.py
```

### Q3: 如何查看生成的教程？

生成的教程在：

```
_workspace/YYYYMMDD_HHMMSS_主题/07_final.md
```

### Q4: 可以修改 CLAUDE.md 吗？

可以！CLAUDE.md 是灵活的指南，你可以根据需要调整流程、写作规范等。但建议保留核心原则（真实性、案例完整性、Think Aloud 等）。

---

## 📊 依赖说明

- **chromadb**: 向量数据库
- **dashscope**: 阿里云 DashScope SDK（用于 Qwen Embedding）
- **python-docx**: 解析 Word 文档
- **PyYAML**: 配置文件解析
- **pypdf**: PDF 文件解析（可选）

---

## 🔄 版本历史

### v3.0 (2026.02)
- 从自动化 pipeline 模式升级为 Claude Code 协调模式
- 简化依赖，移除 OpenAI/Gemini 相关依赖
- 优化工具接口，提升易用性

### v2.0 (2026.01)
- 生成引擎从 Qwen 迁移到 Gemini
- 新增案例驱动生成功能
- 实施分步温度控制策略

---

## 📄 许可证

MIT License

---

## 🤝 贡献

欢迎提交 Issue 和 Pull Request！

---

**现在，开始使用 Claude Code 生成你的第一篇 ABACUS 教程吧！** 🚀
