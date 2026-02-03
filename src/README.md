# knowledge_manager.py - 知识库管理工具

## 功能说明

`knowledge_manager.py` 是 AutoTutorial 3.0 的知识库管理工具，用于构建和维护 ABACUS 知识向量数据库。

**核心功能：**
- 从文档文件构建向量数据库（支持 .md/.txt/.docx/.pdf）
- 使用阿里云 DashScope Embedding API 生成向量
- 支持两种模式：重建（rebuild）和增量追加（append）

## 使用方法

### 1. 重建知识库（首次使用或完全重建）

将所有知识源文件放入 `data/knowledge_source/` 目录，然后运行：

```bash
python src/knowledge_manager.py
```

选择 `1` 进入重建模式，确认后将：
- 清空现有数据库
- 读取 `data/knowledge_source/` 中所有文件
- 生成向量并存入 `data/chroma_db/`

### 2. 增量添加知识（添加新文档）

将新文档放入 `data/knowledge_add/` 目录，然后运行：

```bash
python src/knowledge_manager.py
```

选择 `2` 进入增量模式，将：
- 读取 `data/knowledge_add/` 中的新文件
- 生成向量并追加到现有数据库
- 自动将成功入库的文件归档到 `data/knowledge_source/`

## 配置要求

确保 `config.yaml` 中配置了阿里云 API Key：

```yaml
llm:
  api_key: "sk-xxxxxxxxxxxx"  # 阿里云 DashScope API Key
```

## 支持的文件格式

- `.md` - Markdown 文件（按标题切分）
- `.txt` - 纯文本文件（按段落切分）
- `.docx` - Word 文档
- `.pdf` - PDF 文档

## 注意事项

- 重建模式会清空现有数据库，请谨慎使用
- 增量模式会自动移动文件，确保 `knowledge_add/` 中的文件可以被移动
- Embedding 失败的文件会被跳过，不会污染索引
