# Tools 使用说明

## retriever.py - RAG检索工具

### 功能
从知识库检索相关文档，用于教程生成的知识调研阶段。

### 使用方法

**基本用法：**
```bash
python tools/retriever.py --query "查询词" --top_k 10
```

**参数说明：**
- `--query`：检索查询词（必需）
- `--top_k`：返回文档数量（默认：5）
- `--quiet`：安静模式，只输出文档内容，不输出格式化信息

### 使用示例

**示例1：检索物理原理**
```bash
python tools/retriever.py --query "弹性常数计算 物理原理 参数设置" --top_k 10
```

**示例2：检索常见问题**
```bash
python tools/retriever.py --query "ABACUS 收敛问题 注意事项" --top_k 5
```

**示例3：检索案例示例**
```bash
python tools/retriever.py --query "应力计算 示例 INPUT STRU" --top_k 5
```

**示例4：安静模式（适合脚本调用）**
```bash
python tools/retriever.py --query "DFT计算" --top_k 3 --quiet
```

### 输出格式

**正常模式：**
```
================================================================================
检索查询: DFT计算
返回文档数: 3
================================================================================

文档1:
来源: xxx.md
内容:
[文档内容]

--------------------------------------------------------------------------------

文档2:
来源: yyy.md
内容:
[文档内容]

--------------------------------------------------------------------------------
```

**安静模式：**
```
[来源: xxx.md]
[文档内容]

[来源: yyy.md]
[文档内容]
```

### 前置要求

1. **知识库已初始化**
   - 知识库路径：`data/chroma_db/`
   - 如未初始化，运行：`python src/knowledge_manager.py`

2. **API Key已配置**
   - 在 `config.yaml` 中配置 `api_key`（DashScope API Key）
   - 或设置环境变量：`DASHSCOPE_API_KEY`

### 故障排查

**错误：无法连接到知识库**
```
[错误] 无法连接到知识库
[提示] 请确保知识库已初始化，路径: data/chroma_db/
```
解决方法：运行 `python src/knowledge_manager.py` 初始化知识库

**错误：未找到API Key**
```
[错误] 未找到DashScope API Key
[提示] 请在config.yaml中配置api_key
```
解决方法：在 `config.yaml` 中添加：
```yaml
llm:
  api_key: "your-dashscope-api-key"
```

**错误：编码问题**
如果在Windows上遇到编码问题，工具会自动处理，使用UTF-8编码输出。

---

## case_parser.py - 案例解析工具

### 功能
解析用户提供的案例文件（docx/md/txt），自动提取关键信息：
- 文件结构（INPUT、STRU、KPT等）
- 关键参数及其值
- 计算流程
- 特殊设置和注意事项

### 使用方法

**基本用法：**
```bash
python tools/case_parser.py --input "path/to/case.docx"
```

**参数说明：**
- `--input`：案例文件路径（必需）
- 支持格式：`.docx`、`.md`、`.txt`

### 使用示例

**示例1：解析docx案例**
```bash
python tools/case_parser.py --input "data/input/example.docx"
```

**示例2：解析markdown案例**
```bash
python tools/case_parser.py --input "data/input/test_case.md"
```

### 输出格式

```
================================================================================
案例解析结果
================================================================================

文件: data/input/test_case.md

## 文件结构
- INPUT
- STRU
- KPT
- 赝势文件 (1个)

## 关键参数
calculation = scf
basis_type = pw
ecutwfc = 80
ecutrho = 320
scf_thr = 1e-8
cal_stress = 1
stress_thr = 0.01
smearing_method = gaussian
mixing_type = pulay
ks_solver = genelpa

## 计算流程
1. 结构优化
2. 自洽计算
3. 应力计算

## 特殊设置
- 高精度
- 自定义阈值
- 说明: 本案例使用了高精度设置
```

### 解析能力

**文件结构识别：**
- 自动识别常见ABACUS文件：INPUT、STRU、KPT、CONTROL等
- 识别赝势文件（.upf）和轨道文件（.orb）
- 统计文件数量

**参数提取：**
自动提取30+个常见ABACUS参数，包括：
- 计算类型：calculation, basis_type
- 截断能：ecutwfc, ecutrho
- 收敛参数：scf_thr, scf_nmax
- 应力/力：cal_stress, cal_force, stress_thr, force_thr
- 弛豫：relax_nmax
- 展宽：smearing_method, smearing_sigma
- 混合：mixing_type, mixing_beta
- 求解器：ks_solver
- 其他：nbands, nspin, symmetry等

**计算流程识别：**
自动识别常见计算类型：
- 结构优化（relax, optimization）
- 自洽计算（scf）
- 非自洽计算（nscf）
- 能带计算（band）
- 态密度计算（dos）
- 应力/力计算（stress, force）
- 分子动力学（md）

**特殊设置识别：**
- 高精度设置
- 自定义阈值
- 特殊k点设置
- 自旋极化
- 范德华修正
- DFT+U
- 提取注释中的说明

### 前置要求

**依赖库：**
- `python-docx`（仅解析docx文件时需要）

安装：
```bash
pip install python-docx
```

**注意：** 如果只解析md/txt文件，不需要安装python-docx

### 故障排查

**错误：文件不存在**
```
[错误] 文件不存在: path/to/case.docx
```
解决方法：检查文件路径是否正确

**错误：不支持的文件格式**
```
[错误] 不支持的文件格式: .pdf
[提示] 支持的格式: .docx, .md, .txt
```
解决方法：将案例转换为支持的格式

**错误：缺少python-docx库**
```
[错误] 缺少python-docx库
[提示] 请安装: pip install python-docx
```
解决方法：安装python-docx库

### 解析结果说明

**如果某项未识别：**
- 文件结构：显示"未识别到文件结构"
- 参数：显示"未识别到参数，可能需要手动提取"
- 计算流程：显示"未识别到明确的计算流程"
- 特殊设置：显示"未识别到特殊设置"

这是正常的，说明案例文件中可能没有明确提到这些信息，需要手动补充。

---

## 在Claude Code中使用

这些工具设计用于在Claude Code执行流程时调用。

**在CLAUDE.md中的使用示例：**

```markdown
### Step 1: 知识调研

使用Bash执行：
```bash
python tools/retriever.py --query "主题 物理原理 参数设置" --top_k 10
```

读取输出结果，保存到 01_research.md
```

**Claude Code会：**
1. 执行Bash命令
2. 读取终端输出
3. 分析检索结果质量
4. 保存到对应文件
5. Think Aloud说明检索情况
