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

---

## orbital_validator.py - 轨道文件名验证工具

### 功能
扫描教程 Markdown 文件中的轨道文件名，验证并自动修正两类问题：
1. **轨道文件名错误**：如 `Si_gga_7au_100Ry_2s2p1d.orb` → `Si_gga_8au_100Ry_2s2p1d.orb`
2. **INPUT 参数不兼容**：如 `nbands auto` → 删除（ABACUS v3.10.x 不支持）

### 使用方法

**只检查，不修改：**
```bash
python tools/orbital_validator.py path/to/tutorial.md
```

**检查并自动修正（覆盖原文件）：**
```bash
python tools/orbital_validator.py path/to/tutorial.md --fix
```

**检查并输出到新文件：**
```bash
python tools/orbital_validator.py path/to/tutorial.md --fix --output path/to/fixed.md
```

### 使用场景
在教程生成的 Step 7（最终输出）自动调用，确保轨道文件名和参数正确。

---

## fix_stru.py - STRU 文件格式修复工具

### 功能
修复 STRU 文件中的两类问题：
1. 将赝势和轨道文件路径改为相对路径（解决 Bohrium 容器环境变量问题）
2. 将旧版 4 字段格式转为 ABACUS v3.x 新格式（3 字段 + NUMERICAL_ORBITAL 块）

### 使用方法

```bash
python tools/fix_stru.py --stru path/to/STRU
```

### 使用场景
在计算测试的 Step 2（输入文件验证）自动调用。

---

## static_checker.py - 静态质量检查工具

### 功能
检查 7 项可客观判定的教程质量指标（用于 judgeCLAUDE 评分系统）：
1. 信息精准与一致性（max 5）
2. 文件来源明确度（max 3）
3. 时效与兼容性（max 4）
4. 可导航性与步骤连贯性（max 1）
5. 结构完整性（max 5）
6. 环境条件明确度（max 3）
7. 图文/代码邻近性（max 3）

### 使用方法

```bash
python tools/static_checker.py --input tutorial.md --output scores.json
```

### 输出格式
生成 JSON 文件，包含各项得分和扣分原因。

---

## test_framework_integrated.py - 计算测试框架

### 功能
整合所有测试阶段，支持插件化测试。包含 14 个计算类型插件。

### 使用方法

**Phase 1：准备输入文件**
```bash
python tools/test_framework_integrated.py tutorial.md --test-dir ./test_dir --phase prepare
```

**继续后处理（任务完成后）**
```bash
python tools/test_framework_integrated.py continue ./test_dir
```

### 支持的计算类型
- relax（结构优化）
- elastic（弹性常数）
- band（能带结构）
- dos（态密度）
- dftu（DFT+U 强关联）
- optic（光学性质）
- solvation（隐式溶剂）
- phonopy（声子谱）
- sdft（SDFT/MDFT）
- elf（电子局域函数）
- tddft（RT-TDDFT）
- neb（NEB 过渡态）
- scf（SCF 自洽计算）

### 插件开发
参考 `testCLAUDE/plugin_dev_guide.md` 和 `test_plugins/PLUGIN_REGISTRY.md`。

---

## testCLAUDE/ 模块文档

### bohrium_setup.md
Bohrium 配置流程 + CLI 命令速查表。

### plugin_dev_guide.md
插件自主创建流程（testCLAUDE.md Step 1.5）。

### plugins_history.md
插件开发历史表，记录所有插件的添加时间和首次测试教程。

### troubleshooting.md
8 类故障排除清单（Step 4.4/6.3 出错时加载）。

### file_formats.md
job.json/analysis.json/STRU 格式参考（Step 2 时加载）。

---

## 工具依赖关系

```
教程生成流程（CLAUDE.md）
├── retriever.py          # Step 1, 3
├── case_parser.py        # Step 1
└── orbital_validator.py  # Step 7

计算测试流程（testCLAUDE.md）
├── test_framework_integrated.py  # 主框架
│   ├── test_framework_phase1_analyzer.py
│   └── test_framework_phase3_7_impl.py
├── fix_stru.py           # Step 2
├── orbital_validator.py  # Step 2
└── test_plugins/         # 14 个插件
    ├── base_plugin.py
    ├── relax_plugin.py
    ├── elastic_plugin.py
    ├── band_plugin.py
    ├── dos_plugin.py
    ├── dftu_plugin.py
    ├── optic_plugin.py
    ├── solvation_plugin.py
    ├── phonopy_plugin.py
    ├── sdft_plugin.py
    ├── elf_plugin.py
    ├── tddft_plugin.py
    ├── neb_plugin.py
    └── scf_plugin.py

质量评估
└── static_checker.py     # judgeCLAUDE 评分
```
