# 新测试插件自主创建指南

> 本文件由 testCLAUDE.md Step 1.4 按需 read，仅在 `detected_types = []` 时加载。
> 完成后返回主流程继续 Step 2。

---

### Step 1.5：创建新测试插件（⭐ 自主扩展框架）

**目标：** 为新的计算类型创建标准测试插件，使 `analysis.json` 能正确输出。

**本步骤是自主执行步骤，无需等待用户确认。**

---

**1.5.1 收集插件所需信息**

在动手写代码前，必须明确以下 6 项：

```
1. calc_type（计算类型标识）：全小写，如 dftu / tddft / md / hse
2. can_handle 关键词：教程 INPUT 代码块中的唯一特征（如 dft_plus_u = 1）
3. INPUT 内容：从教程代码块提取，还是需要构造模板？
4. STRU 内容：教程是否给出完整 STRU（含晶格矢量）？
   → 若无，查 RAG 数据库 / ABACUS GitHub tests / Bohrium 数据集
5. KPT：教程是否给出？若无，使用 4×4×4 Gamma
6. expected_results：从教程中提取所有给出的数值
```

**STRU 不完整时的补全策略：**

| 情况 | 推荐方法 |
|------|----------|
| 教程给出体系和晶格参数 | 从 RAG 检索 POSCAR/STRU，验证键长后使用 |
| 教程提到数据集路径（如 Bohrium URL）| 从数据集文档或 ABACUS GitHub 对应测试找 STRU |
| 教程有 Python 生成脚本 | 执行脚本生成 CIF → abacustest 转换为 STRU |
| 以上均失败 | 用 ASE 自动生成，并在报告中注明"近似结构" |

**验证容差选择原则：**

| 物理量 | 推荐容差 | 理由 |
|--------|----------|------|
| 总能量（大负数如 -9255 eV）| 相对 0.1% | 绝对差几 eV 物理上微小，相对容差更合理 |
| 能隙 | 相对 15% | 晶格矢量微小差异会导致能隙~10% 变化 |
| 总磁矩（应为 0）| 绝对 0.1 | 直接判断是否近零 |
| 绝对磁矩 | 相对 5% | 结构相近时磁矩稳定 |
| 弹性常数 | 相对 5% | 标准弹性计算精度 |

---

**1.5.2 参照已有插件写代码**

以 `tools/test_plugins/relax_plugin.py` 或 `tools/test_plugins/dftu_plugin.py` 为模板，创建新文件 `tools/test_plugins/<type>_plugin.py`。

**必须实现的 6 个方法：**

```python
class XxxPlugin(BaseTestPlugin):

    @property
    def plugin_name(self) -> str: ...     # 中文名，如 "DFT+U 强关联体系测试"

    @property
    def calc_type(self) -> str: ...       # 如 "dftu"

    def can_handle(self, tutorial_content) -> bool: ...   # 关键词检测

    def extract_test_info(self, tutorial_content) -> TestInfo: ...
    # 提取 INPUT/STRU/KPT，设置 expected_results / pseudopotentials / orbitals

    def prepare_inputs(self, test_info, work_dir) -> List[Path]: ...
    # 写文件，下载赝势/轨道，修复路径（_fix_input_file / fix_stru_paths）

    def submit_jobs(self, input_dirs) -> List[str]: ...
    # 通常直接调用 self.job_manager.create_job_config + submit_job

    def validate_results(self, job_ids, work_dir, test_info) -> ValidationResult: ...
    # 从 running_*.log 提取实际值，逐项对比

    def generate_report_section(self, validation) -> str: ...
    # Markdown 表格，同 RelaxPlugin 格式
```

---

**1.5.3 注册新插件**

在 `tools/test_framework_integrated.py` 中：

```python
# 1. 添加 import
from test_plugins.<type>_plugin import <Type>Plugin

# 2. 在 self.plugins 列表末尾添加
self.plugins: List[BaseTestPlugin] = [
    RelaxPlugin(...),
    ElasticPlugin(...),
    BandPlugin(...),
    DOSPlugin(...),
    DFTUPlugin(...),
    <Type>Plugin(self.job_manager, self.pp_manager),  # ← 新增
]
```

---

**1.5.4 验证新插件**

重新运行 prepare：
```bash
python tools/test_framework_integrated.py "$tutorial_path" --test-dir "$test_dir" --phase prepare
```

预期输出：
```
[OK] 检测到: <新插件名称>
[OK] 提取测试信息: 案例 <体系名>
[OK] 分析结果已保存: 01_analysis.json
```

验证 `01_analysis.json` 不再为空后，继续 Step 2。

---

**Think Aloud（贯穿整个 Step 1.5）：**
- 说明判断的计算类型和依据
- 说明如何处理 STRU 不完整的情况
- 说明选择容差的理由
- 说明插件注册后的验证结果

**完成标志：**
- 新插件文件已创建
- `test_framework_integrated.py` 已更新
- 重跑 prepare 后 `analysis.json` 非空
- 在 `tools/testCLAUDE/plugins_history.md` 末尾追加一行（含插件文件名、计算类型、添加日期、首次测试教程）

**写入 issues_log.json：**

```bash
python -c "
import json, os
log_path = '$test_dir/issues_log.json'
log = json.loads(open(log_path).read()) if os.path.exists(log_path) else {'issues': []}
log['issues'].append({
    'type': 'plugin_created',
    'category': 'user-confirm',
    'description': '教程计算类型无对应插件，已临时创建 <新插件名>',
    'resolution': '已创建并注册新插件',
    'tutorial_keywords': ['<calc_type关键词>'],
    'insertion_note': '> **📌 说明：** 本教程涉及 <计算类型> 功能，为较新特性，使用前请确认 ABACUS 版本 ≥ <最低版本>。',
    'step': 'Step 1.5'
})
open(log_path, 'w', encoding='utf-8').write(json.dumps(log, ensure_ascii=False, indent=2))
print('[记录] 已写入 issues_log.json：创建新插件 <新插件名>')
"
```
（将 `<新插件名>`、`<calc_type关键词>`、`<计算类型>`、`<最低版本>` 替换为本次实际值）
