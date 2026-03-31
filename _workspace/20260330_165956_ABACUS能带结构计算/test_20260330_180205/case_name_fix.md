# 案例名提取逻辑改进说明

## 问题描述

测试框架在提取案例名时，使用了正则表达式 `r'计算(\w+)的'`，会匹配"计算它的能带结构"这样的句子，错误地提取出代词"它"作为案例名。

**问题表现：**
- `01_analysis.json` 中 `case_names` 为 `["它", "它", "它", "Si"]`
- 生成的目录名为 `它_scf`、`它_nscf` 等

## 解决方案

### 1. 在基类中添加通用方法

在 `tools/test_plugins/base_plugin.py` 中添加 `extract_case_name()` 方法：

**改进要点：**
- 添加中文代词黑名单：`{'它', '他', '她', '这', '那', '其', '此', '该'}`
- 优先匹配明确的案例标记："案例：Si"、"材料体系：NiO"
- 从 STRU 文件的 ATOMIC_SPECIES 块提取元素符号
- 对提取结果进行代词过滤

**提取优先级：**
1. 明确的案例标记（`案例：Si`、`案例1：TiO2`）
2. 材料体系标记（`材料体系：NiO`、`**材料体系**：Al`）
3. STRU 文件中的元素符号（单元素返回元素名，多元素组合返回如 NiO）

### 2. 更新所有插件

将以下插件的 `_extract_case_name()` 方法改为调用基类方法：
- `band_plugin.py`
- `dos_plugin.py`
- `elastic_plugin.py`
- `relax_plugin.py`
- `scf_plugin.py`

**修改示例：**
```python
def _extract_case_name(self, content: str) -> str:
    """提取案例名称（使用基类通用方法）"""
    return self.extract_case_name(content)
```

## 验证结果

重新运行 prepare 阶段后：

**修复前：**
```json
{
  "case_names": ["它", "它", "它", "Si"]
}
```

**修复后：**
```json
{
  "case_names": ["Si", "Si", "Si", "Si"]
}
```

**目录命名：**
- 修复前：`它_scf`、`它_nscf`、`它_dos`
- 修复后：`Si_scf`、`Si_nscf`、`Si_dos`

## 改进效果

1. **避免提取代词**：黑名单机制确保不会提取"它"、"这"等代词
2. **提高准确性**：从 STRU 文件提取元素符号，更可靠
3. **统一逻辑**：所有插件使用相同的提取方法，便于维护
4. **向后兼容**：保留原有的案例标记匹配模式

## 测试案例

已在以下教程中验证：
- ✅ ABACUS能带结构计算（Si）- 从"计算它的能带结构"正确提取为"Si"
- ✅ 其他教程待测试

## 相关文件

- `tools/test_plugins/base_plugin.py` - 新增 `extract_case_name()` 方法
- `tools/test_plugins/band_plugin.py` - 更新
- `tools/test_plugins/dos_plugin.py` - 更新
- `tools/test_plugins/elastic_plugin.py` - 更新
- `tools/test_plugins/relax_plugin.py` - 更新
- `tools/test_plugins/scf_plugin.py` - 更新

---

**修复日期：** 2026-03-31
**修复人员：** Claude (AutoTutorial 3.0)
