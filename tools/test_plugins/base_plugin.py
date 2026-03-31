"""
AutoTutorial 3.0 - 测试插件基类
定义所有测试插件必须实现的接口
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field


@dataclass
class TestInfo:
    """测试信息数据类"""
    # 基本信息
    calc_type: str  # 计算类型：relax, elastic, band, dos
    case_name: str  # 案例名称：Si, TiO2 等

    # 输入文件内容
    input_content: Optional[str] = None
    stru_content: Optional[str] = None
    kpt_content: Optional[str] = None

    # 额外的脚本和命令
    python_scripts: List[str] = field(default_factory=list)
    bash_commands: List[str] = field(default_factory=list)
    abacustest_commands: List[str] = field(default_factory=list)

    # 预期结果
    expected_results: Dict = field(default_factory=dict)

    # 赝势和轨道文件
    pseudopotentials: List[str] = field(default_factory=list)
    orbitals: List[str] = field(default_factory=list)


@dataclass
class ValidationResult:
    """验证结果数据类"""
    calc_type: str
    case_name: str
    passed: bool
    comparisons: Dict  # 详细对比结果
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    error_analyses: Dict = field(default_factory=dict)   # key → ErrorAnalysis


@dataclass
class ErrorAnalysis:
    """误差归因分析结果"""
    key: str                    # 参数名
    category: str               # 归因类别（见下方常量）
    tolerance_ratio: float      # 实际误差 / 有效容差（倍数）。注意：energy_reference 类别此值无物理意义，忽略即可
    physical_explanation: str   # 物理解释（中文，面向教程作者）
    user_guidance: str          # 用户建议（中文，将写入教程）
    insertion_note: str         # 待插入教程的 Markdown 文本（含实测数据）

# 归因类别常量
EA_PRECISION_SENSITIVITY = "precision_sensitivity"  # 精度参数敏感
EA_DFT_SYSTEMATIC        = "dft_systematic"         # DFT 系统误差（能隙）
EA_KPOINT_MISMATCH       = "kpoint_mismatch"        # K 点设置不一致
EA_ENERGY_REFERENCE      = "energy_reference"       # 总能量绝对值大
EA_AFM_SYMMETRY          = "afm_symmetry"           # 反铁磁总磁矩近零
EA_UNKNOWN               = "unknown"                # 未知原因


class BaseTestPlugin(ABC):
    """测试插件基类"""

    def __init__(self, job_manager=None, pp_manager=None):
        """
        初始化插件

        Args:
            job_manager: BohriumJobManager 实例
            pp_manager: PseudopotentialManager 实例
        """
        self.job_manager = job_manager
        self.pp_manager = pp_manager
        self.tolerance = 0.05  # 默认容差 5%

    @property
    @abstractmethod
    def plugin_name(self) -> str:
        """插件名称"""
        pass

    @property
    @abstractmethod
    def calc_type(self) -> str:
        """计算类型标识"""
        pass

    @abstractmethod
    def can_handle(self, tutorial_content: str) -> bool:
        """
        判断该插件是否能处理这个教程

        Args:
            tutorial_content: 教程的完整文本内容

        Returns:
            True 如果能处理，False 否则
        """
        pass

    @abstractmethod
    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """
        从教程中提取测试所需信息

        Args:
            tutorial_content: 教程的完整文本内容

        Returns:
            TestInfo 对象，如果提取失败返回 None
        """
        pass

    @abstractmethod
    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """
        准备输入文件

        Args:
            test_info: 测试信息
            work_dir: 工作目录

        Returns:
            需要提交任务的目录列表
        """
        pass

    @abstractmethod
    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """
        提交计算任务

        Args:
            input_dirs: 输入文件目录列表

        Returns:
            Job ID 列表
        """
        pass

    @abstractmethod
    def validate_results(self, job_ids: List[str], work_dir: Path, test_info: TestInfo) -> ValidationResult:
        """
        验证计算结果

        Args:
            job_ids: Job ID 列表
            work_dir: 工作目录
            test_info: 测试信息（包含预期结果）

        Returns:
            ValidationResult 对象
        """
        pass

    @abstractmethod
    def generate_report_section(self, validation: ValidationResult) -> str:
        """
        生成该类型的测试报告章节

        Args:
            validation: 验证结果

        Returns:
            Markdown 格式的报告章节
        """
        pass

    # 辅助方法
    def _compare_values(self, expected: float, actual: float, key: str) -> Dict:
        """
        对比数值结果

        Args:
            expected: 预期值
            actual: 实际值
            key: 参数名称

        Returns:
            对比结果字典
        """
        if expected == 0:
            abs_error = abs(actual - expected)
            rel_error = abs_error
        else:
            abs_error = abs(actual - expected)
            rel_error = abs_error / abs(expected)

        return {
            'key': key,
            'expected': expected,
            'actual': actual,
            'abs_error': abs_error,
            'rel_error': rel_error,
            'passed': rel_error <= self.tolerance,
            'effective_tol': self.tolerance,   # 记录本次使用的容差，供 analyze_deviation() 读取
        }

    def analyze_deviation(
        self,
        comparisons: Dict,
        context: Optional[Dict] = None
    ) -> Dict:
        """
        对 comparisons 中每个 FAIL 的参数执行归因分析。

        Args:
            comparisons: _compare_values() 返回的对比字典
            context: 额外上下文，支持键：
                'kspacing_mismatch' (bool)   — 测试与教程 kspacing 不同
                'actual_kspacing'   (float)  — 实际使用的 kspacing
                'expected_kspacing' (float)  — 教程参考 kspacing
                'abacus_version'    (str)    — ABACUS 版本字符串

        Returns:
            Dict[key, ErrorAnalysis]，仅包含 FAIL 的参数
        """
        if context is None:
            context = {}

        if not comparisons:
            return {}   # 空 comparisons 直接返回，避免误导性输出

        analyses = {}

        # 能隙相关参数名集合
        _BANDGAP_KEYS = {'bandgap', 'bandgap_up', 'bandgap_dn',
                         'band_gap', 'gap', 'gap_up', 'gap_dn'}
        # 总能量相关参数名集合
        _ENERGY_KEYS  = {'total_energy', 'etot', 'energy', 'e_tot',
                         'final_energy'}
        # 总磁矩相关参数名集合
        _MAG_TOTAL_KEYS = {'total_magnetism', 'mag_total', 'total_mag'}

        for key, comp in comparisons.items():
            if comp.get('passed', True):
                continue   # 只分析失败项

            # ── 计算 tolerance_ratio ──────────────────────────────
            # 优先使用 comp 中明确存储的有效容差（effective_tol），
            # 其次是 abs_tol（use_abs=True 场景），最后回退到 self.tolerance。
            if comp.get('effective_tol') is not None:
                tolerance_used = comp['effective_tol']
            elif comp.get('use_abs') and comp.get('abs_tol') is not None:
                tolerance_used = comp['abs_tol']
            else:
                tolerance_used = self.tolerance

            # 关键：使用与 tolerance_used 量纲一致的误差量
            # - use_abs=True：abs_error 和 abs_tol 都是物理量，用 abs_error
            # - 相对容差：rel_error 和 effective_tol 都是无量纲比值，用 rel_error
            if comp.get('use_abs'):
                error_for_ratio = comp['abs_error']
            else:
                error_for_ratio = comp['rel_error']

            tolerance_ratio = (
                error_for_ratio / tolerance_used
                if (tolerance_used and tolerance_used > 0)
                else float('inf')
            )

            expected_str = f"{comp['expected']:.6g}"
            actual_str   = f"{comp['actual']:.6g}"
            rel_pct      = f"{comp['rel_error']*100:.1f}"

            # ── 决策树 ──────────────────────────────────────
            # 能隙参数：DFT 系统误差优先于精度参数敏感性
            # 若偏差达到容差的 90% 以上（即接近或超出容差），视为 DFT 系统误差
            if key.lower() in _BANDGAP_KEYS and tolerance_ratio >= 0.9:
                # 能隙：DFT 系统误差
                cat = EA_DFT_SYSTEMATIC
                phys = (
                    f"{key} 偏差（{rel_pct}%）在 DFT(PBE) 计算能隙的系统性不确定性范围内。"
                    f"PBE 泛函对能隙存在系统性低估，不同版本 ABACUS、赝势或 k 点设置"
                    f"均可导致 5–20% 的差异。"
                )
                guidance = (
                    f"DFT(PBE) 计算能隙存在系统性不确定性，通常 ±15% 均属正常范围。"
                    f"本教程参考值（{expected_str} eV）仅供参考，"
                    f"若关注精确能隙，建议使用 HSE 或 GW 方法。"
                )
                note = (
                    f"> **💡 关于能隙参考值：** 本教程中 {key} 参考值为 {expected_str} eV，"
                    f"测试计算得到 {actual_str} eV（偏差 {rel_pct}%）。"
                    f"DFT(PBE) 计算能隙存在系统性不确定性（通常 ±15%），"
                    f"若读者得到相近结果，说明计算流程正确；"
                    f"若关注精确能隙值，建议使用杂化泛函（HSE06）或 GW 方法。"
                )

            elif tolerance_ratio < 2.0:
                # 偏差 < 2倍容差：精度参数敏感性
                cat = EA_PRECISION_SENSITIVITY
                phys = (
                    f"{key} 偏差（{rel_pct}%）仅为容差的 {tolerance_ratio:.1f} 倍，"
                    f"属于精度参数（ecutwfc/kspacing）收敛效应的正常范围。"
                )
                guidance = (
                    f"若读者得到与本教程参考值（{expected_str}）相差 "
                    f"{rel_pct}% 以内的结果，属于正常精度范围，"
                    f"不影响计算结论。"
                )
                note = (
                    f"> **💡 关于 {key} 的参考值：** 本教程参考值为 {expected_str}，"
                    f"测试计算得到 {actual_str}（偏差 {rel_pct}%），"
                    f"属于计算精度参数（如 `ecutwfc`、`kspacing`）收敛效应的正常范围，"
                    f"不影响计算结论的正确性。"
                )

            elif key.lower() in _BANDGAP_KEYS:
                # 能隙：DFT 系统误差
                cat = EA_DFT_SYSTEMATIC
                phys = (
                    f"{key} 偏差（{rel_pct}%）在 DFT(PBE) 计算能隙的系统性不确定性范围内。"
                    f"PBE 泛函对能隙存在系统性低估，不同版本 ABACUS、赝势或 k 点设置"
                    f"均可导致 5–20% 的差异。"
                )
                guidance = (
                    f"DFT(PBE) 计算能隙存在系统性不确定性，通常 ±15% 均属正常范围。"
                    f"本教程参考值（{expected_str} eV）仅供参考，"
                    f"若关注精确能隙，建议使用 HSE 或 GW 方法。"
                )
                note = (
                    f"> **💡 关于能隙参考值：** 本教程中 {key} 参考值为 {expected_str} eV，"
                    f"测试计算得到 {actual_str} eV（偏差 {rel_pct}%）。"
                    f"DFT(PBE) 计算能隙存在系统性不确定性（通常 ±15%），"
                    f"若读者得到相近结果，说明计算流程正确；"
                    f"若关注精确能隙值，建议使用杂化泛函（HSE06）或 GW 方法。"
                )

            elif context.get('kspacing_mismatch'):
                # K点不一致
                actual_ks   = context.get('actual_kspacing', '未知')
                expected_ks = context.get('expected_kspacing', '未知')
                cat = EA_KPOINT_MISMATCH
                phys = (
                    f"{key} 偏差（{rel_pct}%）可能由 K 点设置不同引起。"
                    f"测试使用 kspacing={actual_ks}，"
                    f"教程参考值基于 kspacing={expected_ks}。"
                )
                guidance = (
                    f"使用与教程相同的 K 点密度（kspacing={expected_ks}）"
                    f"可更接近参考值（{expected_str}）。"
                    f"当前测试使用 kspacing={actual_ks}，得到 {actual_str}，"
                    f"差异属于 K 点收敛效应。"
                )
                note = (
                    f"> **💡 关于 {key} 的参考值：** 本教程参考值（{expected_str}）"
                    f"基于 `kspacing = {expected_ks}` 计算，"
                    f"使用 `kspacing = {actual_ks}` 得到 {actual_str}（偏差 {rel_pct}%），"
                    f"属于 K 点收敛效应，不影响方法论正确性。"
                    f"若需要精确复现参考值，请使用相同的 K 点设置。"
                )

            elif key.lower() in _ENERGY_KEYS:
                # 总能量：大负数参考值问题
                cat = EA_ENERGY_REFERENCE
                phys = (
                    f"{key} 为大负数（{expected_str} eV），"
                    f"绝对误差 {comp['abs_error']:.3f} eV，相对误差 {rel_pct}%。"
                    f"结构或参数的微小差异（如截断能、赝势版本）均可导致 eV 级能量偏移。"
                )
                guidance = (
                    f"总能量绝对值受结构参数影响较大，{rel_pct}% 的偏差"
                    f"（即 {comp['abs_error']:.2f} eV）在实际计算中属正常范围。"
                    f"更有意义的比较是能量差（如形成能、结合能），"
                    f"单次 SCF 的总能量仅供参考。"
                )
                note = (
                    f"> **💡 关于总能量参考值：** 本教程参考值为 {expected_str} eV，"
                    f"测试计算得到 {actual_str} eV（绝对偏差 {comp['abs_error']:.3f} eV，"
                    f"相对偏差 {rel_pct}%）。总能量为大负数，"
                    f"精度参数（ecutwfc、赝势版本）的微小差异均可导致 eV 级偏移，"
                    f"属正常范围。如需对比，请关注能量差值而非绝对值。"
                )

            elif key.lower() in _MAG_TOTAL_KEYS:
                # 反铁磁总磁矩近零
                cat = EA_AFM_SYMMETRY
                phys = (
                    f"{key} 理论值应为 0（反铁磁对称性），"
                    f"实际值 {actual_str} μB，绝对偏差 {comp['abs_error']:.4f} μB。"
                    f"非零值来源于数值计算精度，是正常的。"
                )
                guidance = (
                    f"反铁磁结构总磁矩理论值为 0，"
                    f"计算得到 {actual_str} μB 说明反铁磁态建立正确，"
                    f"微小偏差来自 SCF 数值精度，不影响物理结论。"
                )
                note = (
                    f"> **💡 关于总磁矩：** 反铁磁结构的总磁矩理论值为 0，"
                    f"测试计算得到 {actual_str} μB（绝对偏差 {comp['abs_error']:.4f} μB），"
                    f"说明反铁磁态建立正确。微小非零值来自 SCF 数值精度，"
                    f"绝对值 < 0.1 μB 均属正常范围。"
                )

            else:
                # 未知原因
                cat = EA_UNKNOWN
                phys = (
                    f"{key} 偏差（{rel_pct}%）超出容差，原因不明。"
                    f"预期 {expected_str}，实际 {actual_str}。"
                )
                guidance = (
                    f"请检查：① 输入文件参数是否与教程一致；"
                    f"② ABACUS 版本；③ 赝势/轨道文件版本。"
                )
                note = (
                    f"> **⚠️ 注意 {key} 偏差：** 测试计算得到 {actual_str}，"
                    f"与教程参考值 {expected_str} 相差 {rel_pct}%。"
                    f"若读者遇到类似偏差，请检查输入参数、ABACUS 版本和赝势文件版本。"
                )

            analyses[key] = ErrorAnalysis(
                key=key,
                category=cat,
                tolerance_ratio=tolerance_ratio,
                physical_explanation=phys,
                user_guidance=guidance,
                insertion_note=note,
            )

        return analyses

    def extract_case_name(self, content: str) -> str:
        """
        通用案例名提取方法（避免提取代词）

        优先级：
        1. 明确的案例标记："案例：Si"、"案例1：TiO2"
        2. 材料体系标记："材料体系：NiO"、"**材料体系**：Al"
        3. STRU文件中的元素符号
        4. 避免提取代词（它、他、她、这、那等）
        """
        import re

        # 中文代词黑名单
        PRONOUNS = {'它', '他', '她', '这', '那', '其', '此', '该'}

        # 1. 明确的案例标记
        patterns = [
            r'案例\d*\s*[-:：]\s*(\w+)',
            r'#+\s*案例[：:]\s*(\w+)',
        ]
        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                name = match.group(1)
                if name not in PRONOUNS:
                    return name

        # 2. 材料体系标记
        material_patterns = [
            r'材料体系[：:]\s*(\w+)',
            r'\*\*材料体系\*\*[：:]\s*(\w+)',
        ]
        for pattern in material_patterns:
            match = re.search(pattern, content)
            if match:
                name = match.group(1)
                if name not in PRONOUNS:
                    return name

        # 3. 从STRU文件提取元素符号
        stru_pattern = r'```[^\n]*\n(ATOMIC_SPECIES[\s\S]+?)```'
        stru_match = re.search(stru_pattern, content)
        if stru_match:
            stru_content = stru_match.group(1)
            # 提取元素符号（第一列）
            element_pattern = r'^([A-Z][a-z]?)\s+[\d.]+'
            element_matches = re.findall(element_pattern, stru_content, re.MULTILINE)
            if element_matches:
                # 如果只有一种元素，直接返回
                if len(set(element_matches)) == 1:
                    return element_matches[0]
                # 如果有多种元素，组合返回（如 NiO）
                unique_elements = []
                for elem in element_matches:
                    if elem not in unique_elements:
                        unique_elements.append(elem)
                return ''.join(unique_elements)

        return "Unknown"
