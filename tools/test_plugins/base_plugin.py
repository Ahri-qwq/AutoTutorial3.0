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
            'passed': rel_error <= self.tolerance
        }
