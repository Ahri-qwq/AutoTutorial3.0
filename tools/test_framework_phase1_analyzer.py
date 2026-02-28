"""
AutoTutorial 3.0 - 测试脚本框架
Phase 1: 文章分析模块(从Markdown教程中提取计算案例信息)
"""

import re
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class TutorialAnalysis:
    """文章分析结果"""
    # 基本信息
    title: str
    topic: str

    # 依赖分析
    needs_ase: bool = False
    needs_abacustest: bool = False
    needs_other_tools: List[str] = None

    # 文件分析
    has_complete_input: bool = False
    has_complete_stru: bool = False
    has_complete_kpt: bool = False

    # 提取的内容
    input_content: Optional[str] = None
    stru_content: Optional[str] = None
    kpt_content: Optional[str] = None

    # 命令序列
    python_scripts: List[str] = None
    abacustest_commands: List[str] = None
    bash_commands: List[str] = None

    # 预期结果
    expected_results: Dict = None


class TutorialAnalyzer:
    """文章分析器"""

    def __init__(self, tutorial_path: str):
        self.tutorial_path = Path(tutorial_path)
        self.content = self._read_tutorial()

    def _read_tutorial(self) -> str:
        """读取教程文件"""
        with open(self.tutorial_path, 'r', encoding='utf-8') as f:
            return f.read()

    def analyze(self) -> TutorialAnalysis:
        """分析教程"""
        analysis = TutorialAnalysis(
            title=self._extract_title(),
            topic=self._extract_topic(),
            python_scripts=[],
            abacustest_commands=[],
            bash_commands=[],
            expected_results={}
        )

        # 1. 检查工具依赖
        analysis.needs_ase = self._check_ase_dependency()
        analysis.needs_abacustest = self._check_abacustest_dependency()
        analysis.needs_other_tools = self._check_other_tools()

        # 2. 检查文件完整性
        analysis.has_complete_input = self._check_input_file()
        analysis.has_complete_stru = self._check_stru_file()
        analysis.has_complete_kpt = self._check_kpt_file()

        # 3. 提取内容
        if analysis.has_complete_input:
            analysis.input_content = self._extract_input()
        if analysis.has_complete_stru:
            analysis.stru_content = self._extract_stru()
        if analysis.has_complete_kpt:
            analysis.kpt_content = self._extract_kpt()

        # 4. 提取命令
        analysis.python_scripts = self._extract_python_scripts()
        analysis.abacustest_commands = self._extract_abacustest_commands()
        analysis.bash_commands = self._extract_bash_commands()

        # 5. 提取预期结果
        analysis.expected_results = self._extract_expected_results()

        return analysis

    def _extract_title(self) -> str:
        """提取标题"""
        match = re.search(r'^#\s+(.+)$', self.content, re.MULTILINE)
        return match.group(1) if match else "Unknown"

    def _extract_topic(self) -> str:
        """提取主题"""
        # 从文件名或标题推断
        return self.tutorial_path.stem

    def _check_ase_dependency(self) -> bool:
        """检查是否需要ASE"""
        patterns = [
            r'from ase',
            r'import ase',
            r'ase\.build',
            r'bulk\(',
        ]
        return any(re.search(p, self.content) for p in patterns)

    def _check_abacustest_dependency(self) -> bool:
        """检查是否需要abacustest"""
        patterns = [
            r'abacustest',
            r'model inputs',
            r'model elastic',
        ]
        return any(re.search(p, self.content) for p in patterns)

    def _check_other_tools(self) -> List[str]:
        """检查其他工具依赖"""
        tools = []
        if re.search(r'phonopy', self.content, re.IGNORECASE):
            tools.append('phonopy')
        if re.search(r'dpgen', self.content, re.IGNORECASE):
            tools.append('dpgen')
        return tools

    def _check_input_file(self) -> bool:
        """检查是否有完整的INPUT文件"""
        # 查找INPUT_PARAMETERS代码块
        pattern = r'```[^\n]*\n(INPUT_PARAMETERS[\s\S]+?)```'
        match = re.search(pattern, self.content)
        if match:
            input_content = match.group(1)
            # 检查是否包含关键参数
            required = ['calculation', 'ecutwfc', 'basis_type']
            return all(param in input_content for param in required)
        return False

    def _check_stru_file(self) -> bool:
        """检查是否有完整的STRU文件"""
        pattern = r'```[^\n]*\n(ATOMIC_SPECIES[\s\S]+?ATOMIC_POSITIONS[\s\S]+?)```'
        return bool(re.search(pattern, self.content))

    def _check_kpt_file(self) -> bool:
        """检查是否有完整的KPT文件"""
        pattern = r'```[^\n]*\n(K_POINTS[\s\S]+?)```'
        return bool(re.search(pattern, self.content))

    def _extract_input(self) -> Optional[str]:
        """提取INPUT文件内容"""
        pattern = r'```[^\n]*\n(INPUT_PARAMETERS[\s\S]+?)```'
        match = re.search(pattern, self.content)
        return match.group(1).strip() if match else None

    def _extract_stru(self) -> Optional[str]:
        """提取STRU文件内容"""
        pattern = r'```[^\n]*\n(ATOMIC_SPECIES[\s\S]+?)```'
        match = re.search(pattern, self.content)
        return match.group(1).strip() if match else None

    def _extract_kpt(self) -> Optional[str]:
        """提取KPT文件内容"""
        pattern = r'```[^\n]*\n(K_POINTS[\s\S]+?)```'
        match = re.search(pattern, self.content)
        return match.group(1).strip() if match else None

    def _extract_python_scripts(self) -> List[str]:
        """提取Python脚本"""
        pattern = r'```[Pp]ython\n([\s\S]+?)```'
        return re.findall(pattern, self.content)

    def _extract_abacustest_commands(self) -> List[str]:
        """提取abacustest命令"""
        pattern = r'abacustest\s+[^\n]+'
        return re.findall(pattern, self.content)

    def _extract_bash_commands(self) -> List[str]:
        """提取bash命令"""
        pattern = r'```[Bb]ash\n([\s\S]+?)```'
        matches = re.findall(pattern, self.content)
        commands = []
        for match in matches:
            commands.extend([line.strip() for line in match.split('\n') if line.strip()])
        return commands

    def _extract_expected_results(self) -> Dict:
        """提取预期结果"""
        results = {}

        # 提取弹性常数
        elastic_pattern = r'C[₁₂₃₄₅₆]{1,2}\s*=\s*([\d.]+)\s*GPa'
        elastic_matches = re.findall(elastic_pattern, self.content)
        if elastic_matches:
            results['elastic_constants'] = elastic_matches

        # 提取模量
        modulus_patterns = {
            'bulk_modulus': r'体模量.*?([\d.]+)\s*GPa',
            'shear_modulus': r'剪切模量.*?([\d.]+)\s*GPa',
            'young_modulus': r'杨氏模量.*?([\d.]+)\s*GPa',
            'poisson_ratio': r'泊松比.*?([\d.]+)',
        }
        for key, pattern in modulus_patterns.items():
            match = re.search(pattern, self.content)
            if match:
                results[key] = float(match.group(1))

        return results

    def print_summary(self, analysis: TutorialAnalysis):
        """打印分析摘要"""
        print("=" * 60)
        print(f"教程分析: {analysis.title}")
        print("=" * 60)
        print(f"\n主题: {analysis.topic}")

        print(f"\n工具依赖:")
        print(f"  - ASE: {'YES' if analysis.needs_ase else 'NO'}")
        print(f"  - abacustest: {'YES' if analysis.needs_abacustest else 'NO'}")
        if analysis.needs_other_tools:
            print(f"  - 其他: {', '.join(analysis.needs_other_tools)}")

        print(f"\n文件完整性:")
        print(f"  - INPUT: {'YES' if analysis.has_complete_input else 'NO'}")
        print(f"  - STRU: {'YES' if analysis.has_complete_stru else 'NO'}")
        print(f"  - KPT: {'YES' if analysis.has_complete_kpt else 'NO'}")

        print(f"\n提取的命令:")
        print(f"  - Python脚本: {len(analysis.python_scripts)} 个")
        print(f"  - abacustest命令: {len(analysis.abacustest_commands)} 个")
        print(f"  - Bash命令: {len(analysis.bash_commands)} 个")

        print(f"\n预期结果:")
        for key, value in analysis.expected_results.items():
            print(f"  - {key}: {value}")

        print("\n" + "=" * 60)


# 使用示例
if __name__ == "__main__":
    analyzer = TutorialAnalyzer("_workspace/20260203_105918_弹性常数计算 copy/07_final.md")
    analysis = analyzer.analyze()
    analyzer.print_summary(analysis)
