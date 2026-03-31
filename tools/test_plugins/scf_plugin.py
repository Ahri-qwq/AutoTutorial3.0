"""
AutoTutorial 3.0 - 普通 SCF 测试插件

适用教程类型：
  - ABACUS 入门 / 输入输出体系
  - 纯自洽场计算（calculation scf），不含 DFT+U / TDDFT 等特殊功能
  - PW 或 LCAO 基组均可
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


class SCFPlugin(BaseTestPlugin):
    """普通 SCF 测试插件"""

    @property
    def plugin_name(self) -> str:
        return "SCF 自洽场计算测试"

    @property
    def calc_type(self) -> str:
        return "scf"

    def can_handle(self, tutorial_content: str) -> bool:
        """
        判断是否是普通 SCF 教程。
        条件：
          1. 教程 INPUT 代码块中含 `calculation  scf`（或 `calculation=scf`）
          2. 不含 dft_plus_u（避免与 DFTUPlugin 冲突）
          3. 不含 esolver_type = tddft（避免与 TDDFTPlugin 冲突）
          4. 不含 esolver_type = sdft（避免与 SDFTPlugin 冲突）
          5. 教程主题是输入输出体系/三文件入门/收敛性测试 等基础入门类
        """
        # 排除特殊 esolver_type
        if re.search(r'esolver_type\s*=?\s*(tddft|sdft)', tutorial_content, re.IGNORECASE):
            return False
        # 排除 DFT+U
        if re.search(r'dft_plus_u\s*=?\s*1', tutorial_content, re.IGNORECASE):
            return False
        # 排除 phonopy/NEB 等教程（其 INPUT 也有 scf）
        exclude_patterns = [
            r'phonopy', r'声子谱', r'FORCE_SETS',
            r'AbacusNEB', r'AbacusAutoNEB', r'ATST-Tools',
        ]
        if any(re.search(p, tutorial_content, re.IGNORECASE) for p in exclude_patterns):
            return False

        # 必须在 INPUT 代码块中明确写 calculation  scf
        # 用正向匹配：INPUT_PARAMETERS 代码块内含 calculation ... scf
        input_blocks = re.findall(
            r'```[^\n]*\nINPUT_PARAMETERS[\s\S]+?```',
            tutorial_content
        )
        for block in input_blocks:
            if re.search(r'calculation\s+scf', block, re.IGNORECASE):
                return True

        return False

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取 SCF 测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name=self._extract_case_name(tutorial_content)
        )

        test_info.input_content = self._extract_input(tutorial_content)
        test_info.stru_content = self._extract_stru(tutorial_content)
        test_info.kpt_content = self._extract_kpt(tutorial_content)

        test_info.expected_results = self._extract_expected_results(tutorial_content)
        test_info.pseudopotentials = self._extract_pseudopotentials(tutorial_content)
        test_info.orbitals = self._extract_orbitals(tutorial_content)

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备输入文件"""
        input_dir = work_dir / f"{test_info.case_name}_scf"
        input_dir.mkdir(exist_ok=True)

        has_complete_files = (
            test_info.input_content and
            test_info.stru_content and
            test_info.kpt_content
        )

        if has_complete_files:
            print(f"  [策略] 从教程提取文件")
            (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
            (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
            (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        else:
            print(f"  [ERROR] 教程未提供完整 INPUT/STRU/KPT，无法准备输入文件")
            return []

        # 下载赝势和轨道文件
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                shutil.copy(pp_file, input_dir / pp)
            for orb in test_info.orbitals:
                orb_file = self.pp_manager.get_file(orb, "orbital")
                shutil.copy(orb_file, input_dir / orb)

        # 修复 INPUT：确保 pseudo_dir / orbital_dir 存在
        input_file = input_dir / "INPUT"
        if input_file.exists():
            print(f"    [修复] INPUT文件 - 检查目录路径")
            self._fix_input_file(input_file)

        # 修复路径格式
        print(f"    [修复] 路径格式 - 改为相对路径")
        from fix_stru import fix_stru_paths, fix_input_paths
        stru_file = input_dir / "STRU"
        if stru_file.exists():
            fix_stru_paths(stru_file, verbose=False)
        if input_file.exists():
            fix_input_paths(input_file, verbose=False)

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 SCF 任务"""
        job_ids = []
        for input_dir in input_dirs:
            job_name = f"{input_dir.name}"
            job_config = self.job_manager.create_job_config(job_name, input_dir)
            job_id = self.job_manager.submit_job(job_config, input_dir)
            if job_id:
                job_ids.append(job_id)
                print(f"  [OK] 提交任务: {job_name} (Job ID: {job_id})")
            else:
                print(f"  [ERROR] 提交失败: {job_name}")
        return job_ids

    def validate_results(
        self, job_ids: List[str], work_dir: Path, test_info: TestInfo
    ) -> ValidationResult:
        """验证 SCF 结果"""
        validation = ValidationResult(
            calc_type=self.calc_type,
            case_name=test_info.case_name,
            passed=False,
            comparisons={}
        )

        results_dir = work_dir / "results"
        results_dir.mkdir(exist_ok=True)

        for job_id in job_ids:
            self.job_manager.download_job_result(job_id, results_dir)

        actual_results = self._extract_actual_results(results_dir)

        if 'final_energy' in test_info.expected_results and 'final_energy' in actual_results:
            comparison = self._compare_values(
                test_info.expected_results['final_energy'],
                actual_results['final_energy'],
                'final_energy'
            )
            validation.comparisons['final_energy'] = comparison
        else:
            validation.warnings.append(
                "未找到可对比的总能量（expected 或 actual 缺失）"
            )

        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.passed = False

        validation.error_analyses = self.analyze_deviation(validation.comparisons)

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"

        if validation.comparisons:
            report += "| 参数 | 预期值 (eV) | 实际值 (eV) | 相对误差 | 状态 |\n"
            report += "|------|-------------|-------------|----------|------|\n"
            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                report += (
                    f"| {key} | {comp['expected']:.6f} | {comp['actual']:.6f} "
                    f"| {comp['rel_error']*100:.4f}% | {status} |\n"
                )

        if validation.errors:
            report += "\n**错误：**\n"
            for error in validation.errors:
                report += f"- {error}\n"

        if validation.warnings:
            report += "\n**警告：**\n"
            for warning in validation.warnings:
                report += f"- {warning}\n"

        if validation.error_analyses:
            report += "\n**误差归因分析：**\n\n"
            for key, ea in validation.error_analyses.items():
                report += f"- **{key}**（{ea.category}）：{ea.physical_explanation}\n"
                report += f"  - 用户建议：{ea.user_guidance}\n"

        overall = "✅ 通过" if validation.passed else "❌ 失败"
        report += f"\n**总体结果：** {overall}\n\n"

        return report

    # ── 辅助方法 ──────────────────────────────────────

    def _extract_case_name(self, content: str) -> str:
        """提取案例名称（使用基类通用方法）"""
        return self.extract_case_name(content)

    def _extract_input(self, content: str) -> Optional[str]:
        """
        提取 INPUT 文件。
        优先提取含 calculation scf 且含 smearing_method 的完整 Al INPUT（3.6节）。
        """
        blocks = re.findall(
            r'```[^\n]*\n(INPUT_PARAMETERS[\s\S]+?)```',
            content
        )
        # 优先选含 smearing_method mp（Al 金属示例）的 block
        for block in blocks:
            if re.search(r'smearing_method\s+mp', block, re.IGNORECASE):
                return block.strip()
        # 其次选含 calculation scf 的 block
        for block in blocks:
            if re.search(r'calculation\s+scf', block, re.IGNORECASE):
                return block.strip()
        return None

    def _extract_stru(self, content: str) -> Optional[str]:
        """
        提取 STRU 文件。
        选取包含完整结构（ATOMIC_SPECIES + LATTICE_CONSTANT + LATTICE_VECTORS +
        ATOMIC_POSITIONS）且为 PW 写法（不含 NUMERICAL_ORBITAL）的块。
        """
        # 匹配整个代码块（含多行）
        blocks = re.findall(
            r'```[^\n]*\n(ATOMIC_SPECIES[\s\S]*?ATOMIC_POSITIONS[\s\S]+?)```',
            content
        )
        for block in blocks:
            # PW 写法：无 NUMERICAL_ORBITAL 块
            if 'NUMERICAL_ORBITAL' not in block and 'Al' in block:
                return block.strip()
        # fallback：任意含 Al 的完整块
        for block in blocks:
            if 'Al' in block:
                return block.strip()
        return blocks[0].strip() if blocks else None

    def _extract_kpt(self, content: str) -> Optional[str]:
        """
        提取 KPT 文件。
        优先提取 8×8×8 Gamma 网格（Al 金属推荐值）。
        """
        blocks = re.findall(
            r'```[^\n]*\n(K_POINTS[\s\S]+?)```',
            content
        )
        for block in blocks:
            if '8  8  8' in block or '8 8 8' in block:
                return block.strip()
        # fallback：第一个 MP 网格（非 Line 格式）
        for block in blocks:
            if 'Line' not in block:
                return block.strip()
        return blocks[0].strip() if blocks else None

    def _extract_expected_results(self, content: str) -> Dict:
        """
        提取预期总能量。
        教程附录中有：!FINAL_ETOT_IS  -215.847068 eV
        """
        results = {}
        patterns = [
            r'FINAL_ETOT_IS\s+([+-]?\d+\.\d+)\s*eV',
            r'总能量.*?([+-]?\d{2,3}\.\d+)\s*eV',
            r'!FINAL_ETOT_IS\s+([+-]?\d+\.\d+)',
        ]
        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                results['final_energy'] = float(match.group(1))
                break
        return results

    def _extract_pseudopotentials(self, content: str) -> List[str]:
        """提取赝势文件名"""
        pattern = r'(\w+_[A-Z]+_[A-Z0-9]+-[\d.]+\.upf)'
        return list(set(re.findall(pattern, content)))

    def _extract_orbitals(self, content: str) -> List[str]:
        """提取轨道文件名（PW 计算通常不需要，但提取备用）"""
        pattern = r'(\w+_\w+_\d+au_\d+Ry_[\w\d]+\.orb)'
        return list(set(re.findall(pattern, content)))

    def _fix_input_file(self, input_path: Path):
        """确保 INPUT 含 pseudo_dir"""
        content = input_path.read_text(encoding='utf-8')
        lines = content.split('\n')
        has_pseudo_dir = any('pseudo_dir' in line.lower() for line in lines)
        if has_pseudo_dir:
            return
        new_lines = []
        added = False
        for line in lines:
            new_lines.append(line)
            if line.strip() == 'INPUT_PARAMETERS' and not added:
                new_lines.append('pseudo_dir      ./')
                added = True
        input_path.write_text('\n'.join(new_lines), encoding='utf-8')

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从 running_scf.log 中提取实际总能量"""
        results = {}
        log_files = list(results_dir.rglob("running_scf.log"))
        for log_file in log_files:
            content = log_file.read_text(encoding='utf-8', errors='ignore')
            # ABACUS 输出格式：!FINAL_ETOT_IS  -215.847068 eV
            match = re.search(
                r'!?FINAL_ETOT_IS\s+([+-]?\d+\.\d+)\s*eV',
                content,
                re.IGNORECASE
            )
            if match:
                results['final_energy'] = float(match.group(1))
                break
        return results
