"""
AutoTutorial 3.0 - 态密度测试插件
"""

import re
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


class DOSPlugin(BaseTestPlugin):
    """态密度测试插件"""

    @property
    def plugin_name(self) -> str:
        return "态密度测试"

    @property
    def calc_type(self) -> str:
        return "dos"

    def can_handle(self, tutorial_content: str) -> bool:
        """判断是否包含态密度计算"""
        patterns = [
            r'态密度',
            r'DOS',
            r'density\s*of\s*states',
            r'out_dos\s*=\s*1',
            r'费米能级',
            r'Fermi\s*level'
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取态密度测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name=self._extract_case_name(tutorial_content)
        )

        # 提取输入文件
        test_info.input_content = self._extract_input(tutorial_content)
        test_info.stru_content = self._extract_stru(tutorial_content)
        test_info.kpt_content = self._extract_kpt(tutorial_content)

        # 提取预期结果
        test_info.expected_results = self._extract_expected_results(tutorial_content)

        # 提取赝势和轨道文件
        test_info.pseudopotentials = self._extract_pseudopotentials(tutorial_content)
        test_info.orbitals = self._extract_orbitals(tutorial_content)

        return test_info if test_info.input_content else None

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备态密度计算输入"""
        input_dir = work_dir / f"{test_info.case_name}_dos"
        input_dir.mkdir(exist_ok=True)

        # 写入输入文件
        if test_info.input_content:
            (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        if test_info.stru_content:
            (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        if test_info.kpt_content:
            (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')

        # 下载赝势和轨道文件
        if self.pp_manager:
            import shutil
            for pp in test_info.pseudopotentials:
                pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                shutil.copy(pp_file, input_dir / pp)
            for orb in test_info.orbitals:
                orb_file = self.pp_manager.get_file(orb, "orbital")
                shutil.copy(orb_file, input_dir / orb)

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交态密度计算任务"""
        job_ids = []

        for input_dir in input_dirs:
            job_name = input_dir.name
            job_config = self.job_manager.create_job_config(job_name, input_dir)
            job_id = self.job_manager.submit_job(job_config, input_dir)

            if job_id:
                job_ids.append(job_id)
                print(f"  [OK] 提交任务: {job_name} (Job ID: {job_id})")
            else:
                print(f"  [ERROR] 提交失败: {job_name}")

        return job_ids

    def validate_results(self, job_ids: List[str], work_dir: Path, test_info: TestInfo) -> ValidationResult:
        """验证态密度结果"""
        validation = ValidationResult(
            calc_type=self.calc_type,
            case_name=test_info.case_name,
            passed=False,
            comparisons={}
        )

        # 下载结果
        results_dir = work_dir / "results"
        results_dir.mkdir(exist_ok=True)

        for job_id in job_ids:
            self.job_manager.download_job_result(job_id, results_dir)

        # 提取实际结果
        actual_results = self._extract_actual_results(results_dir)

        # 对比费米能级
        if 'fermi_level' in test_info.expected_results and 'fermi_level' in actual_results:
            comparison = self._compare_values(
                test_info.expected_results['fermi_level'],
                actual_results['fermi_level'],
                'fermi_level'
            )
            validation.comparisons['fermi_level'] = comparison

        # 判断是否通过
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的费米能级")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"

        if validation.comparisons:
            report += "| 参数 | 预期值 (eV) | 实际值 (eV) | 相对误差 | 状态 |\n"
            report += "|------|-------------|-------------|----------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                report += f"| {key} | {comp['expected']:.4f} | {comp['actual']:.4f} | {comp['rel_error']*100:.2f}% | {status} |\n"

        if validation.errors:
            report += "\n**错误：**\n"
            for error in validation.errors:
                report += f"- {error}\n"

        if validation.warnings:
            report += "\n**警告：**\n"
            for warning in validation.warnings:
                report += f"- {warning}\n"

        overall = "✅ 通过" if validation.passed else "❌ 失败"
        report += f"\n**总体结果：** {overall}\n\n"

        return report

    # 辅助方法
    def _extract_case_name(self, content: str) -> str:
        """提取案例名称"""
        patterns = [
            r'案例\d*\s*[-:：]\s*(\w+)',
            r'#+\s*案例[：:]\s*(\w+)',
            r'计算(\w+)的',
        ]

        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                return match.group(1)

        return "Unknown"

    def _extract_input(self, content: str) -> Optional[str]:
        """提取 INPUT 文件"""
        pattern = r'```[^\n]*\n(INPUT_PARAMETERS[\s\S]+?)```'
        match = re.search(pattern, content)
        return match.group(1).strip() if match else None

    def _extract_stru(self, content: str) -> Optional[str]:
        """提取 STRU 文件"""
        pattern = r'```[^\n]*\n(ATOMIC_SPECIES[\s\S]+?)```'
        match = re.search(pattern, content)
        return match.group(1).strip() if match else None

    def _extract_kpt(self, content: str) -> Optional[str]:
        """提取 KPT 文件"""
        pattern = r'```[^\n]*\n(K_POINTS[\s\S]+?)```'
        match = re.search(pattern, content)
        return match.group(1).strip() if match else None

    def _extract_expected_results(self, content: str) -> Dict:
        """提取预期结果"""
        results = {}

        # 提取费米能级
        fermi_patterns = [
            r'费米能级.*?([+-]?\d+\.?\d*)\s*eV',
            r'Fermi\s*level.*?([+-]?\d+\.?\d*)\s*eV',
            r'E_Fermi.*?([+-]?\d+\.?\d*)'
        ]
        for pattern in fermi_patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                results['fermi_level'] = float(match.group(1))
                break

        return results

    def _extract_pseudopotentials(self, content: str) -> List[str]:
        """提取赝势文件列表"""
        pp_pattern = r'(\w+_[\w-]+\.upf)'
        return list(set(re.findall(pp_pattern, content)))

    def _extract_orbitals(self, content: str) -> List[str]:
        """提取轨道文件列表"""
        orb_pattern = r'(\w+_gga_[\w.]+\.orb)'
        return list(set(re.findall(orb_pattern, content)))

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从计算结果中提取实际值"""
        results = {}

        # 查找 running_scf.log 文件
        log_files = list(results_dir.rglob("running_scf.log"))

        if log_files:
            log_content = log_files[0].read_text(encoding='utf-8', errors='ignore')

            # 提取费米能级
            fermi_match = re.search(r'EFERMI\s*=\s*([+-]?\d+\.?\d*)', log_content)
            if fermi_match:
                results['fermi_level'] = float(fermi_match.group(1))

        return results
