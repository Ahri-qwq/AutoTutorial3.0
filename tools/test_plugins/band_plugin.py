"""
AutoTutorial 3.0 - 能带结构测试插件
"""

import re
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


class BandPlugin(BaseTestPlugin):
    """能带结构测试插件"""

    @property
    def plugin_name(self) -> str:
        return "能带结构测试"

    @property
    def calc_type(self) -> str:
        return "band"

    def can_handle(self, tutorial_content: str) -> bool:
        """判断是否包含能带计算"""
        patterns = [
            r'能带',
            r'band\s*structure',
            r'calculation\s*=\s*nscf',
            r'带隙',
            r'band\s*gap',
            r'out_band\s*=\s*1'
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取能带测试信息"""
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
        """准备能带计算输入"""
        # 能带计算通常需要两步：SCF + NSCF
        scf_dir = work_dir / f"{test_info.case_name}_scf"
        nscf_dir = work_dir / f"{test_info.case_name}_nscf"

        scf_dir.mkdir(exist_ok=True)
        nscf_dir.mkdir(exist_ok=True)

        # 准备 SCF 输入
        if test_info.input_content:
            # 修改为 SCF 计算
            scf_input = test_info.input_content.replace('calculation nscf', 'calculation scf')
            scf_input = scf_input.replace('calculation = nscf', 'calculation = scf')
            (scf_dir / "INPUT").write_text(scf_input, encoding='utf-8')

        if test_info.stru_content:
            (scf_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')

        # SCF 使用常规 K 点
        if test_info.kpt_content:
            # 提取或生成 SCF K 点
            scf_kpt = self._generate_scf_kpt(test_info.kpt_content)
            (scf_dir / "KPT").write_text(scf_kpt, encoding='utf-8')

        # 准备 NSCF 输入（能带计算）
        if test_info.input_content:
            (nscf_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        if test_info.stru_content:
            (nscf_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        if test_info.kpt_content:
            (nscf_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')

        # 下载赝势和轨道文件
        if self.pp_manager:
            import shutil
            for pp in test_info.pseudopotentials:
                pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                shutil.copy(pp_file, scf_dir / pp)
                shutil.copy(pp_file, nscf_dir / pp)
            for orb in test_info.orbitals:
                orb_file = self.pp_manager.get_file(orb, "orbital")
                shutil.copy(orb_file, scf_dir / orb)
                shutil.copy(orb_file, nscf_dir / orb)

        return [scf_dir, nscf_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交能带计算任务"""
        job_ids = []

        # 先提交 SCF，再提交 NSCF
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
        """验证能带结果"""
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

        # 对比带隙
        if 'band_gap' in test_info.expected_results and 'band_gap' in actual_results:
            comparison = self._compare_values(
                test_info.expected_results['band_gap'],
                actual_results['band_gap'],
                'band_gap'
            )
            validation.comparisons['band_gap'] = comparison

        # 判断是否通过
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的带隙值")

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

        # 提取带隙
        patterns = [
            r'带隙[：:]\s*(\d+\.?\d*)\s*eV',
            r'band\s*gap[：:]\s*(\d+\.?\d*)\s*eV',
            r'Eg\s*=\s*(\d+\.?\d*)\s*eV'
        ]

        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                results['band_gap'] = float(match.group(1))
                break

        return results

    def _extract_pseudopotentials(self, content: str) -> List[str]:
        """提取赝势文件列表"""
        pattern = r'(\w+_[A-Z]+_[A-Z0-9]+-[\d.]+\.upf)'
        return list(set(re.findall(pattern, content)))

    def _extract_orbitals(self, content: str) -> List[str]:
        """提取轨道文件列表"""
        pattern = r'(\w+_gga_\w+\.orb)'
        return list(set(re.findall(pattern, content)))

    def _generate_scf_kpt(self, nscf_kpt: str) -> str:
        """从 NSCF KPT 生成 SCF KPT"""
        # 简单处理：如果是 Line 模式，改为 Gamma 或 MP
        if 'Line' in nscf_kpt or 'line' in nscf_kpt:
            return """K_POINTS
0
Gamma
8 8 8 0 0 0
"""
        return nscf_kpt

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从计算结果中提取带隙"""
        results = {}

        # 查找 running_scf.log 或 running_nscf.log
        log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            content = log_file.read_text(encoding='utf-8', errors='ignore')

            # 提取带隙
            match = re.search(r'E_bandgap\s*=\s*(\d+\.?\d*)\s*eV', content)
            if match:
                results['band_gap'] = float(match.group(1))
                break

        return results
