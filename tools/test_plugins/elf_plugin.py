"""
AutoTutorial 3.0 - ELF 电子局域函数测试插件
处理 ABACUS 输出 ELF.cube 文件的计算（三案例：H2O-PW、H2O-LCAO、Fe-BCC）
首次测试教程：ABACUS 使用教程｜电子局域函数（ELF）计算与可视化
"""

import re
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


class ELFPlugin(BaseTestPlugin):
    """ELF 电子局域函数测试插件"""

    @property
    def plugin_name(self) -> str:
        return "ELF 电子局域函数测试"

    @property
    def calc_type(self) -> str:
        return "elf"

    def can_handle(self, tutorial_content: str) -> bool:
        """判断是否包含 ELF 计算"""
        patterns = [
            r'out_elf\s+1',
            r'ELF\.cube',
            r'电子局域函数',
            r'Electron\s+Localization\s+Function'
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """
        提取 ELF 测试信息
        将所有案例数据存入 expected_results['cases']，供 prepare_inputs 使用
        """
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="ELF_multi_cases"
        )

        # 提取所有案例数据
        cases = self._extract_all_cases(tutorial_content)

        if not cases:
            return None

        # 将案例数据存入 expected_results 供 prepare_inputs 使用
        test_info.expected_results['cases'] = cases
        test_info.expected_results['case_count'] = len(cases)

        # 用第一个案例填充标准字段（满足框架需求）
        first = cases[0]
        test_info.input_content = first['input']
        test_info.stru_content = first['stru']
        test_info.kpt_content = first.get('kpt', '')
        test_info.pseudopotentials = first['pseudopotentials']
        test_info.orbitals = first.get('orbitals', [])

        print(f"  [ELF] 提取到 {len(cases)} 个案例：{[c['name'] for c in cases]}")

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """
        准备所有 ELF 案例的输入文件
        从 test_info.expected_results['cases'] 中读取案例数据
        """
        input_dirs = []
        cases = test_info.expected_results.get('cases', [])

        for case in cases:
            input_dir = work_dir / case['name']
            input_dir.mkdir(exist_ok=True)

            # 写入输入文件
            (input_dir / "INPUT").write_text(case['input'], encoding='utf-8')
            (input_dir / "STRU").write_text(case['stru'], encoding='utf-8')
            if case.get('kpt'):
                (input_dir / "KPT").write_text(case['kpt'], encoding='utf-8')

            # 下载赝势文件
            if self.pp_manager:
                import shutil
                for pp in case['pseudopotentials']:
                    try:
                        pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                        shutil.copy(pp_file, input_dir / pp)
                    except Exception as e:
                        print(f"  [WARN] 赝势文件 {pp} 下载失败: {e}")

                # 下载轨道文件（LCAO 案例）
                for orb in case.get('orbitals', []):
                    try:
                        orb_file = self.pp_manager.get_file(orb, "orbital")
                        shutil.copy(orb_file, input_dir / orb)
                    except Exception as e:
                        print(f"  [WARN] 轨道文件 {orb} 下载失败: {e}")

            print(f"  [ELF] 已准备案例: {case['name']} → {input_dir}")
            input_dirs.append(input_dir)

        return input_dirs

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 ELF 计算任务（所有案例可并行提交）"""
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
        """
        验证 ELF 结果
        - 检查 ELF.cube 文件是否存在且非空
        - 对于 nspin=2 的案例（Fe），检查 ELF_SPIN1.cube / ELF_SPIN2.cube
        """
        validation = ValidationResult(
            calc_type=self.calc_type,
            case_name="ELF_multi_cases",
            passed=False,
            comparisons={}
        )

        # 下载结果
        results_dir = work_dir / "results"
        results_dir.mkdir(exist_ok=True)

        for job_id in job_ids:
            self.job_manager.download_job_result(job_id, results_dir)

        cases = test_info.expected_results.get('cases', [])
        all_passed = True

        for case in cases:
            case_name = case['name']
            case_results_dir = results_dir / case_name
            is_spin = case.get('nspin', 1) == 2

            # 查找 ELF.cube（在 OUT.*/ELF.cube 路径下）
            elf_cubes = list(case_results_dir.rglob("ELF.cube")) if case_results_dir.exists() else []
            # 也查找根级别
            if not elf_cubes:
                elf_cubes = list(results_dir.rglob(f"*{case_name}*/ELF.cube"))

            cube_ok = len(elf_cubes) > 0 and elf_cubes[0].stat().st_size > 100

            if is_spin:
                spin1_cubes = list(results_dir.rglob(f"*{case_name}*/ELF_SPIN1.cube"))
                spin2_cubes = list(results_dir.rglob(f"*{case_name}*/ELF_SPIN2.cube"))
                spin_ok = bool(spin1_cubes) and bool(spin2_cubes)
                case_passed = cube_ok and spin_ok
                validation.comparisons[case_name] = {
                    'ELF.cube': '✅' if cube_ok else '❌',
                    'ELF_SPIN1.cube': '✅' if spin1_cubes else '❌',
                    'ELF_SPIN2.cube': '✅' if spin2_cubes else '❌',
                    'passed': case_passed
                }
            else:
                case_passed = cube_ok
                validation.comparisons[case_name] = {
                    'ELF.cube': '✅' if cube_ok else '❌',
                    'passed': case_passed
                }

            if not case_passed:
                all_passed = False
                validation.errors.append(f"{case_name}: ELF cube 文件缺失或为空")

        validation.passed = all_passed
        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成 ELF 测试报告"""
        status = "✅ 通过" if validation.passed else "❌ 失败"

        report = f"### {self.plugin_name}\n\n"
        report += f"**状态：** {status}\n\n"
        report += "**ELF 输出文件检查：**\n\n"
        report += "| 案例 | ELF.cube | ELF_SPIN1.cube | ELF_SPIN2.cube | 状态 |\n"
        report += "|------|----------|----------------|----------------|------|\n"

        for case_name, comp in validation.comparisons.items():
            spin1 = comp.get('ELF_SPIN1.cube', 'N/A')
            spin2 = comp.get('ELF_SPIN2.cube', 'N/A')
            overall = "✅ PASS" if comp.get('passed') else "❌ FAIL"
            report += f"| {case_name} | {comp['ELF.cube']} | {spin1} | {spin2} | {overall} |\n"

        if validation.errors:
            report += "\n**错误：**\n"
            for err in validation.errors:
                report += f"- {err}\n"

        return report

    # ========== 辅助方法 ==========

    def _extract_all_cases(self, content: str) -> List[Dict]:
        """从教程中提取所有 ELF 计算案例"""
        cases = []

        # ── 案例一：H2O PW（第3.1节） ──
        # 截取 3.1 节区域（从 3.1 到 3.2）
        sec_31 = self._extract_section_text(content, r'### 3\.1', r'### 3\.2')
        if sec_31:
            input_text = self._extract_code_block(sec_31, 'INPUT_PARAMETERS')
            stru_text = self._extract_code_block(sec_31, 'ATOMIC_SPECIES')
            kpt_text = self._extract_code_block(sec_31, 'K_POINTS')

            if input_text and stru_text:
                cases.append({
                    'name': 'H2O_PW',
                    'input': input_text,
                    'stru': stru_text,
                    'kpt': kpt_text,
                    'nspin': 1,
                    'pseudopotentials': ['H_ONCV_PBE-1.0.upf', 'O_ONCV_PBE-1.0.upf'],
                    'orbitals': []
                })

        # ── 案例二：H2O LCAO（第3.2节） ──
        sec_32 = self._extract_section_text(content, r'### 3\.2', r'---\s*\n## 四')
        if sec_32:
            input_text = self._extract_code_block(sec_32, 'INPUT_PARAMETERS')
            stru_text = self._extract_code_block(sec_32, 'ATOMIC_SPECIES')

            if input_text and stru_text:
                cases.append({
                    'name': 'H2O_LCAO',
                    'input': input_text,
                    'stru': stru_text,
                    'kpt': 'K_POINTS\n0\nGamma\n1 1 1 0 0 0',
                    'nspin': 1,
                    'pseudopotentials': ['H_ONCV_PBE-1.0.upf', 'O_ONCV_PBE-1.0.upf'],
                    'orbitals': ['H_gga_6au_100Ry_2s1p.orb', 'O_gga_7au_100Ry_2s2p1d.orb']
                })

        # ── 案例三：BCC Fe（第4章） ──
        sec_4 = self._extract_section_text(content, r'## 四、案例二', r'## 五、')
        if sec_4:
            input_text = self._extract_code_block(sec_4, 'INPUT_PARAMETERS')
            stru_text = self._extract_code_block(sec_4, 'ATOMIC_SPECIES')
            kpt_text = self._extract_code_block(sec_4, 'K_POINTS')

            if input_text and stru_text:
                cases.append({
                    'name': 'Fe_BCC',
                    'input': input_text,
                    'stru': stru_text,
                    'kpt': kpt_text,
                    'nspin': 2,
                    'pseudopotentials': ['Fe_ONCV_PBE-1.0.upf'],
                    'orbitals': []
                })

        return cases

    def _extract_section_text(self, content: str, start_pattern: str, end_pattern: str) -> Optional[str]:
        """提取两个标题之间的文本"""
        start_match = re.search(start_pattern, content)
        if not start_match:
            return None
        start = start_match.start()

        end_match = re.search(end_pattern, content[start:])
        if end_match:
            return content[start:start + end_match.start()]
        else:
            return content[start:]

    def _extract_code_block(self, text: str, start_keyword: str) -> Optional[str]:
        """从文本中提取以指定关键词开头的代码块"""
        # 匹配 ``` 开头（可以有语言标识），内容以 start_keyword 开头
        pattern = r'```[^\n]*\n(' + re.escape(start_keyword) + r'[\s\S]*?)```'
        match = re.search(pattern, text)
        return match.group(1).strip() if match else None
