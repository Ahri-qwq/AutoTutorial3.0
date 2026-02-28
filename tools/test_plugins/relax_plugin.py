"""
AutoTutorial 3.0 - 结构优化测试插件
"""

import re
import subprocess
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


class RelaxPlugin(BaseTestPlugin):
    """结构优化测试插件"""

    @property
    def plugin_name(self) -> str:
        return "结构优化测试"

    @property
    def calc_type(self) -> str:
        return "relax"

    def can_handle(self, tutorial_content: str) -> bool:
        """判断是否包含结构优化计算"""
        patterns = [
            r'calculation\s*=\s*(relax|cell-relax)',
            r'结构优化',
            r'晶格优化',
            r'原子弛豫'
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取结构优化测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name=self._extract_case_name(tutorial_content)
        )

        # 提取输入文件
        test_info.input_content = self._extract_input(tutorial_content)
        test_info.stru_content = self._extract_stru(tutorial_content)
        test_info.kpt_content = self._extract_kpt(tutorial_content)

        # 提取Python脚本
        test_info.python_scripts = self._extract_python_scripts(tutorial_content)

        # 提取预期结果
        test_info.expected_results = self._extract_expected_results(tutorial_content)

        # 提取赝势和轨道文件
        test_info.pseudopotentials = self._extract_pseudopotentials(tutorial_content)
        test_info.orbitals = self._extract_orbitals(tutorial_content)

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备输入文件"""
        input_dir = work_dir / f"{test_info.case_name}_relax"
        input_dir.mkdir(exist_ok=True)

        # 策略1: 尝试从教程提取完整文件
        has_complete_files = (test_info.input_content and
                             test_info.stru_content and
                             test_info.kpt_content)

        if has_complete_files:
            print(f"  [策略] 从教程提取文件")
            # 写入输入文件
            (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
            (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
            (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        else:
            # 策略2: 执行Python脚本 + abacustest
            print(f"  [策略] 执行脚本 + abacustest")

            # 执行Python脚本生成结构
            if test_info.python_scripts:
                for i, script in enumerate(test_info.python_scripts):
                    script_file = input_dir / f"generate_structure_{i}.py"
                    script_file.write_text(script, encoding='utf-8')

                    print(f"    [执行] Python脚本 {i}")
                    result = subprocess.run(
                        ['python', str(script_file)],
                        cwd=str(input_dir),
                        capture_output=True,
                        text=True
                    )

                    if result.returncode == 0:
                        print(f"    [OK] 结构生成成功")
                    else:
                        print(f"    [ERROR] 结构生成失败: {result.stderr}")
                        return []

            # 查找生成的CIF文件
            cif_files = list(input_dir.glob("*.cif"))
            if cif_files:
                cif_file = cif_files[0]
                print(f"    [找到] 结构文件: {cif_file.name}")

                # 调用abacustest生成输入文件
                print(f"    [调用] abacustest model inputs")
                result = subprocess.run(
                    ['abacustest', 'model', 'inputs',
                     '-f', str(cif_file),
                     '--ftype', 'cif',
                     '--jtype', 'cell-relax',
                     '--lcao',
                     '--folder-syntax', test_info.case_name],
                    cwd=str(input_dir),
                    capture_output=True,
                    text=True
                )

                if result.returncode == 0:
                    print(f"    [OK] abacustest 执行成功")

                    # 将生成的文件移到input_dir
                    gen_dir = input_dir / test_info.case_name
                    if gen_dir.exists():
                        import shutil
                        for file in gen_dir.glob("*"):
                            shutil.move(str(file), str(input_dir / file.name))
                        gen_dir.rmdir()
                else:
                    print(f"    [ERROR] abacustest 失败: {result.stderr}")
                    return []
            else:
                print(f"    [ERROR] 未找到CIF文件")
                return []

        # 下载赝势和轨道文件
        if self.pp_manager:
            import shutil
            for pp in test_info.pseudopotentials:
                pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                shutil.copy(pp_file, input_dir / pp)
            for orb in test_info.orbitals:
                orb_file = self.pp_manager.get_file(orb, "orbital")
                shutil.copy(orb_file, input_dir / orb)

        # 修复STRU文件 - 添加赝势和轨道文件路径
        stru_file = input_dir / "STRU"
        if stru_file.exists():
            print(f"    [修复] STRU文件 - 添加文件路径")
            self._fix_stru_file(stru_file, test_info.pseudopotentials, test_info.orbitals)

        # 修复INPUT文件 - 添加pseudo_dir和orbital_dir
        input_file = input_dir / "INPUT"
        if input_file.exists():
            print(f"    [修复] INPUT文件 - 添加目录路径")
            self._fix_input_file(input_file)

        # 修复路径格式 - 将绝对路径改为相对路径（解决Bohrium环境变量问题）
        print(f"    [修复] 路径格式 - 改为相对路径")
        from fix_stru import fix_stru_paths, fix_input_paths
        if stru_file.exists():
            fix_stru_paths(stru_file, verbose=False)
        if input_file.exists():
            fix_input_paths(input_file, verbose=False)

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交结构优化任务"""
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

    def validate_results(self, job_ids: List[str], work_dir: Path, test_info: TestInfo) -> ValidationResult:
        """验证结构优化结果"""
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

        # 对比结果
        if 'final_energy' in test_info.expected_results and 'final_energy' in actual_results:
            comparison = self._compare_values(
                test_info.expected_results['final_energy'],
                actual_results['final_energy'],
                'final_energy'
            )
            validation.comparisons['final_energy'] = comparison

        # 判断是否通过
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的结果")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"

        if validation.comparisons:
            report += "| 参数 | 预期值 | 实际值 | 相对误差 | 状态 |\n"
            report += "|------|--------|--------|----------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                report += f"| {key} | {comp['expected']:.6f} | {comp['actual']:.6f} | {comp['rel_error']*100:.2f}% | {status} |\n"

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
        # 尝试多种格式
        patterns = [
            r'案例\d*\s*[-:：]\s*(\w+)',  # 案例1 - Si 或 案例：Si
            r'#+\s*案例[：:]\s*(\w+)',     # ## 案例：Si
            r'计算(\w+)的',                 # 计算Si的弹性模量
            r'##\s*(\w+)的弹性',           # ## Si的弹性
        ]

        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                return match.group(1)

        # 默认名称
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

    def _extract_python_scripts(self, content: str) -> List[str]:
        """提取Python脚本"""
        pattern = r'```[Pp]ython\n([\s\S]+?)```'
        return re.findall(pattern, content)

    def _extract_expected_results(self, content: str) -> Dict:
        """提取预期结果"""
        results = {}

        # 提取最终能量
        patterns = {
            'final_energy': r'最终能量.*?([+-]?\d+\.\d+)\s*eV',
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                results[key] = float(match.group(1))

        return results

    def _extract_pseudopotentials(self, content: str) -> List[str]:
        """提取赝势文件名"""
        pattern = r'(\w+_[A-Z]+_[A-Z0-9]+-[\d.]+\.upf)'
        return list(set(re.findall(pattern, content)))

    def _extract_orbitals(self, content: str) -> List[str]:
        """提取轨道文件名"""
        pattern = r'(\w+_\w+_\d+au_\d+Ry_[\w\d]+\.orb)'
        return list(set(re.findall(pattern, content)))

    def _fix_stru_file(self, stru_path: Path, pseudopotentials: List[str], orbitals: List[str]):
        """
        修复STRU文件：
        1. 在ATOMIC_SPECIES行添加赝势文件名
        2. 添加NUMERICAL_ORBITAL块（只包含体系中存在的元素的轨道文件）

        Args:
            stru_path: STRU文件路径
            pseudopotentials: 赝势文件列表
            orbitals: 轨道文件列表
        """
        content = stru_path.read_text(encoding='utf-8')
        lines = content.split('\n')
        new_lines = []
        in_atomic_species = False
        numerical_orbital_added = False
        elements_in_system = set()  # 记录体系中的元素

        for line in lines:
            # 检测ATOMIC_SPECIES块
            if line.strip().startswith('ATOMIC_SPECIES'):
                in_atomic_species = True
                new_lines.append(line)
                continue

            # 处理ATOMIC_SPECIES块中的元素行
            if in_atomic_species and line.strip() and not line.strip().startswith('LATTICE') and not line.strip().startswith('NUMERICAL'):
                parts = line.split()
                if len(parts) >= 2 and len(parts) < 3:  # 只有元素和质量
                    element = parts[0]
                    mass = parts[1]
                    elements_in_system.add(element)  # 记录元素

                    # 查找对应的赝势文件
                    pp_file = None
                    for pp in pseudopotentials:
                        if element in pp:
                            pp_file = pp
                            break

                    if pp_file:
                        new_line = f"{element} {mass} {pp_file}"
                        new_lines.append(new_line)
                        print(f"      添加赝势 {element}: {pp_file}")
                    else:
                        new_lines.append(line)
                else:
                    new_lines.append(line)
                continue

            # 检测块结束
            if in_atomic_species and (line.strip().startswith('LATTICE') or line.strip().startswith('NUMERICAL')):
                in_atomic_species = False

            # 在LATTICE_CONSTANT之前添加NUMERICAL_ORBITAL块
            if not numerical_orbital_added and line.strip().startswith('LATTICE_CONSTANT'):
                if 'NUMERICAL_ORBITAL' not in content:
                    new_lines.append('NUMERICAL_ORBITAL')
                    # 只添加体系中存在的元素的轨道文件
                    added_orbitals = []
                    for orb in orbitals:
                        for element in elements_in_system:
                            if element in orb:
                                new_lines.append(f"./{orb}")  # 添加 ./ 前缀
                                added_orbitals.append(orb)
                                break
                    new_lines.append('')  # 空行
                    numerical_orbital_added = True
                    print(f"      添加NUMERICAL_ORBITAL块，包含 {len(added_orbitals)} 个轨道文件: {', '.join(added_orbitals)}")

            new_lines.append(line)

        # 写回文件
        stru_path.write_text('\n'.join(new_lines), encoding='utf-8')

    def _fix_input_file(self, input_path: Path):
        """
        修复INPUT文件，添加pseudo_dir和orbital_dir

        Args:
            input_path: INPUT文件路径
        """
        content = input_path.read_text(encoding='utf-8')
        lines = content.split('\n')

        # 检查是否已经有这些参数
        has_pseudo_dir = any('pseudo_dir' in line.lower() for line in lines)
        has_orbital_dir = any('orbital_dir' in line.lower() for line in lines)

        if has_pseudo_dir and has_orbital_dir:
            return  # 已经有了，不需要修复

        # 在INPUT_PARAMETERS后添加
        new_lines = []
        added = False

        for line in lines:
            new_lines.append(line)
            if line.strip() == 'INPUT_PARAMETERS' and not added:
                if not has_pseudo_dir:
                    new_lines.append('pseudo_dir      ./')
                    print(f"      添加 pseudo_dir")
                if not has_orbital_dir:
                    new_lines.append('orbital_dir     ./')
                    print(f"      添加 orbital_dir")
                added = True

        # 写回文件
        input_path.write_text('\n'.join(new_lines), encoding='utf-8')

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从计算结果中提取实际值"""
        results = {}

        # 查找 running_scf.log 或 running_relax.log
        log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            content = log_file.read_text(encoding='utf-8', errors='ignore')

            # 提取最终能量
            match = re.search(r'final etot is\s+([+-]?\d+\.\d+)\s+eV', content)
            if match:
                results['final_energy'] = float(match.group(1))
                break

        return results
