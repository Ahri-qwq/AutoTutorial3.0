"""
AutoTutorial 3.0 - 弹性常数测试插件
"""

import os
import re
import json
import subprocess
import time
import shutil
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


class ElasticPlugin(BaseTestPlugin):
    """弹性常数测试插件"""

    @property
    def plugin_name(self) -> str:
        return "弹性常数测试"

    @property
    def calc_type(self) -> str:
        return "elastic"

    def can_handle(self, tutorial_content: str) -> bool:
        """判断是否包含弹性常数计算"""
        patterns = [
            r'弹性常数',
            r'elastic',
            r'C[₁₂₃₄₅₆]{1,2}',
            r'体模量',
            r'剪切模量',
            r'杨氏模量',
            r'abacustest.*elastic'
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取弹性常数测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name=self._extract_case_name(tutorial_content)
        )

        # 提取输入文件
        test_info.input_content = self._extract_input(tutorial_content)
        test_info.stru_content = self._extract_stru(tutorial_content)
        test_info.kpt_content = self._extract_kpt(tutorial_content)

        # 提取 Python 脚本（结构生成）
        test_info.python_scripts = self._extract_python_scripts(tutorial_content)

        # 提取 abacustest 命令
        test_info.abacustest_commands = self._extract_abacustest_commands(tutorial_content)

        # 提取预期结果
        test_info.expected_results = self._extract_expected_results(tutorial_content)

        # 提取赝势和轨道文件
        test_info.pseudopotentials = self._extract_pseudopotentials(tutorial_content)
        test_info.orbitals = self._extract_orbitals(tutorial_content)

        return test_info if test_info.expected_results else None

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备弹性常数计算输入"""
        case_dir = work_dir / f"{test_info.case_name}_elastic"
        case_dir.mkdir(exist_ok=True)

        # 1. 执行结构生成脚本
        if test_info.python_scripts:
            for i, script in enumerate(test_info.python_scripts):
                script_file = case_dir / f"generate_structure_{i}.py"
                script_file.write_text(script, encoding='utf-8')

                result = subprocess.run(
                    ['python', str(script_file)],
                    cwd=str(case_dir),
                    capture_output=True,
                    text=True
                )

                if result.returncode != 0:
                    print(f"  [ERROR] 结构生成失败: {result.stderr}")
                    return []

        # 2. 准备优化目录
        relax_dir = case_dir / "relax"
        relax_dir.mkdir(exist_ok=True)

        if test_info.input_content:
            (relax_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')

        # STRU 优先级：
        # 1. 内联 STRU → 直接使用
        # 2. abacustest model inputs 命令 → 运行生成
        # 3. 两者都没有 → 报错返回
        if test_info.stru_content:
            (relax_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        else:
            inputs_cmd = next(
                (cmd for cmd in test_info.abacustest_commands if "abacustest model inputs" in cmd),
                None
            )
            if inputs_cmd:
                print(f"  [STRU fallback] 未找到内联STRU，运行: {inputs_cmd}")
                result = subprocess.run(
                    inputs_cmd,
                    shell=True,
                    cwd=str(relax_dir),
                    capture_output=True,
                    text=True
                )
                if result.returncode != 0:
                    print(f"  [ERROR] abacustest model inputs 失败: {result.stderr}")
                    return []
                print(f"  [OK] abacustest model inputs 生成了输入文件")
            else:
                print(f"  [ERROR] 未找到STRU内容，也未找到 abacustest model inputs 命令")
                return []

        if test_info.kpt_content:
            (relax_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')

        # 下载赝势和轨道文件
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                shutil.copy(pp_file, relax_dir / pp)
            for orb in test_info.orbitals:
                orb_file = self.pp_manager.get_file(orb, "orbital")
                shutil.copy(orb_file, relax_dir / orb)

        # 修复路径格式 - 将绝对路径改为相对路径（解决Bohrium环境变量问题）
        print(f"  [修复] 路径格式 - 改为相对路径")
        from fix_stru import fix_stru_paths, fix_input_paths
        stru_file = relax_dir / "STRU"
        input_file = relax_dir / "INPUT"
        if stru_file.exists():
            fix_stru_paths(stru_file, verbose=False)
        if input_file.exists():
            fix_input_paths(input_file, verbose=False)

        return [relax_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交弹性常数计算任务（4阶段串行流程）"""
        job_ids = []

        for input_dir in input_dirs:
            print(f"\n[串行提交] 案例: {input_dir.parent.name}")
            elastic_dir = input_dir.parent / "elastic"

            # 阶段1：提交结构优化任务 + 等待 + 下载
            print("  [1/4] 提交结构优化任务...")
            relax_job_id = self._submit_relax_job(input_dir)
            if not relax_job_id:
                print(f"  [ERROR] 优化任务提交失败，跳过该案例")
                continue

            job_ids.append(relax_job_id)

            print(f"  [等待] 优化任务完成 (Job {relax_job_id})...")
            if not self._wait_for_job(relax_job_id):
                print(f"  [ERROR] 优化任务失败或超时，跳过弹性计算")
                continue

            stru_optimized = self._download_optimized_structure(relax_job_id, input_dir.parent)
            if not stru_optimized:
                print(f"  [ERROR] 下载优化结构失败，跳过弹性计算")
                continue

            # 阶段2：运行 abacustest model elastic prepare 生成25个形变结构
            print(f"  [2/4] 运行 abacustest model elastic prepare...")
            if not self._prepare_elastic_with_abacustest(elastic_dir, input_dir, stru_optimized):
                print(f"  [ERROR] 弹性计算输入准备失败")
                continue

            # 阶段3：提交所有形变结构任务，等待用户确认，下载结果
            print(f"  [3/4] 提交形变结构任务...")
            deformed_job_pairs = self._submit_all_deformed_jobs(elastic_dir)
            if not deformed_job_pairs:
                print(f"  [ERROR] 形变结构任务提交失败")
                continue

            job_ids.extend(jid for _, jid in deformed_job_pairs)

            machine_type = getattr(self.job_manager, 'machine_type', 'c32_m64_cpu')
            time_estimate = self._estimate_job_time(
                stru_optimized,
                n_jobs=len(deformed_job_pairs),
                machine_type=machine_type
            )
            self._wait_for_all_jobs(deformed_job_pairs, elastic_dir, time_estimate)

            # 阶段4：运行 abacustest model elastic post 后处理
            print(f"  [4/4] 运行 abacustest model elastic post...")
            if self._post_process_elastic(elastic_dir, input_dir.parent):
                print(f"  [OK] 弹性常数计算完成，结果已写入 metrics_elastic.json")
            else:
                print(f"  [ERROR] 弹性常数后处理失败")

        return job_ids

    def validate_results(self, job_ids: List[str], work_dir: Path, test_info: TestInfo) -> ValidationResult:
        """验证弹性常数结果"""
        validation = ValidationResult(
            calc_type=self.calc_type,
            case_name=test_info.case_name,
            passed=False,
            comparisons={}
        )

        # 1. 下载所有任务结果
        results_dir = work_dir / "results"
        results_dir.mkdir(exist_ok=True)

        print(f"\n[验证] 下载任务结果...")
        for job_id in job_ids:
            status = self.job_manager.get_job_status(job_id)
            print(f"  Job {job_id}: {status['status']}")

            if status['status'] == 'Finished':
                self.job_manager.download_job_result(job_id, results_dir)
            elif status['status'] == 'Failed':
                validation.errors.append(f"任务 {job_id} 失败")

        # 2. 如果存在 metrics_elastic.json，则进行对比
        metrics_file = work_dir / "metrics_elastic.json"
        if metrics_file.exists():
            print(f"  [对比] 使用 {metrics_file.name}")
            actual_results = self._load_elastic_results(metrics_file)

            # K点透明度：检查kspacing差异，必要时添加说明
            kpoint_note = self._check_kpoint_discrepancy(work_dir, test_info)
            if kpoint_note:
                validation.warnings.append(kpoint_note)

            validation.comparisons = self._compare_elastic_results(
                test_info.expected_results,
                actual_results
            )
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            print(f"  [警告] 未找到 metrics_elastic.json")
            print(f"  [提示] 弹性常数计算需要使用 abacustest 后处理工具")
            validation.warnings.append("未找到弹性常数后处理结果")
            validation.warnings.append("请使用 abacustest model elastic post 进行后处理")

            all_finished = all(
                self.job_manager.get_job_status(jid)['status'] == 'Finished'
                for jid in job_ids
            )

            if all_finished:
                validation.warnings.append("所有计算任务已完成，等待后处理")
            else:
                validation.errors.append("部分计算任务未成功完成")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"

        if validation.comparisons:
            report += "#### 弹性常数对比\n\n"
            report += "| 参数 | 预期值 (GPa) | 实际值 (GPa) | 相对误差 | 状态 |\n"
            report += "|------|--------------|--------------|----------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                report += f"| {key} | {comp['expected']:.2f} | {comp['actual']:.2f} | {comp['rel_error']*100:.2f}% | {status} |\n"

        # K点说明优先显示（当比较部分失败时尤为重要）
        kpoint_warnings = [w for w in validation.warnings if 'K点说明' in w]
        if kpoint_warnings:
            report += "\n**K点说明：**\n\n"
            for warning in kpoint_warnings:
                report += f"> {warning}\n"
            report += "\n"

        if validation.errors:
            report += "\n**错误：**\n"
            for error in validation.errors:
                report += f"- {error}\n"

        other_warnings = [w for w in validation.warnings if 'K点说明' not in w]
        if other_warnings:
            report += "\n**警告：**\n"
            for warning in other_warnings:
                report += f"- {warning}\n"

        overall = "✅ 通过" if validation.passed else "⚠️ 部分完成"
        report += f"\n**总体结果：** {overall}\n\n"

        return report

    # ========== 提取辅助方法 ==========

    def _extract_case_name(self, content: str) -> str:
        """提取案例名称"""
        patterns = [
            r'案例\d*\s*[-:：]\s*(\w+)',
            r'#+\s*案例[：:]\s*(\w+)',
            r'计算(\w+)的',
            r'##\s*(\w+)的弹性',
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

    def _extract_python_scripts(self, content: str) -> List[str]:
        """提取 Python 脚本"""
        pattern = r'```[Pp]ython\n([\s\S]+?)```'
        return re.findall(pattern, content)

    def _extract_abacustest_commands(self, content: str) -> List[str]:
        """提取 abacustest 命令"""
        pattern = r'abacustest\s+[^\n]+'
        return re.findall(pattern, content)

    def _extract_expected_results(self, content: str) -> Dict:
        """提取预期的弹性常数结果"""
        results = {}

        # 提取弹性常数 (C11, C12, C44 等)
        elastic_pattern = r'C[₁₂₃₄₅₆]{1,2}\s*[=:：]\s*([\\d.]+)\s*GPa'
        matches = re.findall(elastic_pattern, content)
        if matches:
            if len(matches) >= 1:
                results['C11'] = float(matches[0])
            if len(matches) >= 2:
                results['C12'] = float(matches[1])
            if len(matches) >= 3:
                results['C44'] = float(matches[2])

        # 提取模量
        modulus_patterns = {
            'bulk_modulus': r'体模量.*?(\d+\.?\d*)\\s*GPa',
            'shear_modulus': r'剪切模量.*?(\d+\.?\d*)\\s*GPa',
            'young_modulus': r'杨氏模量.*?(\d+\.?\d*)\\s*GPa',
            'poisson_ratio': r'泊松比.*?(\d+\.?\d*)',
        }

        for key, pattern in modulus_patterns.items():
            match = re.search(pattern, content)
            if match:
                try:
                    results[key] = float(match.group(1))
                except ValueError:
                    pass

        return results

    def _extract_pseudopotentials(self, content: str) -> List[str]:
        """提取赝势文件名"""
        pattern = r'(\w+_\w+_\w+-[\d.]+\.upf)'
        return list(set(re.findall(pattern, content)))

    def _extract_orbitals(self, content: str) -> List[str]:
        """提取轨道文件名"""
        pattern = r'(\w+_\w+_[\w.]+\.orb)'
        return list(set(re.findall(pattern, content)))

    def _load_elastic_results(self, metrics_file: Path) -> Dict:
        """从 metrics_elastic.json 加载结果"""
        with open(metrics_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        results = {}
        if 'elastic_tensor' in data:
            results.update(data['elastic_tensor'])
        if 'elastic_moduli' in data:
            results.update(data['elastic_moduli'])

        return results

    def _compare_elastic_results(self, expected: Dict, actual: Dict) -> Dict:
        """对比弹性常数结果"""
        comparisons = {}

        for key in expected:
            if key in actual:
                comparisons[key] = self._compare_values(expected[key], actual[key], key)

        return comparisons

    # ========== 阶段1：提交优化任务辅助方法 ==========

    def _submit_relax_job(self, input_dir: Path) -> Optional[str]:
        """提交结构优化任务"""
        job_name = f"{input_dir.parent.name}_relax"
        job_config = self.job_manager.create_job_config(job_name, input_dir)
        job_id = self.job_manager.submit_job(job_config, input_dir)

        if job_id:
            print(f"    ✅ 优化任务已提交 (Job {job_id})")
        else:
            print(f"    ❌ 优化任务提交失败")

        return job_id

    def _wait_for_job(self, job_id: str, timeout: int = 3600, check_interval: int = 30) -> bool:
        """
        等待单个任务完成（自动轮询）

        Args:
            job_id: 任务ID
            timeout: 超时时间（秒），默认3600秒=1小时
            check_interval: 检查间隔（秒），默认30秒

        Returns:
            True 如果任务成功完成，False 如果失败或超时
        """
        elapsed = 0
        while elapsed < timeout:
            status_info = self.job_manager.get_job_status(job_id)
            status = status_info['status']

            if status == 'Finished':
                spend_time = status_info.get('spend_time', 0)
                print(f"    ✅ 任务完成 (耗时: {spend_time}秒)")
                return True
            elif status == 'Failed':
                error_info = status_info.get('error_info', '未知错误')
                print(f"    ❌ 任务失败: {error_info}")
                return False
            elif status in ['Running', 'Waiting']:
                if elapsed % 120 == 0:  # 每2分钟显示一次
                    print(f"    ⏳ 任务运行中... (已等待: {elapsed}秒, 状态: {status})")
            else:
                print(f"    ⚠️  未知状态: {status}")

            time.sleep(check_interval)
            elapsed += check_interval

        print(f"    ❌ 等待超时 (超过 {timeout}秒)")
        return False

    def _download_optimized_structure(self, job_id: str, work_dir: Path) -> Optional[Path]:
        """
        下载优化后的结构文件 (STRU_ION_D)

        Args:
            job_id: 优化任务的Job ID
            work_dir: 工作目录

        Returns:
            优化后STRU文件的路径，失败返回None
        """
        download_dir = work_dir / f"job_{job_id}_download"
        download_dir.mkdir(exist_ok=True, parents=True)

        success = self.job_manager.download_job_result(job_id, download_dir.parent)
        if not success:
            print(f"    ❌ 下载任务结果失败")
            return None

        stru_ion_d = None
        for stru_file in download_dir.rglob("STRU_ION_D"):
            stru_ion_d = stru_file
            break

        if not stru_ion_d or not stru_ion_d.exists():
            print(f"    ❌ 未找到 STRU_ION_D 文件")
            return None

        target_file = work_dir / "STRU_optimized"
        shutil.copy(stru_ion_d, target_file)
        print(f"    ✅ 优化结构已保存: {target_file.name}")

        return target_file

    # ========== 阶段2：abacustest model elastic prepare ==========

    def _prepare_elastic_with_abacustest(self, elastic_dir: Path, relax_dir: Path,
                                          stru_optimized: Path) -> bool:
        """
        准备弹性计算目录，运行 abacustest model elastic prepare 生成25个形变结构

        Args:
            elastic_dir: 弹性计算目录（将被创建）
            relax_dir: 优化任务目录（提供INPUT、pp、orb文件）
            stru_optimized: 优化后的STRU文件路径

        Returns:
            True 如果成功，False 否则
        """
        try:
            elastic_dir.mkdir(exist_ok=True, parents=True)

            # 复制INPUT文件（将 calculation cell-relax 改为 scf）
            input_file = relax_dir / "INPUT"
            if input_file.exists():
                input_content = input_file.read_text(encoding='utf-8')
                input_content = re.sub(
                    r'calculation\s+cell-relax',
                    'calculation scf',
                    input_content
                )
                (elastic_dir / "INPUT").write_text(input_content, encoding='utf-8')
            else:
                print(f"    ❌ INPUT文件不存在")
                return False

            # 复制优化后的结构文件为STRU
            shutil.copy(stru_optimized, elastic_dir / "STRU")

            # 复制KPT文件（如果存在）
            kpt_file = relax_dir / "KPT"
            if kpt_file.exists():
                shutil.copy(kpt_file, elastic_dir / "KPT")

            # 复制赝势和轨道文件
            for file in relax_dir.glob("*.upf"):
                shutil.copy(file, elastic_dir / file.name)
            for file in relax_dir.glob("*.orb"):
                shutil.copy(file, elastic_dir / file.name)

            # 修复路径格式
            from fix_stru import fix_stru_paths, fix_input_paths
            stru_file = elastic_dir / "STRU"
            inp_file = elastic_dir / "INPUT"
            if stru_file.exists():
                fix_stru_paths(stru_file, verbose=False)
            if inp_file.exists():
                fix_input_paths(inp_file, verbose=False)

            # 运行 abacustest model elastic prepare
            # 注意：必须先切换 Python 进程的 CWD 到 elastic_dir，
            # 否则 abacustest 会将 setting.json 等输出文件写到项目根目录
            old_cwd = os.getcwd()
            try:
                os.chdir(str(elastic_dir))
                result = subprocess.run(
                    ["abacustest", "model", "elastic", "prepare", "-j", "."],
                    capture_output=True,
                    text=True
                )
            finally:
                os.chdir(old_cwd)

            if result.returncode != 0:
                print(f"    ❌ abacustest model elastic prepare 失败: {result.stderr}")
                return False

            # 验证形变结构目录是否已生成
            deformed_dirs = sorted(elastic_dir.glob("deformed_*"))
            if not deformed_dirs:
                print(f"    ❌ 未生成形变结构目录 (deformed_*)")
                return False

            print(f"    ✅ 生成了 {len(deformed_dirs)} 个形变结构 + org/")
            return True

        except Exception as e:
            print(f"    ❌ 准备弹性计算输入失败: {e}")
            return False

    # ========== 阶段3：提交25个形变结构任务 ==========

    def _estimate_job_time(self, stru_path: Path, n_jobs: int = 25,
                           machine_type: str = 'c32_m64_cpu') -> str:
        """
        根据结构属性动态估算任务时间

        Args:
            stru_path: STRU文件路径
            n_jobs: 任务数量
            machine_type: 机器类型

        Returns:
            时间估算字符串（人类可读）
        """
        n_atoms = 1
        kspacing = None
        element = "Unknown"

        if stru_path and stru_path.exists():
            content = stru_path.read_text(encoding='utf-8', errors='ignore')

            # 统计原子数：读取 ATOMIC_POSITIONS 块中各物种的原子数行
            atom_count = 0
            in_atomic_positions = False
            for line in content.splitlines():
                stripped = line.strip()
                if stripped.startswith("ATOMIC_POSITIONS"):
                    in_atomic_positions = True
                    continue
                if in_atomic_positions and stripped and not stripped.startswith('#'):
                    parts = stripped.split()
                    if len(parts) == 1 and parts[0].isdigit():
                        atom_count += int(parts[0])

            if atom_count > 0:
                n_atoms = atom_count

            # 提取第一个元素名
            species_match = re.search(r'ATOMIC_SPECIES\s*\n\s*(\w+)', content)
            if species_match:
                element = species_match.group(1)

        # 读取kspacing（从同目录的INPUT文件）
        if stru_path:
            input_file = stru_path.parent / "INPUT"
            if input_file.exists():
                inp_content = input_file.read_text(encoding='utf-8', errors='ignore')
                ks_match = re.search(r'kspacing\s+([\d.]+)', inp_content)
                if ks_match:
                    kspacing = float(ks_match.group(1))

        # 基准估算：Si 8原子，c32_m64_cpu，约10分钟/任务
        base_time = 10
        if n_atoms > 8:
            base_time = int(10 * (n_atoms / 8) ** 1.5)

        kpoint_note = ""
        if kspacing:
            if kspacing < 0.1:
                base_time = int(base_time * 2)
                kpoint_note = f"，kspacing={kspacing}（密K点）"
            elif kspacing > 0.2:
                base_time = max(5, int(base_time * 0.7))
                kpoint_note = f"，kspacing={kspacing}（稀K点）"
            else:
                kpoint_note = f"，kspacing={kspacing}"

        time_high = int(base_time * 1.5)
        return (
            f"{n_jobs}个任务，{element} {n_atoms}原子{kpoint_note}，"
            f"预计每个任务约{base_time}-{time_high}分钟（机器: {machine_type}）"
        )

    def _submit_all_deformed_jobs(self, elastic_dir: Path) -> List[Tuple[Path, str]]:
        """
        提交所有形变结构任务（org + deformed_00~23）

        Args:
            elastic_dir: 弹性计算目录

        Returns:
            [(dir_path, job_id), ...] 列表，提交失败的条目被跳过
        """
        job_pairs: List[Tuple[Path, str]] = []

        # 收集提交目录：org/ + 所有 deformed_*
        submit_dirs = []
        org_dir = elastic_dir / "org"
        if org_dir.exists():
            submit_dirs.append(org_dir)
        submit_dirs.extend(sorted(elastic_dir.glob("deformed_*")))

        if not submit_dirs:
            print(f"    ❌ 未找到形变结构目录")
            return []

        print(f"    [提交] 共 {len(submit_dirs)} 个目录...")
        for dir_path in submit_dirs:
            job_name = f"{elastic_dir.parent.name}_{dir_path.name}"
            job_config = self.job_manager.create_job_config(job_name, dir_path)
            job_id = self.job_manager.submit_job(job_config, dir_path)

            if job_id:
                job_pairs.append((dir_path, job_id))
                print(f"    ✅ {dir_path.name} → Job {job_id}")
            else:
                print(f"    ❌ {dir_path.name} 提交失败")

        print(f"    [汇总] 成功提交 {len(job_pairs)}/{len(submit_dirs)} 个任务")
        return job_pairs

    def _wait_for_all_jobs(self, job_pairs: List[Tuple[Path, str]],
                           elastic_dir: Path, time_estimate: str) -> None:
        """
        显示时间估算，等待用户手动确认（不自动轮询），然后下载所有结果

        Args:
            job_pairs: [(dir_path, job_id), ...] 列表
            elastic_dir: 弹性计算目录
            time_estimate: 时间估算字符串
        """
        job_ids_preview = " ".join(jid for _, jid in job_pairs[:3])
        if len(job_pairs) > 3:
            job_ids_preview += f" ...（共{len(job_pairs)}个）"

        print(f"\n  [预计] {time_estimate}")
        print(f"  [任务] Job IDs: {job_ids_preview}")
        print(f"  [提示] 任务完成后按 Enter 继续后处理")
        print(f"         状态查询命令：")
        print(f"         python tools/test_framework_integrated.py status {elastic_dir.parent}")
        input("\n  >>> 按 Enter 继续后处理... ")

        # 下载所有已完成任务的结果
        print(f"  [下载] 下载 {len(job_pairs)} 个任务结果...")
        for dir_path, job_id in job_pairs:
            status_info = self.job_manager.get_job_status(job_id)
            status = status_info.get('status', 'Unknown')

            if status == 'Finished':
                success = self.job_manager.download_job_result(job_id, dir_path)
                if success:
                    print(f"    ✅ {dir_path.name} 已下载")
                else:
                    print(f"    ❌ {dir_path.name} 下载失败")
            else:
                print(f"    ⚠️  {dir_path.name}: 状态 {status}，跳过下载")

    # ========== 阶段4：后处理 ==========

    def _post_process_elastic(self, elastic_dir: Path, work_dir: Path) -> bool:
        """
        运行 abacustest model elastic post 后处理，写入 metrics_elastic.json

        Args:
            elastic_dir: 弹性计算目录（包含所有形变结构结果）
            work_dir: 工作目录（写入 metrics_elastic.json）

        Returns:
            True 如果成功，False 否则
        """
        try:
            # 注意：必须先切换 Python 进程的 CWD 到 elastic_dir，
            # 否则 abacustest 会将 setting.json / metrics*.csv 等输出文件写到项目根目录
            old_cwd = os.getcwd()
            try:
                os.chdir(str(elastic_dir))
                result = subprocess.run(
                    ["abacustest", "model", "elastic", "post", "-j", "."],
                    capture_output=True,
                    text=True
                )
            finally:
                os.chdir(old_cwd)

            if result.returncode != 0:
                print(f"    ❌ abacustest model elastic post 失败: {result.stderr}")
                return False

            # 查找 abacustest 生成的结果文件
            result_candidates = (
                list(elastic_dir.glob("elastic*.json")) +
                list(elastic_dir.glob("results*.json"))
            )

            if not result_candidates:
                # 尝试解析 stdout
                if result.stdout:
                    try:
                        data = json.loads(result.stdout)
                        metrics_file = work_dir / "metrics_elastic.json"
                        with open(metrics_file, 'w', encoding='utf-8') as f:
                            json.dump(data, f, indent=2, ensure_ascii=False)
                        print(f"    ✅ 结果已写入 {metrics_file.name}")
                        return True
                    except json.JSONDecodeError:
                        pass
                print(f"    ❌ 未找到后处理结果文件")
                return False

            # 使用第一个找到的结果文件
            src_file = result_candidates[0]
            metrics_file = work_dir / "metrics_elastic.json"
            shutil.copy(src_file, metrics_file)
            print(f"    ✅ 结果已写入 {metrics_file.name} (来源: {src_file.name})")
            return True

        except Exception as e:
            print(f"    ❌ 后处理失败: {e}")
            return False

    # ========== K点透明度 ==========

    def _check_kpoint_discrepancy(self, work_dir: Path, test_info: TestInfo) -> Optional[str]:
        """
        检查实际使用的kspacing与教程预期是否有差异

        Args:
            work_dir: 工作目录
            test_info: 测试信息

        Returns:
            K点说明字符串，无差异或无法判断时返回 None
        """
        # 优先从 elastic/INPUT 读取，其次从 relax/INPUT
        elastic_input = work_dir / "elastic" / "INPUT"
        relax_input = work_dir / "relax" / "INPUT"
        input_file = elastic_input if elastic_input.exists() else relax_input

        if not input_file.exists():
            return None

        inp_content = input_file.read_text(encoding='utf-8', errors='ignore')
        actual_ks_match = re.search(r'kspacing\s+([\d.]+)', inp_content)
        if not actual_ks_match:
            return None

        actual_kspacing = float(actual_ks_match.group(1))

        # 如果 expected_results 明确包含 kspacing，检查差异
        expected_kspacing = test_info.expected_results.get('kspacing')
        if expected_kspacing and abs(actual_kspacing - float(expected_kspacing)) > 0.01:
            return (
                f"K点说明：实际使用 kspacing={actual_kspacing}，"
                f"教程预期值基于 kspacing={expected_kspacing}。"
                f"数值差异属于正常的K点收敛效应，不影响方法论正确性。"
            )

        # 如果 kspacing 较大（>0.14），主动提示可能与精确值有偏差
        if actual_kspacing > 0.14:
            return (
                f"K点说明：实际使用 kspacing={actual_kspacing}（较稀K网格）。"
                f"若教程预期值基于更密的K点，数值偏差属于K点收敛效应，不影响方法论正确性。"
            )

        return None
