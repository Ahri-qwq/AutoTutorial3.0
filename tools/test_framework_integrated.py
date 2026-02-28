"""
AutoTutorial 3.0 - 完整测试流程
整合所有Phase，支持插件化测试
"""

import sys
import io
# Fix Unicode output on Windows (GBK terminal)
if hasattr(sys.stdout, 'buffer'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
if hasattr(sys.stderr, 'buffer'):
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from pathlib import Path
from typing import Dict, List
import subprocess
import json
import time
from datetime import datetime, timedelta

# 导入所有模块
from test_framework_phase1_analyzer import TutorialAnalyzer, TutorialAnalysis
from test_framework_phase3_7_impl import (
    BohriumJobManager,
    PseudopotentialManager,
    ResultComparator,
    ReportGenerator
)

# 导入插件系统
from test_plugins.base_plugin import BaseTestPlugin
from test_plugins.relax_plugin import RelaxPlugin
from test_plugins.elastic_plugin import ElasticPlugin
from test_plugins.band_plugin import BandPlugin
from test_plugins.dos_plugin import DOSPlugin
from test_plugins.dftu_plugin import DFTUPlugin


class FullTestExecutor:
    """完整测试执行器 - 支持插件化测试"""

    def __init__(self, tutorial_path: str, work_dir: str = None, test_dir: str = None):
        self.tutorial_path = Path(tutorial_path).resolve() if tutorial_path else Path(".")

        if test_dir:
            # 直接使用外部传入的目录（testCLAUDE.md 新流程）
            self.test_dir = Path(test_dir).resolve()
            self.test_dir.mkdir(exist_ok=True, parents=True)
        else:
            # 兼容旧流程：在 work_dir 下自动创建时间戳子目录
            wd = Path(work_dir or "./test_workspace").resolve()
            wd.mkdir(exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.test_dir = wd / f"{timestamp}_{self.tutorial_path.stem}"
            self.test_dir.mkdir(exist_ok=True)

        self.analysis: TutorialAnalysis = None
        self.start_time = None
        self.end_time = None

        # 初始化管理器
        self.job_manager = BohriumJobManager()
        self.pp_manager = PseudopotentialManager()
        self.comparator = ResultComparator(tolerance=0.05)
        self.report_generator = ReportGenerator(self.test_dir)

        # 初始化插件系统
        self.plugins: List[BaseTestPlugin] = [
            RelaxPlugin(self.job_manager, self.pp_manager),
            ElasticPlugin(self.job_manager, self.pp_manager),
            BandPlugin(self.job_manager, self.pp_manager),
            DOSPlugin(self.job_manager, self.pp_manager),
            DFTUPlugin(self.job_manager, self.pp_manager),
        ]

        # 存储每个插件的测试信息和任务
        self.plugin_tests: Dict[str, Dict] = {}

    def run_full_test(self):
        """运行完整测试流程"""
        self.start_time = datetime.now()

        print("\n" + "=" * 70)
        print("AutoTutorial 3.0 - 自动化测试流程")
        print("=" * 70)

        try:
            # Phase 1: 分析文章
            self.phase1_analyze()

            # Phase 2: 准备输入文件
            self.phase2_prepare_inputs()

            # Phase 3: 提交任务
            self.phase3_submit_jobs()

            # Phase 4: 智能监控（1分钟后检查）
            self.phase4_smart_monitor()

            # 等待用户反馈
            print("\n" + "=" * 70)
            print("[等待] 请稍后检查任务状态")
            print("=" * 70)
            print("\n提示：当任务完成后，运行以下命令继续：")
            print(f"  python tools/test_framework_continue.py {self.test_dir}")

            # 保存状态
            self._save_state()

        except Exception as e:
            print(f"\n[ERROR] 测试失败: {e}")
            import traceback
            traceback.print_exc()

    def run_prepare_only(self):
        """仅运行 Phase 1+2：解析教程 + 准备输入文件（不提交、不监控）"""
        self.start_time = datetime.now()

        print("\n" + "=" * 70)
        print("AutoTutorial 3.0 - 准备输入文件（--phase prepare）")
        print("=" * 70)

        try:
            # Phase 1: 分析文章
            self.phase1_analyze()

            # Phase 2: 准备输入文件
            self.phase2_prepare_inputs()

            print("\n" + "=" * 70)
            print("[OK] 准备完成，请按教程内容决定任务提交顺序")
            print("=" * 70)
            print(f"\n测试目录：{self.test_dir}")
            print(f"解析结果：{self.test_dir / '01_analysis.json'}")
            print("\n下一步：")
            print("  1. 读取教程内容，理解计算流程和任务依赖关系")
            print("  2. 按顺序用 bohr job submit 提交任务（串行/并行由你决定）")
            print(f"  3. 所有任务完成后运行：")
            print(f"     python tools/test_framework_integrated.py continue {self.test_dir}")

        except Exception as e:
            print(f"\n[ERROR] 准备失败: {e}")
            import traceback
            traceback.print_exc()

    def continue_after_jobs_done(self):
        """任务完成后继续执行"""
        print("\n" + "=" * 70)
        print("继续测试流程")
        print("=" * 70)

        try:
            # 加载状态
            self._load_state()

            # Phase 5: 下载结果
            self.phase5_download_results()

            # Phase 6: 对比验证
            comparison = self.phase6_compare_results()

            # Phase 7: 生成报告
            self.end_time = datetime.now()
            self.phase7_generate_report(comparison)

            print("\n" + "=" * 70)
            print("[OK] 测试完成！")
            print("=" * 70)

        except Exception as e:
            print(f"\n[ERROR] 测试失败: {e}")
            import traceback
            traceback.print_exc()

    def phase1_analyze(self):
        """Phase 1: 分析文章并识别测试类型"""
        print("\n" + "-" * 70)
        print("Phase 1: 分析文章")
        print("-" * 70)

        # 读取教程内容
        tutorial_content = self.tutorial_path.read_text(encoding='utf-8')

        # 识别可以处理的插件
        print("\n[检测] 识别计算类型...")
        active_plugins = []
        for plugin in self.plugins:
            if plugin.can_handle(tutorial_content):
                print(f"  [OK] 检测到: {plugin.plugin_name}")
                active_plugins.append(plugin)

        if not active_plugins:
            print("  [WARNING] 未检测到支持的计算类型")
            return

        # 为每个插件提取测试信息
        print("\n[提取] 提取测试信息...")
        for plugin in active_plugins:
            test_info = plugin.extract_test_info(tutorial_content)
            if test_info:
                self.plugin_tests[plugin.calc_type] = {
                    'plugin': plugin,
                    'test_info': test_info,
                    'job_ids': [],
                    'validation': None
                }
                print(f"  [OK] {plugin.plugin_name}: 案例 {test_info.case_name}")
            else:
                print(f"  [ERROR] {plugin.plugin_name}: 提取失败")

        # 保存分析结果
        analysis_file = self.test_dir / "01_analysis.json"
        with open(analysis_file, 'w', encoding='utf-8') as f:
            json.dump({
                'detected_types': [pt['test_info'].calc_type for pt in self.plugin_tests.values()],
                'case_names': [pt['test_info'].case_name for pt in self.plugin_tests.values()]
            }, f, indent=2, ensure_ascii=False)

        print(f"\n[OK] 分析结果已保存: {analysis_file}")

    def phase2_prepare_inputs(self):
        """Phase 2: 准备输入文件"""
        print("\n" + "-" * 70)
        print("Phase 2: 准备输入文件")
        print("-" * 70)

        for calc_type, test_data in self.plugin_tests.items():
            plugin = test_data['plugin']
            test_info = test_data['test_info']

            print(f"\n[准备] {plugin.plugin_name} - {test_info.case_name}")
            input_dirs = plugin.prepare_inputs(test_info, self.test_dir)

            if input_dirs:
                test_data['input_dirs'] = input_dirs
                print(f"  [OK] 准备完成: {len(input_dirs)} 个目录")
            else:
                print(f"  [ERROR] 准备失败")

        print("\n[OK] 输入文件准备完成")

    def phase3_submit_jobs(self):
        """Phase 3: 提交任务"""
        print("\n" + "-" * 70)
        print("Phase 3: 提交任务")
        print("-" * 70)

        for calc_type, test_data in self.plugin_tests.items():
            if 'input_dirs' not in test_data:
                continue

            plugin = test_data['plugin']
            input_dirs = test_data['input_dirs']

            print(f"\n[提交] {plugin.plugin_name}")
            job_ids = plugin.submit_jobs(input_dirs)

            if job_ids:
                test_data['job_ids'] = job_ids
                print(f"  [OK] 提交成功: {len(job_ids)} 个任务")
            else:
                print(f"  [ERROR] 提交失败")

        print("\n[OK] 任务提交完成")

    def phase4_smart_monitor(self):
        """Phase 4: 智能监控"""
        print("\n" + "-" * 70)
        print("Phase 4: 智能监控")
        print("-" * 70)

        # 收集所有 Job IDs
        all_job_ids = []
        for test_data in self.plugin_tests.values():
            all_job_ids.extend(test_data.get('job_ids', []))

        if not all_job_ids:
            print("[跳过] 没有任务需要监控")
            return

        print("[等待] 等待1分钟后检查任务状态...")
        time.sleep(60)

        print("\n[检查] 任务状态...")
        status_summary = self.job_manager.check_jobs_running(all_job_ids)

        print(f"  - Running: {status_summary['Running']}")
        print(f"  - Finished: {status_summary['Finished']}")
        print(f"  - Failed: {status_summary['Failed']}")

        if status_summary['Running'] > 0:
            # 估算完成时间（假设每个任务12分钟）
            estimated_minutes = 12
            estimated_time = datetime.now() + timedelta(minutes=estimated_minutes)
            print(f"\n[OK] 任务正常运行")
            print(f"[预计] 完成时间: {estimated_time.strftime('%H:%M:%S')}")
        elif status_summary['Failed'] > 0:
            print(f"\n[ERROR] 有任务失败，请检查")
        elif status_summary['Finished'] == len(all_job_ids):
            print(f"\n[OK] 所有任务已完成")

    def phase5_download_results(self):
        """Phase 5: 下载结果"""
        print("\n" + "-" * 70)
        print("Phase 5: 下载结果")
        print("-" * 70)

        results_dir = self.test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        for calc_type, test_data in self.plugin_tests.items():
            job_ids = test_data.get('job_ids', [])
            if not job_ids:
                continue

            plugin = test_data['plugin']
            print(f"\n[下载] {plugin.plugin_name}")

            for job_id in job_ids:
                print(f"  下载 Job {job_id}...")
                self.job_manager.download_job_result(job_id, results_dir)

        print("\n[OK] 结果下载完成")

    def phase6_compare_results(self) -> Dict:
        """Phase 6: 对比验证"""
        print("\n" + "-" * 70)
        print("Phase 6: 对比验证")
        print("-" * 70)

        all_validations = {}

        for calc_type, test_data in self.plugin_tests.items():
            job_ids = test_data.get('job_ids', [])
            if not job_ids:
                continue

            plugin = test_data['plugin']
            test_info = test_data['test_info']

            print(f"\n[验证] {plugin.plugin_name}")
            validation = plugin.validate_results(job_ids, self.test_dir, test_info)

            test_data['validation'] = validation
            all_validations[calc_type] = validation

            status = "[PASS]" if validation.passed else "[FAIL]"
            print(f"  结果: {status}")

        print("\n[OK] 验证完成")
        return all_validations

    def phase7_generate_report(self, validations: Dict):
        """Phase 7: 生成测试报告"""
        print("\n" + "-" * 70)
        print("Phase 7: 生成测试报告")
        print("-" * 70)

        report_file = self.test_dir.parent / "test_report.md"

        # 生成报告头部
        report = "# AutoTutorial 3.0 - 测试报告\n\n"
        report += f"**测试时间：** {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        report += f"**教程文件：** {self.tutorial_path.name}\n\n"
        report += "---\n\n"

        # 生成每个插件的报告章节
        report += "## 测试结果\n\n"
        for calc_type, validation in validations.items():
            test_data = self.plugin_tests[calc_type]
            plugin = test_data['plugin']
            report += plugin.generate_report_section(validation)

        # 生成总结
        report += "## 总结\n\n"
        total_tests = len(validations)
        passed_tests = sum(1 for v in validations.values() if v.passed)
        report += f"- 总测试数：{total_tests}\n"
        report += f"- 通过数：{passed_tests}\n"
        report += f"- 失败数：{total_tests - passed_tests}\n"
        report += f"- 通过率：{passed_tests/total_tests*100:.1f}%\n\n"

        if passed_tests == total_tests:
            report += "[OK] **所有测试通过！**\n"
        else:
            report += "[FAIL] **部分测试失败，请检查详细结果。**\n"

        # 保存报告
        report_file.write_text(report, encoding='utf-8')
        print(f"\n[OK] 测试报告已生成: {report_file}")

    def _save_state(self):
        """保存测试状态"""
        state_file = self.test_dir / "test_state.json"
        plugin_states = {}
        for calc_type, test_data in self.plugin_tests.items():
            plugin_states[calc_type] = {
                'job_ids': test_data.get('job_ids', []),
                'case_name': test_data['test_info'].case_name
            }

        with open(state_file, 'w', encoding='utf-8') as f:
            json.dump({
                'tutorial_path': str(self.tutorial_path),
                'plugin_states': plugin_states,
                'start_time': self.start_time.isoformat() if self.start_time else None,
            }, f, indent=2)

    def _load_state(self):
        """加载测试状态"""
        state_file = self.test_dir / "test_state.json"
        if state_file.exists():
            with open(state_file, 'r', encoding='utf-8') as f:
                state = json.load(f)
                plugin_states = state.get('plugin_states', {})
                # 重新加载插件状态
                for calc_type, pstate in plugin_states.items():
                    if calc_type in self.plugin_tests:
                        self.plugin_tests[calc_type]['job_ids'] = pstate.get('job_ids', [])
                if state.get('start_time'):
                    self.start_time = datetime.fromisoformat(state['start_time'])


# 使用示例
if __name__ == "__main__":
    # 手动解析 --test-dir 和 --phase 参数（避免引入额外依赖）
    test_dir_arg = None
    phase_arg = None
    remaining = []
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '--test-dir' and i + 1 < len(sys.argv):
            test_dir_arg = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--phase' and i + 1 < len(sys.argv):
            phase_arg = sys.argv[i + 1]
            i += 2
        else:
            remaining.append(sys.argv[i])
            i += 1

    if remaining and remaining[0] == "continue":
        # 继续模式：python tools/test_framework_integrated.py continue <test_dir>
        td = remaining[1] if len(remaining) > 1 else test_dir_arg
        if td:
            executor = FullTestExecutor("", test_dir=td)
            executor.continue_after_jobs_done()
        else:
            print("用法: python tools/test_framework_integrated.py continue <test_dir>")
    else:
        # 正常模式：python tools/test_framework_integrated.py <tutorial_path> [--test-dir <dir>] [--phase prepare]
        tutorial_path = remaining[0] if remaining else "_workspace/elastic_tutorial.md"
        executor = FullTestExecutor(
            tutorial_path=tutorial_path,
            test_dir=test_dir_arg
        )
        if phase_arg == "prepare":
            executor.run_prepare_only()
        else:
            executor.run_full_test()
