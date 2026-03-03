"""
AutoTutorial 3.0 - 隐式溶剂模型测试插件
处理含 imp_sol = 1 的 PW-SCF 计算（H₂ 分子在水溶液中）
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────
# 硬编码的 H₂@水溶液 算例文件内容
# ─────────────────────────────────────────────

# H₂ 分子置于 10 Å 立方超胞中心
# 晶格常数 = 1.0 Bohr（使用直接矢量指定 Bohr 单位晶格）
# H-H 键长 0.74 Å = 1.3984 Bohr，分数坐标各 ±0.037
# 超胞 18.8973 Bohr（= 10 Å），H-H 分数偏移 = 0.6992/18.8973 ≈ 0.037
_H2_STRU_TEMPLATE = """\
ATOMIC_SPECIES
H 1.008 H_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
18.8973  0.0000  0.0000
 0.0000 18.8973  0.0000
 0.0000  0.0000 18.8973

ATOMIC_POSITIONS
Direct

H
0.0
2
0.500000  0.500000  0.463000  0  0  0
0.500000  0.500000  0.537000  0  0  0
"""

# INPUT 来自教程案例（完整，所有参数原样保留）
_H2_INPUT_TEMPLATE = """\
INPUT_PARAMETERS
#Parameters (1.General)
suffix              H2
calculation         scf
ntype               1
nbands              2
symmetry            0
pseudo_dir          ./

#Parameters (2.Iteration)
ecutwfc             60
scf_thr             1e-6
scf_nmax            100

#Parameters (3.Basis)
basis_type          pw

#Parameters (Solvation Model)
imp_sol             1
eb_k                80
tau                 0.000010798
sigma_k             0.6
nc_k                0.00037
"""

_H2_KPT = """\
K_POINTS
0
Gamma
1 1 1 0 0 0
"""


class SolvationPlugin(BaseTestPlugin):
    """隐式溶剂模型 PW-SCF 测试插件"""

    def __init__(self, job_manager=None, pp_manager=None):
        super().__init__(job_manager, pp_manager)

    @property
    def plugin_name(self) -> str:
        return "隐式溶剂模型 SCF 测试"

    @property
    def calc_type(self) -> str:
        return "solvation"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测是否包含 imp_sol 隐式溶剂相关参数"""
        patterns = [
            r'imp_sol\s*[=\s]+1',
            r'隐式溶剂',
            r'Implicit.*Solvation',
            r'implicit.*solvation',
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """从教程提取隐式溶剂 SCF 测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="H2"
        )

        # ── 尝试从教程提取 INPUT，失败则使用模板 ──
        extracted_input = self._extract_input(tutorial_content)
        if extracted_input and 'imp_sol' in extracted_input:
            test_info.input_content = extracted_input
            print("    [INPUT] 从教程提取成功（含 imp_sol 参数）")
        else:
            test_info.input_content = _H2_INPUT_TEMPLATE
            print("    [INPUT] 使用内置 H₂ 隐式溶剂模板")

        # ── STRU：H₂ 分子置于 10 Å 立方超胞（教程未提供，使用近似结构）──
        test_info.stru_content = _H2_STRU_TEMPLATE
        print("    [STRU] H₂ 分子，10 Å 立方超胞，H-H 键长 0.74 Å（近似结构）")

        # ── KPT：Γ 点（分子体系）──
        test_info.kpt_content = _H2_KPT
        print("    [KPT] Γ 点（分子体系，无需 k 点采样）")

        # ── 预期结果：定性验证（符号检查）──
        # 教程仅给出定性期望：E_sol_el < 0，E_sol_cav > 0
        # 使用特殊标记值 -1/+1 表示符号约束（validate_results 中特殊处理）
        test_info.expected_results = {
            'e_sol_el_negative': -1.0,   # 期望 E_sol_el < 0
            'e_sol_cav_positive': 1.0,   # 期望 E_sol_cav > 0
            'scf_converged': 1.0,        # 期望 SCF 收敛
        }

        # ── 赝势（PW 基组，无轨道文件）──
        test_info.pseudopotentials = ['H_ONCV_PBE-1.0.upf']
        test_info.orbitals = []

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备 H₂ 隐式溶剂 SCF 输入文件"""
        input_dir = work_dir / f"{test_info.case_name}_solvation"
        input_dir.mkdir(exist_ok=True)

        # 写入输入文件
        (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        print(f"  [写入] INPUT / STRU / KPT → {input_dir}")

        # 下载赝势文件（PW 基组无轨道文件）
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                try:
                    pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                    shutil.copy(pp_file, input_dir / pp)
                    print(f"  [PP]  {pp}")
                except Exception as e:
                    print(f"  [WARN] 赝势下载失败 {pp}: {e}")

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交隐式溶剂 SCF 任务"""
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

    def validate_results(self, job_ids: List[str], work_dir: Path,
                         test_info: TestInfo) -> ValidationResult:
        """验证隐式溶剂 SCF 结果（定性验证：符号检查）"""
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
        actual = self._extract_actual_results(results_dir)

        # ── 验证 SCF 收敛 ──
        scf_conv = actual.get('scf_converged', 0.0)
        validation.comparisons['scf_converged'] = {
            'key': 'scf_converged',
            'expected': 1.0,
            'actual': scf_conv,
            'abs_error': abs(scf_conv - 1.0),
            'rel_error': abs(scf_conv - 1.0),
            'passed': scf_conv >= 1.0,
            'use_abs': True,
            'description': 'SCF 是否收敛'
        }

        # ── 验证 E_sol_el < 0（静电项为负）──
        e_sol_el = actual.get('e_sol_el', None)
        if e_sol_el is not None:
            passed_el = e_sol_el < 0
            validation.comparisons['e_sol_el_negative'] = {
                'key': 'e_sol_el_negative',
                'expected': -1.0,    # 符号标记
                'actual': e_sol_el,
                'abs_error': 0.0 if passed_el else abs(e_sol_el),
                'rel_error': 0.0 if passed_el else 1.0,
                'passed': passed_el,
                'use_abs': True,
                'description': f'E_sol_el 应 < 0（实际: {e_sol_el:.6f} eV）'
            }
        else:
            validation.warnings.append("未找到 E_sol_el 输出值")

        # ── 验证 E_sol_cav > 0（空腔项为正）──
        e_sol_cav = actual.get('e_sol_cav', None)
        if e_sol_cav is not None:
            passed_cav = e_sol_cav > 0
            validation.comparisons['e_sol_cav_positive'] = {
                'key': 'e_sol_cav_positive',
                'expected': 1.0,    # 符号标记
                'actual': e_sol_cav,
                'abs_error': 0.0 if passed_cav else abs(e_sol_cav),
                'rel_error': 0.0 if passed_cav else 1.0,
                'passed': passed_cav,
                'use_abs': True,
                'description': f'E_sol_cav 应 > 0（实际: {e_sol_cav:.6f} eV）'
            }
        else:
            validation.warnings.append("未找到 E_sol_cav 输出值")

        # ── 汇总 ──
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的结果（日志文件可能不存在）")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"
        report += "> 验证方式：定性符号检查（E_sol_el < 0，E_sol_cav > 0，SCF 收敛）\n\n"

        if validation.comparisons:
            report += "| 验证项 | 期望 | 实际值 | 状态 |\n"
            report += "|--------|------|--------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                desc = comp.get('description', key)
                actual_str = f"{comp['actual']:.6f}" if isinstance(comp['actual'], float) else str(comp['actual'])
                report += f"| {desc} | — | {actual_str} | {status} |\n"

        if validation.errors:
            report += "\n**错误：**\n"
            for err in validation.errors:
                report += f"- {err}\n"

        if validation.warnings:
            report += "\n**警告：**\n"
            for warn in validation.warnings:
                report += f"- {warn}\n"

        overall = "✅ 通过" if validation.passed else "❌ 失败"
        report += f"\n**总体结果：** {overall}\n\n"
        return report

    # ── 辅助方法 ──────────────────────────────

    def _extract_input(self, content: str) -> Optional[str]:
        """从教程提取 INPUT 文件内容"""
        pattern = r'```[^\n]*\n(INPUT_PARAMETERS[\s\S]+?)```'
        match = re.search(pattern, content)
        return match.group(1).strip() if match else None

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从 ABACUS 输出文件提取溶剂化计算结果"""
        results = {}

        log_files = list(results_dir.rglob("running_scf.log"))
        if not log_files:
            log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            text = log_file.read_text(encoding='utf-8', errors='ignore')

            # 检查 SCF 收敛
            if 'convergence is achieved' in text or 'SCF CONVERGED' in text.upper():
                results['scf_converged'] = 1.0
            else:
                results['scf_converged'] = 0.0

            # 提取 E_sol_el（静电项）
            # ABACUS 输出格式：E_sol_el = -X.XXXXXX eV  或  E_sol_el      -X.XXXXXX
            m = re.search(r'E_sol_el\s*[=:]?\s*([-+]?\d+\.\d+)', text)
            if m:
                results['e_sol_el'] = float(m.group(1))

            # 提取 E_sol_cav（空腔项）
            m = re.search(r'E_sol_cav\s*[=:]?\s*([-+]?\d+\.\d+)', text)
            if m:
                results['e_sol_cav'] = float(m.group(1))

            # 总能量（用于完整性记录）
            m = re.search(r'FINAL_ETOT_IS\s+([-\d.]+)\s+eV', text)
            if m:
                results['total_energy'] = float(m.group(1))

            if results:
                break

        return results
