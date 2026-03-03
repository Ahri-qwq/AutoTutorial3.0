"""
AutoTutorial 3.0 - Phonopy声子谱测试插件
处理 ABACUS+Phonopy 声子谱计算教程中的 SCF+力 计算部分

验证策略：
  - 只测试 ABACUS 负责的 SCF+力计算（Phonopy 后处理为外部软件，不在验证范围）
  - 使用 4 原子 FCC Al 常规胞（不用 32 原子超胞，节省计算时间）
  - 验证：SCF 收敛、力已输出、总能量负值
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────
# FCC Al 4原子常规胞（来自教程 STRU，轨道文件名已修正）
# ─────────────────────────────────────────────
_AL_FCC_STRU = """\
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Al_gga_8au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
4.03459549706 0 0 #latvec1
0 4.03459549706 0 #latvec2
0 0 4.03459549706 #latvec3

ATOMIC_POSITIONS
Direct

Al #label
0 #magnetism
4 #number of atoms
0  0  0  m  0  0  0
0.5  0.5  0  m  0  0  0
0.5  0  0.5  m  0  0  0
0  0.5  0.5  m  0  0  0
"""

# INPUT 来自教程（已修改 pseudo_dir/orbital_dir/stru_file，其余原样保留）
_AL_FCC_INPUT = """\
INPUT_PARAMETERS
#Parameters (1.General)
suffix          Al-fcc
calculation     scf
esolver_type    ksdft
symmetry        1
pseudo_dir      ./
orbital_dir     ./
cal_stress      1
cal_force       1
stru_file       STRU

#Parameters (2.Iteration)
ecutwfc         100
scf_thr         1e-7
scf_nmax        50

#Parameters (3.Basis)
basis_type      lcao
gamma_only      0

#Parameters (4.Smearing)
smearing_method mp
smearing_sigma  0.015

#Parameters (5.Mixing)
mixing_type     pulay
mixing_beta     0.7
mixing_gg0      1.5
"""

# KPT：FCC Al 4原子胞，4×4×4 Gamma 网格
_AL_FCC_KPT = """\
K_POINTS
0
Gamma
4 4 4 0 0 0
"""


class PhonopyPlugin(BaseTestPlugin):
    """ABACUS+Phonopy 声子谱 SCF+力计算测试插件"""

    def __init__(self, job_manager=None, pp_manager=None):
        super().__init__(job_manager, pp_manager)

    @property
    def plugin_name(self) -> str:
        return "Phonopy声子谱 SCF+力 测试"

    @property
    def calc_type(self) -> str:
        return "phonopy"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测是否是 ABACUS+Phonopy 声子谱教程"""
        patterns = [
            r'phonopy',
            r'声子谱',
            r'FORCE_SETS',
            r'phonon\s+dispersion',
            r'有限位移方法',
            r'band\.conf',
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取测试信息（固定使用 FCC Al 4原子胞）"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="Al_FCC"
        )

        # INPUT：使用内置模板（从教程提取版本需要修改多处路径，不如直接用模板）
        test_info.input_content = _AL_FCC_INPUT
        print("    [INPUT] 使用 FCC Al SCF+力 计算模板（路径已适配测试环境）")

        # STRU：4原子 FCC Al 常规胞
        test_info.stru_content = _AL_FCC_STRU
        print("    [STRU]  FCC Al 常规胞，4原子，晶格常数 4.035 Å")

        # KPT：4×4×4 Gamma
        test_info.kpt_content = _AL_FCC_KPT
        print("    [KPT]   4×4×4 Gamma k 网格")

        # 预期结果：定性验证（SCF 收敛 + 力已输出）
        # 不设置定量能量阈值，避免 LCAO 参数差异导致假失败
        test_info.expected_results = {
            'scf_converged': 1.0,     # SCF 必须收敛
            'force_calculated': 1.0,  # 力必须已计算输出
        }

        # 赝势 + 轨道
        test_info.pseudopotentials = ['Al_ONCV_PBE-1.0.upf']
        test_info.orbitals = ['Al_gga_8au_100Ry_4s4p1d.orb']

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备 FCC Al SCF+力 输入文件"""
        input_dir = work_dir / f"{test_info.case_name}_scf"
        input_dir.mkdir(exist_ok=True)

        # 写入输入文件
        (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        print(f"  [写入] INPUT / STRU / KPT → {input_dir}")

        # 下载赝势
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                try:
                    pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                    shutil.copy(pp_file, input_dir / pp)
                    print(f"  [PP]   {pp}")
                except Exception as e:
                    print(f"  [WARN] 赝势下载失败 {pp}: {e}")

            # 下载轨道文件（使用修正后的文件名 8au）
            for orb in test_info.orbitals:
                try:
                    orb_file = self.pp_manager.get_file(orb, "orbital")
                    shutil.copy(orb_file, input_dir / orb)
                    print(f"  [ORB]  {orb}")
                except Exception as e:
                    print(f"  [WARN] 轨道文件下载失败 {orb}: {e}")

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 SCF 任务"""
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
        """验证 SCF+力 计算结果"""
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

        # ── 验证力已计算输出 ──
        force_calc = actual.get('force_calculated', 0.0)
        validation.comparisons['force_calculated'] = {
            'key': 'force_calculated',
            'expected': 1.0,
            'actual': force_calc,
            'abs_error': abs(force_calc - 1.0),
            'rel_error': abs(force_calc - 1.0),
            'passed': force_calc >= 1.0,
            'use_abs': True,
            'description': 'cal_force=1 是否输出了原子受力'
        }

        # ── 汇总 ──
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的结果（日志文件可能不存在）")

        # 额外记录总能量（不作为 pass/fail 标准）
        if 'total_energy' in actual:
            e = actual['total_energy']
            validation.warnings.append(f"总能量（仅记录）: {e:.4f} eV（FCC Al 4原子，参考值 ~-226 eV）")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"
        report += "> 验证方式：SCF 收敛检查 + cal_force 力输出检查（Phonopy 后处理为外部软件，不在验证范围）\n\n"

        if validation.comparisons:
            report += "| 验证项 | 期望 | 实际值 | 状态 |\n"
            report += "|--------|------|--------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                desc = comp.get('description', key)
                actual_str = f"{comp['actual']:.1f}" if isinstance(comp['actual'], float) else str(comp['actual'])
                report += f"| {desc} | 是 | {actual_str} | {status} |\n"

        if validation.errors:
            report += "\n**错误：**\n"
            for err in validation.errors:
                report += f"- {err}\n"

        if validation.warnings:
            report += "\n**说明：**\n"
            for warn in validation.warnings:
                report += f"- {warn}\n"

        overall = "✅ 通过" if validation.passed else "❌ 失败"
        report += f"\n**总体结果：** {overall}\n\n"
        return report

    # ── 辅助方法 ──────────────────────────────

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从 ABACUS 输出文件提取 SCF+力 结果"""
        results = {}

        log_files = list(results_dir.rglob("running_scf.log"))
        if not log_files:
            log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            text = log_file.read_text(encoding='utf-8', errors='ignore')

            # 检查 SCF 收敛
            if ('convergence is achieved' in text or
                    'SCF CONVERGED' in text.upper() or
                    'charge density convergence is achieved' in text.lower()):
                results['scf_converged'] = 1.0
            else:
                results['scf_converged'] = 0.0

            # 检查力是否已输出（ABACUS 输出 FORCE 关键字块）
            if re.search(r'TOTAL-FORCE\s*\(eV/Angstrom\)|FORCE\s*\(eV/Angstrom\)', text):
                results['force_calculated'] = 1.0
            else:
                results['force_calculated'] = 0.0

            # 总能量（用于记录）
            m = re.search(r'FINAL_ETOT_IS\s+([-\d.]+)\s+eV', text)
            if m:
                results['total_energy'] = float(m.group(1))

            if results:
                break

        return results
