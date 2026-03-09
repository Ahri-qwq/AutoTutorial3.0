"""
AutoTutorial 3.0 - NEB/AutoNEB 过渡态搜索测试插件
处理 ATST-Tools + ABACUS 的 NEB 计算教程验证

测试策略：
- NEB 全过程（多步迭代 + 多ABACUS调用）成本极高，不直接运行
- 改为对 Case 1（Li in Si）的初始态做单点 LCAO-SCF，验证：
    1. 赝势（Li_ONCV_PBE-1.2.upf, Si_ONCV_PBE-1.2.upf）可用
    2. 轨道文件（Li_gga_8au_100Ry_4s1p.orb, Si_gga_8au_100Ry_2s2p1d.orb）可用
    3. ABACUS 参数（basis_type=lcao, kpts=[2,2,2], nspin=1, symmetry=0）有效
- Case 2（AutoNEB Cy-Pt@graphene）因 nspin=2 + vdw + dipole + 10-image 过于昂贵，跳过
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────
# Li in Si 初始态近似结构
# 使用 Si 金刚石结构常规晶胞（8原子）+ Li 填隙
# Si 金刚石：a = 5.431 Å = 10.262 Bohr
# Li 置于八面体间隙 (0.5, 0.5, 0.5)（直角坐标分数位）
# 实际 NEB 算例使用更大超胞，此处仅用于验证参数可用性
# ─────────────────────────────────────────────

_LI_SI_STRU = """\
ATOMIC_SPECIES
Si  28.0855  Si_ONCV_PBE-1.2.upf
Li   6.941   Li_ONCV_PBE-1.2.upf

NUMERICAL_ORBITAL
Si_gga_8au_100Ry_2s2p1d.orb
Li_gga_8au_100Ry_4s1p.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
10.262  0.000  0.000
 0.000 10.262  0.000
 0.000  0.000 10.262

ATOMIC_POSITIONS
Direct

Si
0.0
8
0.000  0.000  0.000  0  0  0
0.500  0.500  0.000  0  0  0
0.500  0.000  0.500  0  0  0
0.000  0.500  0.500  0  0  0
0.250  0.250  0.250  0  0  0
0.750  0.750  0.250  0  0  0
0.750  0.250  0.750  0  0  0
0.250  0.750  0.750  0  0  0

Li
0.0
1
0.500  0.500  0.500  0  0  0
"""

# INPUT：参照 NEB 教程中 Li-Si 案例的 ABACUS 参数
# 去除 NEB 专用参数（k, fmax 等），保留 ABACUS 侧参数
_LI_SI_INPUT = """\
INPUT_PARAMETERS
#Parameters (1.General)
suffix              Li-Si-NEB-test
calculation         scf
ntype               2
nspin               1
symmetry            0
pseudo_dir          ./
orbital_dir         ./

#Parameters (2.Iteration)
ecutwfc             100
scf_thr             1e-6
scf_nmax            100

#Parameters (3.Basis)
basis_type          lcao
"""

_LI_SI_KPT = """\
K_POINTS
0
Gamma
2 2 2 0 0 0
"""


class NEBPlugin(BaseTestPlugin):
    """NEB/AutoNEB 过渡态搜索测试插件（ATST-Tools + ABACUS）

    实际测试内容：Li in Si IS 单点 SCF（验证赝势/轨道/参数可用）
    """

    def __init__(self, job_manager=None, pp_manager=None):
        super().__init__(job_manager, pp_manager)

    @property
    def plugin_name(self) -> str:
        return "NEB/AutoNEB 过渡态搜索测试"

    @property
    def calc_type(self) -> str:
        return "neb"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测是否为 NEB/ATST-Tools 相关教程"""
        patterns = [
            r'AbacusNEB',
            r'AbacusAutoNEB',
            r'ATST-Tools',
            r'atst.tools',
            r'dyneb_run',
            r'autoneb_run',
            r'neb_run\.py',
            r'AutoNEB.*ABACUS',
            r'ABACUS.*AutoNEB',
            r'NEB.*ATST',
            r'ATST.*NEB',
            r'过渡态搜索.*ATST',
            r'ATST.*过渡态搜索',
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """从教程提取 NEB 测试信息（仅 Case 1 Li-Si IS SCF）"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="Li-in-Si-IS"
        )

        # ── INPUT：使用内置模板（含教程中 Li-Si 案例的 ABACUS 参数）──
        test_info.input_content = _LI_SI_INPUT
        print("    [INPUT] Li-Si IS SCF 模板（basis_type=lcao, kpts=2×2×2, nspin=1）")

        # ── STRU：Si 金刚石常规晶胞 + Li 八面体间隙（近似初始态）──
        test_info.stru_content = _LI_SI_STRU
        print("    [STRU] Si 金刚石 8 原子胞 + Li 间隙（近似 NEB 初始态）")

        # ── KPT ──
        test_info.kpt_content = _LI_SI_KPT
        print("    [KPT] 2×2×2 Gamma 中心 k 网格")

        # ── 预期结果（定性验证）──
        test_info.expected_results = {
            'scf_converged': 1.0,         # SCF 必须收敛
            'energy_negative': -1.0,       # 总能量应为负（符号标记 -1）
        }

        # ── 赝势 + 轨道文件（来自教程 Case 1）──
        test_info.pseudopotentials = [
            'Li_ONCV_PBE-1.2.upf',
            'Si_ONCV_PBE-1.2.upf',
        ]
        test_info.orbitals = [
            'Li_gga_8au_100Ry_4s1p.orb',
            'Si_gga_8au_100Ry_2s2p1d.orb',
        ]

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备 Li-Si IS SCF 输入文件"""
        input_dir = work_dir / f"{test_info.case_name}_neb"
        input_dir.mkdir(exist_ok=True)

        # 写入 ABACUS 输入文件
        (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        print(f"  [写入] INPUT / STRU / KPT → {input_dir}")

        # 下载赝势文件
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                try:
                    pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                    shutil.copy(pp_file, input_dir / pp)
                    print(f"  [PP]  {pp}")
                except Exception as e:
                    print(f"  [WARN] 赝势下载失败 {pp}: {e}")

            # 下载轨道文件
            for orb in test_info.orbitals:
                try:
                    orb_file = self.pp_manager.get_file(orb, "orbital")
                    shutil.copy(orb_file, input_dir / orb)
                    print(f"  [ORB] {orb}")
                except Exception as e:
                    print(f"  [WARN] 轨道文件下载失败 {orb}: {e}")

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 Li-Si IS SCF 任务"""
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
        """验证 Li-Si IS SCF 结果（SCF 收敛 + 能量为负）"""
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

        # ── 验证总能量为负（物理合理性检查）──
        total_e = actual.get('total_energy', None)
        if total_e is not None:
            passed_e = total_e < 0
            validation.comparisons['energy_negative'] = {
                'key': 'energy_negative',
                'expected': -1.0,  # 符号标记
                'actual': total_e,
                'abs_error': 0.0 if passed_e else abs(total_e),
                'rel_error': 0.0 if passed_e else 1.0,
                'passed': passed_e,
                'use_abs': True,
                'description': f'总能量应 < 0（实际: {total_e:.4f} eV）'
            }
        else:
            validation.warnings.append("未找到总能量输出（FINAL_ETOT_IS）")

        # ── 汇总 ──
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的结果（日志文件可能不存在）")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成测试报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"
        report += "> 验证方式：Li in Si 初始态单点 SCF（验证赝势/轨道文件/LCAO参数可用性）\n"
        report += "> 注：完整 NEB 迭代成本过高，仅验证参数层面的正确性。\n\n"

        if validation.comparisons:
            report += "| 验证项 | 期望 | 实际值 | 状态 |\n"
            report += "|--------|------|--------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                desc = comp.get('description', key)
                if isinstance(comp['actual'], float):
                    actual_str = f"{comp['actual']:.4f}"
                else:
                    actual_str = str(comp['actual'])
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

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从 ABACUS 输出文件提取 SCF 结果"""
        results = {}

        log_files = list(results_dir.rglob("running_scf.log"))
        if not log_files:
            log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            text = log_file.read_text(encoding='utf-8', errors='ignore')

            # SCF 收敛
            if 'convergence is achieved' in text or 'SCF CONVERGED' in text.upper():
                results['scf_converged'] = 1.0
            else:
                results['scf_converged'] = 0.0

            # 总能量
            m = re.search(r'FINAL_ETOT_IS\s+([-\d.]+)\s+eV', text)
            if m:
                results['total_energy'] = float(m.group(1))

            if results:
                break

        return results
