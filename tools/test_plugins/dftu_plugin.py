"""
AutoTutorial 3.0 - DFT+U 强关联体系测试插件
处理含 dft_plus_u 的 LCAO-SCF 计算（如反铁磁 NiO）
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────
# 硬编码的 NiO 反铁磁算例文件内容
# ─────────────────────────────────────────────

# 晶格矢量来源：ABACUS-Agent-tools RAG 数据（已验证 Ni-O 键长 ≈ 2.09 Å）
# Type-II AFM NiO：Ni1(↑) 和 Ni2(↓) 位于体对角线两端
_NIO_STRU_TEMPLATE = """\
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf
Ni2 58.693 Ni_ONCV_PBE-1.0.upf
O   15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Ni_gga_9au_100Ry_4s2p2d1f.orb
Ni_gga_9au_100Ry_4s2p2d1f.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
9.6226  0.0000  0.0000    #latteice vector A1 (Bohr)
7.9999  5.3468  0.0000    #latteice vector A2 (Bohr)
7.9999  2.4270  4.7629    #latteice vector A3 (Bohr)

ATOMIC_POSITIONS
Direct

Ni1
2.0
1
0.000000000  0.000000000  0.000000000  0  0  0

Ni2
-2.0
1
0.500000000  0.500000000  0.500000000  0  0  0

O
0.0
2
0.250000000  0.250000000  0.250000000  0  0  0
0.750000000  0.750000000  0.750000000  0  0  0
"""

_NIO_INPUT_TEMPLATE = """\
INPUT_PARAMETERS
#Parameters (General)
suffix                  NiO
calculation             scf
esolver_type            ksdft
symmetry                0

#Parameters (Basis)
basis_type              lcao
ecutwfc                 100

#Parameters (Spin)
nspin                   2

#Parameters (Iteration)
scf_thr                 1e-7
scf_nmax                200

#Parameters (Smearing)
smearing_method         gauss
smearing_sigma          0.01

#Parameters (Mixing)
mixing_type             broyden
mixing_beta             0.4

#Parameters (Output)
out_bandgap             1
out_mul                 1
out_chg                 1

#Parameters (DFT+U)
dft_plus_u              1
orbital_corr            2 2 -1
hubbard_u               5.0 5.0 0.0
"""

_NIO_KPT = """\
K_POINTS
0
Gamma
4 4 4 0 0 0
"""


class DFTUPlugin(BaseTestPlugin):
    """DFT+U 强关联体系测试插件"""

    def __init__(self, job_manager=None, pp_manager=None):
        super().__init__(job_manager, pp_manager)
        self.tolerance = 0.001   # 0.1% 相对容差（用于能量等绝对量）
        self.abs_tol_energy = 1.0   # 允许能量绝对偏差 1 eV
        self.abs_tol_mag = 0.1      # 允许磁矩绝对偏差 0.1

    @property
    def plugin_name(self) -> str:
        return "DFT+U 强关联体系测试"

    @property
    def calc_type(self) -> str:
        return "dftu"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测是否包含 DFT+U 相关参数"""
        patterns = [
            r'dft_plus_u\s*[=\s]+1',
            r'DFT\+U',
            r'dft\+u',
            r'Hubbard.*U',
            r'强关联.*DFT',
        ]
        return any(re.search(p, tutorial_content, re.IGNORECASE) for p in patterns)

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """从教程提取 DFT+U 测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="NiO"
        )

        # ── 尝试从教程提取 INPUT，失败则使用模板 ──
        extracted_input = self._extract_input(tutorial_content)
        if extracted_input and 'dft_plus_u' in extracted_input:
            test_info.input_content = extracted_input
            print("    [INPUT] 从教程提取成功")
        else:
            test_info.input_content = _NIO_INPUT_TEMPLATE
            print("    [INPUT] 使用内置 NiO 模板")

        # ── STRU：使用完整的 NiO 反铁磁结构（教程中仅给出示意坐标）──
        test_info.stru_content = _NIO_STRU_TEMPLATE
        print("    [STRU] 使用 NiO 反铁磁结构模板（Ni-O 键长 2.09 Å）")

        # ── KPT：4×4×4 Gamma ──
        test_info.kpt_content = _NIO_KPT
        print("    [KPT] 4×4×4 Gamma 均匀 k 网格")

        # ── 预期结果（来自 Bohrium 算例 abacus-magnetic-eu2y/v4）──
        test_info.expected_results = {
            'total_energy':       -9255.7279034240546025,  # eV
            'bandgap_up':         0.205369322748,           # eV (spin-up)
            'bandgap_dn':         2.794192983776,           # eV (spin-down)
            'total_magnetism':    0.0,                      # Bohr mag/cell
            'absolute_magnetism': 3.35321634,               # Bohr mag/cell
        }

        # ── 赝势与轨道文件 ──
        test_info.pseudopotentials = [
            'Ni_ONCV_PBE-1.0.upf',
            'O_ONCV_PBE-1.0.upf',
        ]
        test_info.orbitals = [
            'Ni_gga_9au_100Ry_4s2p2d1f.orb',
            'O_gga_7au_100Ry_2s2p1d.orb',
        ]

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备 NiO DFT+U SCF 输入文件"""
        input_dir = work_dir / f"{test_info.case_name}_dftu"
        input_dir.mkdir(exist_ok=True)

        # 写入输入文件
        (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        print(f"  [写入] INPUT / STRU / KPT → {input_dir}")

        # 下载赝势和轨道文件
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                try:
                    pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                    shutil.copy(pp_file, input_dir / pp)
                    print(f"  [PP]  {pp}")
                except Exception as e:
                    print(f"  [WARN] 赝势下载失败 {pp}: {e}")

            for orb in test_info.orbitals:
                try:
                    orb_file = self.pp_manager.get_file(orb, "orbital")
                    shutil.copy(orb_file, input_dir / orb)
                    print(f"  [ORB] {orb}")
                except Exception as e:
                    print(f"  [WARN] 轨道下载失败 {orb}: {e}")

        # 修复 INPUT 文件路径（pseudo_dir / orbital_dir）
        input_file = input_dir / "INPUT"
        self._fix_input_file(input_file)

        # 修复 STRU 文件路径格式（改为相对路径）
        stru_file = input_dir / "STRU"
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent.parent))
            from fix_stru import fix_stru_paths, fix_input_paths
            fix_stru_paths(stru_file, verbose=False)
            fix_input_paths(input_file, verbose=False)
        except Exception as e:
            print(f"  [WARN] 路径修复跳过: {e}")

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 DFT+U SCF 任务"""
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
        """验证 DFT+U SCF 结果"""
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

        # 对比各指标
        exp = test_info.expected_results

        # 总能量（相对容差 0.1%，比绝对值更合理——结构微小差异会导致 eV 级能量偏移）
        if 'total_energy' in exp and 'total_energy' in actual:
            comp = self._compare_values(exp['total_energy'], actual['total_energy'],
                                        'total_energy')
            # 使用 0.1% 相对容差（对应 9255 eV 体系约 9.3 eV，远宽于 SCF 精度）
            comp['passed'] = comp['rel_error'] <= 0.001
            validation.comparisons['total_energy'] = comp

        # spin-up 能隙（相对容差 15%）
        if 'bandgap_up' in exp and 'bandgap_up' in actual:
            comp = self._compare_values(exp['bandgap_up'], actual['bandgap_up'], 'bandgap_up')
            comp['passed'] = comp['rel_error'] <= 0.15
            validation.comparisons['bandgap_up'] = comp

        # spin-down 能隙（相对容差 15%）
        if 'bandgap_dn' in exp and 'bandgap_dn' in actual:
            comp = self._compare_values(exp['bandgap_dn'], actual['bandgap_dn'], 'bandgap_dn')
            comp['passed'] = comp['rel_error'] <= 0.15
            validation.comparisons['bandgap_dn'] = comp

        # 绝对磁矩（相对容差 5%）
        if 'absolute_magnetism' in exp and 'absolute_magnetism' in actual:
            comp = self._compare_values(exp['absolute_magnetism'],
                                        actual['absolute_magnetism'], 'absolute_magnetism')
            validation.comparisons['absolute_magnetism'] = comp

        # 总磁矩（应为 0，用绝对容差 0.1 Bohr mag）
        if 'total_magnetism' in exp and 'total_magnetism' in actual:
            comp = self._compare_abs(exp['total_magnetism'], actual['total_magnetism'],
                                     'total_magnetism', self.abs_tol_mag)
            validation.comparisons['total_magnetism'] = comp

        # 汇总
        if validation.comparisons:
            validation.passed = all(c['passed'] for c in validation.comparisons.values())
        else:
            validation.warnings.append("未找到可对比的结果（日志文件可能不存在）")

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成报告章节"""
        report = f"### {self.plugin_name} - {validation.case_name}\n\n"

        if validation.comparisons:
            report += "| 参数 | 预期值 | 实际值 | 误差 | 状态 |\n"
            report += "|------|--------|--------|------|------|\n"

            for key, comp in validation.comparisons.items():
                status = "✅ PASS" if comp['passed'] else "❌ FAIL"
                if comp.get('use_abs', False):
                    err_str = f"±{comp['abs_error']:.4f}（绝对）"
                else:
                    err_str = f"{comp['rel_error']*100:.2f}%（相对）"
                report += (f"| {key} | {comp['expected']:.6f} "
                           f"| {comp['actual']:.6f} | {err_str} | {status} |\n")

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

    def _fix_input_file(self, input_path: Path):
        """添加 pseudo_dir / orbital_dir（若缺失）"""
        text = input_path.read_text(encoding='utf-8')
        lines = text.split('\n')
        has_pseudo = any('pseudo_dir' in l.lower() for l in lines)
        has_orb = any('orbital_dir' in l.lower() for l in lines)
        if has_pseudo and has_orb:
            return

        new_lines = []
        added = False
        for line in lines:
            new_lines.append(line)
            if line.strip() == 'INPUT_PARAMETERS' and not added:
                if not has_pseudo:
                    new_lines.append('pseudo_dir      ./')
                if not has_orb:
                    new_lines.append('orbital_dir     ./')
                added = True
        input_path.write_text('\n'.join(new_lines), encoding='utf-8')

    def _compare_abs(self, expected: float, actual: float,
                     key: str, abs_tol: float) -> Dict:
        """使用绝对容差对比数值"""
        abs_error = abs(actual - expected)
        rel_error = abs_error / max(abs(expected), 1e-10)
        return {
            'key': key,
            'expected': expected,
            'actual': actual,
            'abs_error': abs_error,
            'rel_error': rel_error,
            'passed': abs_error <= abs_tol,
            'use_abs': True,
        }

    def _extract_actual_results(self, results_dir: Path) -> Dict:
        """从 ABACUS 输出文件提取计算结果"""
        results = {}

        log_files = list(results_dir.rglob("running_scf.log"))
        if not log_files:
            log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            text = log_file.read_text(encoding='utf-8', errors='ignore')

            # 总能量
            m = re.search(r'FINAL_ETOT_IS\s+([-\d.]+)\s+eV', text)
            if m:
                results['total_energy'] = float(m.group(1))

            # 能隙（最后一次出现的 E_bandgap 行）
            bg_matches = re.findall(
                r'E_bandgap\s+([+\-\d.]+)\s+([+\-\d.]+)', text)
            if bg_matches:
                up, dn = bg_matches[-1]
                results['bandgap_up'] = float(up)
                results['bandgap_dn'] = float(dn)

            # 磁矩（最后一次出现）
            tm_matches = re.findall(
                r'total magnetism \(Bohr mag/cell\) =\s*([-\d.]+)', text)
            if tm_matches:
                results['total_magnetism'] = float(tm_matches[-1])

            am_matches = re.findall(
                r'absolute magnetism \(Bohr mag/cell\) =\s*([\d.]+)', text)
            if am_matches:
                results['absolute_magnetism'] = float(am_matches[-1])

            if results:
                break   # 找到一个有效日志即可

        return results
