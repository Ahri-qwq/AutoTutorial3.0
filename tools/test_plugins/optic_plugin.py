"""
AutoTutorial 3.0 - ABACUS + PYATB 光学性质测试插件
处理 ABACUS LCAO-SCF + PYATB OPTICAL_CONDUCTIVITY 两步计算流程
首次测试教程：光学性质（介电函数/吸收谱）计算教程（晶态 SiO₂）
"""

import re
import json
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────
# 硬编码的 β-cristobalite SiO₂ 案例文件内容
# 参数来源：教程案例原文（ABACUS SCF 已验证运行，NBASE=312）
# ─────────────────────────────────────────────

_SIO2_INPUT = """\
INPUT_PARAMETERS

suffix                  silica
calculation             scf
esolver_type            ksdft
symmetry                0
init_chg                atomic

pseudo_dir              ./
orbital_dir             ./

basis_type              lcao
ks_solver               genelpa
smearing_method         gaussian
smearing_sigma          0.01
mixing_type             broyden
mixing_beta             0.1
ecutwfc                 100
scf_thr                 1e-8
mixing_gg0              1.5
mixing_ndim             20

out_chg                 1
out_mat_hs2             1
out_mat_r               1
"""

_SIO2_STRU = """\
ATOMIC_SPECIES
Si 28.086  Si_ONCV_PBE-1.0.upf
O 15.999   O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_7au_100Ry_2s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261246257702

LATTICE_VECTORS
7.1199998856 0.0 0.0
0.0 7.1199998856 0.0
0.0 0.0 7.1199998856

ATOMIC_POSITIONS
Cartesian
Si
0.0
8
0.000000000000 0.000000000000 0.000000000000 1 1 1
0.000000000000 3.559999943000 3.559999943000 1 1 1
3.559999943000 3.559999943000 0.000000000000 1 1 1
3.559999943000 0.000000000000 3.559999943000 1 1 1
5.339999914000 1.779999971000 5.339999914000 1 1 1
1.779999971000 1.779999971000 1.779999971000 1 1 1
1.779999971000 5.339999914000 5.339999914000 1 1 1
5.339999914000 5.339999914000 1.779999971000 1 1 1
O
0.0
16
0.889999986000 0.889999986000 0.889999986000 1 1 1
6.229999900000 2.669999957000 4.449999928000 1 1 1
2.669999957000 4.449999928000 6.229999900000 1 1 1
4.449999928000 6.229999900000 2.669999957000 1 1 1
0.889999986000 4.449999928000 4.449999928000 1 1 1
6.229999900000 6.229999900000 0.889999986000 1 1 1
2.669999957000 0.889999986000 2.669999957000 1 1 1
4.449999928000 2.669999957000 6.229999900000 1 1 1
4.449999928000 0.889999986000 4.449999928000 1 1 1
2.669999957000 2.669999957000 0.889999986000 1 1 1
6.229999900000 4.449999928000 2.669999957000 1 1 1
0.889999986000 6.229999900000 6.229999900000 1 1 1
4.449999928000 4.449999928000 0.889999986000 1 1 1
2.669999957000 6.229999900000 4.449999928000 1 1 1
6.229999900000 0.889999986000 6.229999900000 1 1 1
0.889999986000 2.669999957000 2.669999957000 1 1 1
"""

_SIO2_KPT = """\
K_POINTS
0
Gamma
6 6 6 0 0 0
"""

_SIO2_PYATB_INPUT = """\
INPUT_PARAMETERS
{
nspin               1
package             ABACUS
fermi_energy        5.5385382545
fermi_energy_unit   eV
HR_route            data-HR-sparse_SPIN0.csr
SR_route            data-SR-sparse_SPIN0.csr
rR_route            data-rR-sparse.csr
HR_unit             Ry
rR_unit             Bohr
}

LATTICE
{
lattice_constant        1.8897261246257702
lattice_constant_unit   Bohr
lattice_vector
7.1199998856 0.0 0.0
0.0 7.1199998856 0.0
0.0 0.0 7.1199998856
}

OPTICAL_CONDUCTIVITY
{
occ_band    64
omega       0 30
domega      0.01
eta         0.1
grid        20 20 20
}
"""

# 单任务完成 ABACUS SCF → 复制矩阵文件 → PYATB 光学计算
# 使用镜像：abacus-pyatb-open（内置 ABACUS + PYATB）
_JOB_COMMAND = (
    "export ABACUS_PP_PATH=./ && export ABACUS_ORB_PATH=./ && "
    "OMP_NUM_THREADS=1 mpirun -np 16 abacus && "
    "mkdir -p pyatb && "
    "cp OUT.silica/data-HR-sparse_SPIN0.csr "
    "OUT.silica/data-SR-sparse_SPIN0.csr "
    "OUT.silica/data-rR-sparse.csr pyatb/ && "
    "cd pyatb && "
    "OMP_NUM_THREADS=1 mpirun -np 16 pyatb"
)

# 专用镜像（含 ABACUS + PYATB，来自教程计算环境说明）
_PYATB_IMAGE = "registry.dp.tech/dptech/prod-19853/abacus-pyatb-open:v0.0.1"


class OpticPlugin(BaseTestPlugin):
    """ABACUS + PYATB 光学性质测试插件（介电函数/吸收谱）"""

    def __init__(self, job_manager=None, pp_manager=None):
        super().__init__(job_manager, pp_manager)
        self.tolerance = 0.05   # 默认相对容差 5%

    @property
    def plugin_name(self) -> str:
        return "ABACUS+PYATB 光学性质测试"

    @property
    def calc_type(self) -> str:
        return "optic"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测是否为 ABACUS + PYATB 光学计算教程"""
        patterns = [
            r'out_mat_hs2',           # 输出哈密顿量矩阵（PYATB 必需）
            r'out_mat_r',             # 输出偶极矩阵（PYATB 必需）
            r'OPTICAL_CONDUCTIVITY',  # PYATB 光学电导率模块
            r'pyatb.*optical',        # PYATB 光学关键词
        ]
        # 必须同时满足：有 out_mat_hs2/out_mat_r（ABACUS侧）且有 OPTICAL_CONDUCTIVITY（PYATB侧）
        has_abacus_output = any(
            re.search(p, tutorial_content, re.IGNORECASE)
            for p in patterns[:2]
        )
        has_pyatb = any(
            re.search(p, tutorial_content, re.IGNORECASE)
            for p in patterns[2:]
        )
        return has_abacus_output and has_pyatb

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """从教程提取 SiO₂ 光学计算测试信息"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="SiO2_optic"
        )

        # INPUT：使用硬编码模板（已从教程案例提取，完全一致）
        test_info.input_content = _SIO2_INPUT
        print("    [INPUT] 使用 SiO₂ SCF 模板（basis_type=lcao, out_mat_hs2=1, out_mat_r=1）")

        # STRU：β-cristobalite SiO₂，8 Si + 16 O，晶格常数 7.12 Å
        test_info.stru_content = _SIO2_STRU
        print("    [STRU] 使用 β-方英石 SiO₂ 结构（24 个原子，立方超胞）")

        # KPT：6×6×6 Gamma 中心（与教程一致）
        test_info.kpt_content = _SIO2_KPT
        print("    [KPT] 6×6×6 Gamma 均匀 k 网格")

        # 预期结果（来自教程案例 SCF 输出，GE13 步收敛）
        test_info.expected_results = {
            'total_energy':    -7835.176,    # eV（来自 GE13: -7.835176e+03）
            'occupied_bands':  64.0,         # 8×4+16×6=128 电子，nspin=1 占 64 带
            'fermi_energy':    5.5385382545, # eV（来自 E_Fermi 第二列）
            'pyatb_completed': 1.0,          # 1=PYATB 输出文件存在，0=不存在
        }

        # 赝势与轨道文件（ONCV PBE，7au 截断——已在 Bohrium abacus-pyatb-open 镜像中验证）
        test_info.pseudopotentials = [
            'Si_ONCV_PBE-1.0.upf',
            'O_ONCV_PBE-1.0.upf',
        ]
        test_info.orbitals = [
            'Si_gga_7au_100Ry_2s2p1d.orb',
            'O_gga_7au_100Ry_2s2p1d.orb',
        ]

        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备 ABACUS SCF + PYATB 光学计算输入文件"""
        input_dir = work_dir / test_info.case_name
        input_dir.mkdir(exist_ok=True)

        # ── 写入 ABACUS 输入文件 ──
        (input_dir / "INPUT").write_text(test_info.input_content, encoding='utf-8')
        (input_dir / "STRU").write_text(test_info.stru_content, encoding='utf-8')
        (input_dir / "KPT").write_text(test_info.kpt_content, encoding='utf-8')
        print(f"  [写入] ABACUS INPUT / STRU / KPT → {input_dir}")

        # ── 创建 pyatb/ 子目录并写入 PYATB Input 文件 ──
        pyatb_dir = input_dir / "pyatb"
        pyatb_dir.mkdir(exist_ok=True)
        (pyatb_dir / "Input").write_text(_SIO2_PYATB_INPUT, encoding='utf-8')
        print(f"  [写入] PYATB Input → {pyatb_dir / 'Input'}")

        # ── 下载赝势文件 ──
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

        return [input_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 ABACUS+PYATB 两步计算任务（使用专用 abacus-pyatb-open 镜像）"""
        job_ids = []
        for input_dir in input_dirs:
            job_name = input_dir.name

            # 构造自定义 job config（覆盖默认镜像和命令）
            job_config = {
                "job_name": job_name,
                "command": _JOB_COMMAND,
                "log_file": "OUT.silica/running_scf.log",
                "backward_files": [],
                "project_id": self.job_manager.project_id,
                "platform": "ali",
                "job_type": "container",
                "machine_type": "c16_m32_cpu",
                "image_address": _PYATB_IMAGE,
            }

            job_id = self.job_manager.submit_job(job_config, input_dir)

            if job_id:
                job_ids.append(job_id)
                print(f"  [OK] 提交任务: {job_name} (Job ID: {job_id})")
                print(f"       镜像: {_PYATB_IMAGE}")
            else:
                print(f"  [ERROR] 提交失败: {job_name}")

        return job_ids

    def validate_results(self, job_ids: List[str], work_dir: Path,
                         test_info: TestInfo) -> ValidationResult:
        """验证 ABACUS SCF + PYATB 计算结果"""
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

        exp = test_info.expected_results

        # 1. 总能量（相对容差 0.1%，对应 -7835 eV 体系约 8 eV）
        if 'total_energy' in exp and 'total_energy' in actual:
            comp = self._compare_values(exp['total_energy'], actual['total_energy'],
                                        'total_energy')
            comp['passed'] = comp['rel_error'] <= 0.001
            validation.comparisons['total_energy'] = comp

        # 2. 占据能带数（精确匹配，允许绝对误差 < 0.5）
        if 'occupied_bands' in exp and 'occupied_bands' in actual:
            abs_err = abs(actual['occupied_bands'] - exp['occupied_bands'])
            validation.comparisons['occupied_bands'] = {
                'key': 'occupied_bands',
                'expected': exp['occupied_bands'],
                'actual': actual['occupied_bands'],
                'abs_error': abs_err,
                'rel_error': abs_err / max(abs(exp['occupied_bands']), 1),
                'passed': abs_err < 0.5,
            }

        # 3. 费米能（相对容差 5%）
        if 'fermi_energy' in exp and 'fermi_energy' in actual:
            comp = self._compare_values(exp['fermi_energy'], actual['fermi_energy'],
                                        'fermi_energy')
            comp['passed'] = comp['rel_error'] <= 0.05
            validation.comparisons['fermi_energy'] = comp

        # 4. PYATB 输出文件是否存在（布尔判断）
        if 'pyatb_completed' in exp:
            pyatb_output = self._check_pyatb_output(results_dir)
            validation.comparisons['pyatb_completed'] = {
                'key': 'pyatb_completed',
                'expected': 1.0,
                'actual': 1.0 if pyatb_output else 0.0,
                'abs_error': 0.0 if pyatb_output else 1.0,
                'rel_error': 0.0 if pyatb_output else 1.0,
                'passed': pyatb_output,
            }
            if not pyatb_output:
                validation.warnings.append(
                    "PYATB 输出文件不存在（dielectric_function_*.dat），"
                    "请检查 PYATB 是否正常运行"
                )

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
                elif key == 'occupied_bands':
                    err_str = f"±{comp['abs_error']:.1f}（绝对）"
                elif key == 'pyatb_completed':
                    err_str = "布尔" if comp['passed'] else "未完成"
                else:
                    err_str = f"{comp['rel_error']*100:.3f}%（相对）"
                report += (f"| {key} | {comp['expected']:.4f} "
                           f"| {comp['actual']:.4f} | {err_str} | {status} |\n")

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
        """从 ABACUS 和 PYATB 输出文件提取计算结果"""
        results = {}

        # 查找 running_scf.log
        log_files = list(results_dir.rglob("running_scf.log"))
        if not log_files:
            log_files = list(results_dir.rglob("running_*.log"))

        for log_file in log_files:
            text = log_file.read_text(encoding='utf-8', errors='ignore')

            # 总能量（最后一步 SCF）
            m = re.search(r'FINAL_ETOT_IS\s+([-\d.eE+]+)\s+eV', text)
            if m:
                results['total_energy'] = float(m.group(1))

            # 占据能带数
            m = re.search(r'occupied bands\s*=\s*(\d+)', text)
            if m:
                results['occupied_bands'] = float(m.group(1))

            # 费米能（最后一次出现的 E_Fermi，取第二列 eV 值）
            fermi_matches = re.findall(
                r'E_Fermi\s+([-\d.eE+]+)\s+([-\d.eE+]+)', text)
            if fermi_matches:
                results['fermi_energy'] = float(fermi_matches[-1][1])

            if results:
                break

        return results

    def _check_pyatb_output(self, results_dir: Path) -> bool:
        """检查 PYATB 介电函数输出文件是否存在"""
        # PYATB 输出路径（在 pyatb/ 目录下运行后生成）
        patterns = [
            "**/dielectric_function_imag_part.dat",
            "**/dielectric_function_real_part.dat",
        ]
        for pattern in patterns:
            if list(results_dir.rglob(pattern.replace("**/", ""))):
                return True
        # 也检查不带子目录的情况
        for pattern in patterns:
            if list(results_dir.glob(pattern)):
                return True
        return False
