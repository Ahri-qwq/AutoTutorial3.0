"""
AutoTutorial 3.0 - SDFT/MDFT 随机密度泛函理论测试插件
处理含 esolver_type = sdft 的平面波计算（SCF、MD、DOS）

案例：
  case1 - pw_Si2       : Si 金刚石，MDFT SCF，T=0.6 Ry（≈8.16 eV）
  case2 - pw_md_Al     : Al FCC 16原子，纯SDFT MD，T=7.35 Ry（≈100 eV）
  case3 - pw_sdos_Si   : Si 1原子，MDFT+DOS SCF，T=0.6 Ry

验证策略：
  教程未提供具体数值，故验证"计算是否正常完成"：
  - SCF/MD 是否收敛/完成
  - Chebyshev Precision < 1e-8
  - 总能量为有限负数
  - DOS 案例：DOS1_smearing.dat 是否生成
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────────────────────────────
# 硬编码的输入文件内容
# STRU 来源：标准晶体学数据，Si 金刚石 a=10.2 Bohr，Al FCC a=7.65 Bohr
# ─────────────────────────────────────────────────────────────────────

# ── case1: pw_Si2 ── Si 金刚石，2原子，MDFT SCF ──
_SI2_INPUT = """\
INPUT_PARAMETERS
#Parameters (General)
calculation    scf
esolver_type   sdft
pseudo_dir     ./
nbands         4
nbands_sto     64
nche_sto       100
method_sto     1
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
symmetry       1
#Parameters (Smearing)
smearing_method fd
smearing_sigma  0.6
"""

_SI2_STRU = """\
ATOMIC_SPECIES
Si  28.085  Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
1.0  0.0  0.0
0.0  1.0  0.0
0.0  0.0  1.0

ATOMIC_POSITIONS
Direct

Si
0
2
0.00  0.00  0.00  1  1  1
0.25  0.25  0.25  1  1  1
"""

_SI2_KPT = """\
K_POINTS
0
Gamma
2 2 2 0 0 0
"""

# ── case2: pw_md_Al ── Al FCC 16原子，纯SDFT MD ──
# 2×2×1 超胞（conventional FCC a=7.65 Bohr × 2=15.3），16个原子
_AL_INPUT = """\
INPUT_PARAMETERS
#Parameters (General)
calculation    md
esolver_type   sdft
pseudo_dir     ./
nbands         0
nbands_sto     64
nche_sto       20
method_sto     2
#Parameters (Accuracy)
ecutwfc        50
scf_nmax       20
scf_thr        1e-6
symmetry       0
#Parameters (Smearing)
smearing_method fd
smearing_sigma  7.34986072
#Parameters (MD)
md_tfirst      1160400
md_dt          0.2
md_nstep       3
"""

# 16 Al 原子 FCC 结构，2×2×1 超胞（conventional FCC a=7.65 Bohr）
# 原FCC 4原子 → 4×4=16个原子
_AL_STRU = """\
ATOMIC_SPECIES
Al  26.982  Al_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
15.30  0.00   0.00
0.00   15.30  0.00
0.00   0.00   7.65

ATOMIC_POSITIONS
Direct

Al
0
16
0.00  0.00  0.00  1  1  1
0.50  0.00  0.00  1  1  1
0.00  0.50  0.00  1  1  1
0.50  0.50  0.00  1  1  1
0.25  0.25  0.00  1  1  1
0.75  0.25  0.00  1  1  1
0.25  0.75  0.00  1  1  1
0.75  0.75  0.00  1  1  1
0.25  0.00  0.50  1  1  1
0.75  0.00  0.50  1  1  1
0.25  0.50  0.50  1  1  1
0.75  0.50  0.50  1  1  1
0.00  0.25  0.50  1  1  1
0.50  0.25  0.50  1  1  1
0.00  0.75  0.50  1  1  1
0.50  0.75  0.50  1  1  1
"""

_AL_KPT = """\
K_POINTS
0
Gamma
1 1 1 0 0 0
"""

# ── case3: pw_sdos_Si ── Si 1原子，MDFT+DOS SCF ──
# FCC 原胞（1 Si），晶格矢量 a=(0,5.1,5.1),(5.1,0,5.1),(5.1,5.1,0) Bohr
_SDOS_INPUT = """\
INPUT_PARAMETERS
#Parameters (1.General)
suffix          autotest
calculation     scf
esolver_type    sdft
method_sto      2
nbands          10
nbands_sto      10
nche_sto        120
emax_sto        0
emin_sto        0
seed_sto        20000
pseudo_dir      ./
symmetry        1
kpar            1
bndpar          1
#Parameters (2.Iteration)
ecutwfc         20
scf_thr         1e-6
scf_nmax        20
#Parameters (3.Basis)
basis_type      pw
#Parameters (4.Smearing)
smearing_method fd
smearing_sigma  0.6
#Parameters (5.Mixing)
mixing_type     broyden
mixing_beta     0.4
out_dos         1
dos_emin_ev     -20
dos_emax_ev     100
dos_edelta_ev   0.1
dos_sigma       4
dos_nche        240
npart_sto       1
"""

# Si FCC 原胞（1 原子），a_conventional=10.2 Bohr
# FCC 原胞矢量：(0,a/2,a/2),(a/2,0,a/2),(a/2,a/2,0) = (0,5.1,5.1) 等
_SDOS_STRU = """\
ATOMIC_SPECIES
Si  28.085  Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
0.0  5.1  5.1
5.1  0.0  5.1
5.1  5.1  0.0

ATOMIC_POSITIONS
Direct

Si
0
1
0.00  0.00  0.00  1  1  1
"""

_SDOS_KPT = """\
K_POINTS
0
Gamma
4 4 4 0 0 0
"""


class SDFTPlugin(BaseTestPlugin):
    """SDFT/MDFT 随机密度泛函理论测试插件"""

    @property
    def plugin_name(self) -> str:
        return "SDFT/MDFT 随机密度泛函理论测试"

    @property
    def calc_type(self) -> str:
        return "sdft"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测 esolver_type = sdft（在代码块中）"""
        # 在 INPUT 代码块中搜索 esolver_type = sdft
        input_block_pattern = r'```[\s\S]*?INPUT_PARAMETERS[\s\S]*?```'
        input_blocks = re.findall(input_block_pattern, tutorial_content)
        for block in input_blocks:
            if re.search(r'esolver_type\s*=?\s*sdft', block, re.IGNORECASE):
                return True
        # 也匹配正文中的明确标识
        return bool(re.search(r'esolver_type\s*=?\s*sdft', tutorial_content, re.IGNORECASE))

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取 SDFT 测试信息 —— 固定返回三个案例中的第一个（pw_Si2 SCF）"""
        # 使用 case_name="sdft_si2" 标识 pw_Si2 SCF 案例
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="sdft_si2"
        )
        test_info.input_content = _SI2_INPUT
        test_info.stru_content = _SI2_STRU
        test_info.kpt_content = _SI2_KPT
        test_info.pseudopotentials = ["Si_ONCV_PBE-1.0.upf"]
        test_info.orbitals = []
        # 无数值预期，验证收敛性
        test_info.expected_results = {}
        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备三个案例的输入文件"""
        input_dirs = []

        cases = [
            ("sdft_si2_scf",   _SI2_INPUT,   _SI2_STRU,   _SI2_KPT,   ["Si_ONCV_PBE-1.0.upf"]),
            ("sdft_al_md",     _AL_INPUT,    _AL_STRU,    _AL_KPT,    ["Al_ONCV_PBE-1.0.upf"]),
            ("sdft_si_sdos",   _SDOS_INPUT,  _SDOS_STRU,  _SDOS_KPT,  ["Si_ONCV_PBE-1.0.upf"]),
        ]

        for case_name, inp, stru, kpt, pps in cases:
            case_dir = work_dir / case_name
            case_dir.mkdir(exist_ok=True)

            (case_dir / "INPUT").write_text(inp, encoding="utf-8")
            (case_dir / "STRU").write_text(stru, encoding="utf-8")
            (case_dir / "KPT").write_text(kpt, encoding="utf-8")

            # 下载赝势
            if self.pp_manager:
                for pp in pps:
                    try:
                        pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                        shutil.copy(pp_file, case_dir / pp)
                    except Exception as e:
                        print(f"  [WARN] 赝势下载失败 {pp}: {e}")

            input_dirs.append(case_dir)

        return input_dirs

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交三个 SDFT 计算任务"""
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
        """验证三个案例的计算结果"""
        validation = ValidationResult(
            calc_type=self.calc_type,
            case_name="SDFT三案例",
            passed=False,
            comparisons={}
        )

        results_dir = work_dir / "results"
        results_dir.mkdir(exist_ok=True)

        for job_id in job_ids:
            self.job_manager.download_job_result(job_id, results_dir)

        # ── 验证 pw_Si2 SCF ──
        si2_logs = list(results_dir.rglob("sdft_si2_scf/OUT.*/running_scf.log"))
        if si2_logs:
            si2_result = self._validate_scf_log(si2_logs[0], "sdft_si2_scf")
            validation.comparisons.update(si2_result)
        else:
            validation.warnings.append("未找到 sdft_si2_scf 的 running_scf.log")

        # ── 验证 pw_md_Al MD ──
        al_logs = list(results_dir.rglob("sdft_al_md/OUT.*/running_md.log"))
        if al_logs:
            al_result = self._validate_md_log(al_logs[0], "sdft_al_md")
            validation.comparisons.update(al_result)
        else:
            validation.warnings.append("未找到 sdft_al_md 的 running_md.log")

        # ── 验证 DOS 案例 ──
        dos_files = list(results_dir.rglob("sdft_si_sdos/OUT.autotest/DOS1_smearing.dat"))
        dos_ok = len(dos_files) > 0
        validation.comparisons["sdft_si_sdos_dos_file"] = {
            "key": "DOS1_smearing.dat 是否生成",
            "expected": "存在",
            "actual": "存在" if dos_ok else "不存在",
            "abs_error": 0 if dos_ok else 1,
            "rel_error": 0 if dos_ok else 1,
            "passed": dos_ok
        }
        if not dos_ok:
            validation.warnings.append("未找到 DOS1_smearing.dat，DOS 案例可能未完成")

        # 综合判断
        if validation.comparisons:
            validation.passed = all(c.get("passed", False) for c in validation.comparisons.values())
        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成 Markdown 报告章节"""
        report = f"### {self.plugin_name}\n\n"
        report += "| 案例 | 验证项 | 预期 | 实际 | 状态 |\n"
        report += "|------|--------|------|------|------|\n"

        for key, comp in validation.comparisons.items():
            status = "✅ PASS" if comp.get("passed") else "❌ FAIL"
            report += f"| {key} | {comp.get('key','-')} | {comp.get('expected','-')} | {comp.get('actual','-')} | {status} |\n"

        if validation.warnings:
            report += "\n**警告：**\n"
            for w in validation.warnings:
                report += f"- {w}\n"

        overall = "✅ 通过" if validation.passed else "❌ 失败"
        report += f"\n**总体结果：** {overall}\n\n"
        return report

    # ── 辅助验证方法 ──

    def _validate_scf_log(self, log_path: Path, case_name: str) -> Dict:
        """从 running_scf.log 中验证 SCF 收敛和 Chebyshev Precision"""
        results = {}
        content = log_path.read_text(encoding="utf-8", errors="ignore")

        # 验证 SCF 收敛（ABACUS SDFT 收敛标志）
        scf_ok = ('charge density convergence is achieved' in content or
                  bool(re.search(r'SCF\s+CONVERGE', content, re.IGNORECASE)))
        results[f"{case_name}_scf_converged"] = {
            "key": "SCF 收敛",
            "expected": "是",
            "actual": "是" if scf_ok else "否",
            "abs_error": 0 if scf_ok else 1,
            "rel_error": 0 if scf_ok else 1,
            "passed": scf_ok
        }

        # 验证 Chebyshev Precision（取最后一步的值）
        cheb_matches = re.findall(r'Chebyshev\s+Precision\s*[=:]\s*([0-9eE+\-.]+)', content)
        if cheb_matches:
            last_cheb = float(cheb_matches[-1])
            cheb_ok = last_cheb < 1e-8
            results[f"{case_name}_chebyshev"] = {
                "key": "Chebyshev Precision < 1e-8",
                "expected": "< 1e-8",
                "actual": f"{last_cheb:.2e}",
                "abs_error": 0 if cheb_ok else 1,
                "rel_error": 0 if cheb_ok else 1,
                "passed": cheb_ok
            }
        else:
            results[f"{case_name}_chebyshev"] = {
                "key": "Chebyshev Precision < 1e-8",
                "expected": "< 1e-8",
                "actual": "未找到",
                "abs_error": 1,
                "rel_error": 1,
                "passed": False
            }

        # 提取总能量（两列格式：E_KohnSham  <per-unit>  <total-eV>）
        etot_matches = re.findall(r'E_KohnSham\s+[+-]?[\d\.]+(?:[eE][+-]?\d+)?\s+([+-]?[\d\.]+(?:[eE][+-]?\d+)?)', content)
        if etot_matches:
            etot = float(etot_matches[-1])
            etot_ok = etot < 0
            results[f"{case_name}_etot"] = {
                "key": "总能量为负数",
                "expected": "< 0 eV",
                "actual": f"{etot:.4f} eV",
                "abs_error": 0 if etot_ok else 1,
                "rel_error": 0 if etot_ok else 1,
                "passed": etot_ok
            }

        return results

    def _validate_md_log(self, log_path: Path, case_name: str) -> Dict:
        """验证 MD 日志：步数是否完成"""
        results = {}
        content = log_path.read_text(encoding="utf-8", errors="ignore")

        # 查找 MD 步数（格式：STEP OF MOLECULAR DYNAMICS : N）
        step_matches = re.findall(r'STEP OF MOLECULAR DYNAMICS\s*:\s*(\d+)', content)
        completed_steps = max((int(s) for s in step_matches), default=0)
        md_ok = completed_steps >= 3  # INPUT 设置 md_nstep=3

        results[f"{case_name}_md_steps"] = {
            "key": "MD 步数完成（≥3步）",
            "expected": "≥ 3",
            "actual": str(completed_steps),
            "abs_error": 0 if md_ok else 1,
            "rel_error": 0 if md_ok else 1,
            "passed": md_ok
        }

        return results
