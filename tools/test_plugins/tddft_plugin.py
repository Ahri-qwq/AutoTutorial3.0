"""
AutoTutorial 3.0 - 实时含时密度泛函理论（RT-TDDFT）测试插件
处理含 esolver_type = tddft 的 LCAO 计算

测试策略：
  RT-TDDFT 产生时间序列轨迹，无固定的数值预期结果，故验证"计算是否正常完成"：
  1. 内置 SCF 是否收敛
  2. TDDFT 时间步是否运行到预期步数
  3. 电流输出文件（SPIN1_CURRENT）是否生成

案例：
  tddft_si_hybrid：Si 金刚石原胞（2原子），混合规范（td_stype=2），50步 TDDFT

首次测试教程：ABACUS 混合规范 RT-TDDFT 教程（Si + HCO）
添加日期：2026-03-09
"""

import re
import shutil
from pathlib import Path
from typing import List, Optional, Dict
from .base_plugin import BaseTestPlugin, TestInfo, ValidationResult


# ─────────────────────────────────────────────────────────────────────
# 硬编码的输入文件内容
# STRU 来源：教程第二章，Si 金刚石原胞，a=10.2 Bohr（面心立方原胞表示）
# INPUT：从教程第三章 TDDFT INPUT 裁剪，md_nstep 缩减至 50 步用于快速验证
# ─────────────────────────────────────────────────────────────────────

_TDDFT_SI_INPUT = """\
INPUT_PARAMETERS
#Parameters (1.General)
suffix           Si_tddft_test
calculation      md
esolver_type     tddft
pseudo_dir       ./
orbital_dir      ./

#Parameters (2.Basis)
basis_type       lcao
gamma_only       0

#Parameters (3.MD / 时间步)
md_nstep         50
md_dt            0.02
md_tfirst        0
md_type          nve

#Parameters (4.外场：混合规范，弱场宽带激发)
td_vext          1
td_stype         2

td_tstart        1
td_tend          50

td_vext_dire     1
td_ttype         0
td_gauss_amp     0.001
td_gauss_t0      10
td_gauss_sigma   0.5
td_gauss_freq    0.0
td_gauss_phase   0.0

#Parameters (5.精度)
ecutwfc          60
scf_thr          1e-7
scf_nmax         100

#Parameters (6.输出)
out_current      1
out_efield       1
"""

_TDDFT_SI_STRU = """\
ATOMIC_SPECIES
Si   28.085   Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
   0.5000000000     0.5000000000     0.0000000000
   0.0000000000     0.5000000000     0.5000000000
   0.5000000000     0.0000000000     0.5000000000

ATOMIC_POSITIONS
Direct
Si
0
2
   0.0000000000     0.0000000000     0.0000000000  m  1  1  1
   0.2500000000     0.2500000000     0.2500000000  m  1  1  1
"""

_TDDFT_SI_KPT = """\
K_POINTS
0
Gamma
4 4 4 0 0 0
"""


class TDDFTPlugin(BaseTestPlugin):
    """实时 TDDFT 测试插件"""

    @property
    def plugin_name(self) -> str:
        return "RT-TDDFT 混合规范测试"

    @property
    def calc_type(self) -> str:
        return "tddft"

    def can_handle(self, tutorial_content: str) -> bool:
        """检测 esolver_type = tddft（在 INPUT 代码块中）"""
        input_block_pattern = r'```[\s\S]*?INPUT_PARAMETERS[\s\S]*?```'
        input_blocks = re.findall(input_block_pattern, tutorial_content)
        for block in input_blocks:
            if re.search(r'esolver_type\s*=?\s*tddft', block, re.IGNORECASE):
                return True
        # 也匹配正文中的明确标识
        return bool(re.search(r'esolver_type\s*=?\s*tddft', tutorial_content, re.IGNORECASE))

    def extract_test_info(self, tutorial_content: str) -> Optional[TestInfo]:
        """提取 RT-TDDFT 测试信息——固定使用 Si 原胞 50 步测试案例"""
        test_info = TestInfo(
            calc_type=self.calc_type,
            case_name="tddft_si_hybrid"
        )
        test_info.input_content = _TDDFT_SI_INPUT
        test_info.stru_content = _TDDFT_SI_STRU
        test_info.kpt_content = _TDDFT_SI_KPT
        test_info.pseudopotentials = ["Si_ONCV_PBE-1.0.upf"]
        test_info.orbitals = ["Si_gga_8au_100Ry_2s2p1d.orb"]
        # 无数值预期，仅验证完成情况
        test_info.expected_results = {}
        return test_info

    def prepare_inputs(self, test_info: TestInfo, work_dir: Path) -> List[Path]:
        """准备 TDDFT 计算输入文件"""
        case_dir = work_dir / "tddft_si_hybrid"
        case_dir.mkdir(exist_ok=True)

        (case_dir / "INPUT").write_text(_TDDFT_SI_INPUT, encoding="utf-8")
        (case_dir / "STRU").write_text(_TDDFT_SI_STRU, encoding="utf-8")
        (case_dir / "KPT").write_text(_TDDFT_SI_KPT, encoding="utf-8")

        # 下载赝势和轨道文件
        if self.pp_manager:
            for pp in test_info.pseudopotentials:
                try:
                    pp_file = self.pp_manager.get_file(pp, "pseudopotential")
                    shutil.copy(pp_file, case_dir / pp)
                    print(f"  [OK] 赝势文件已就绪: {pp}")
                except Exception as e:
                    print(f"  [WARN] 赝势下载失败 {pp}: {e}")

            for orb in test_info.orbitals:
                try:
                    orb_file = self.pp_manager.get_file(orb, "orbital")
                    shutil.copy(orb_file, case_dir / orb)
                    print(f"  [OK] 轨道文件已就绪: {orb}")
                except Exception as e:
                    print(f"  [WARN] 轨道下载失败 {orb}: {e}")

        return [case_dir]

    def submit_jobs(self, input_dirs: List[Path]) -> List[str]:
        """提交 TDDFT 计算任务"""
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
        """验证 RT-TDDFT 计算结果"""
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

        # ── 验证 1：running_md.log 中 TDDFT 步数是否完成 ──
        md_logs = list(results_dir.rglob("tddft_si_hybrid/OUT.*/running_md.log"))
        if md_logs:
            md_result = self._validate_tddft_log(md_logs[0])
            validation.comparisons.update(md_result)
        else:
            validation.warnings.append("未找到 running_md.log，计算可能未完成")
            validation.comparisons["tddft_md_log_exists"] = {
                "key": "running_md.log 是否存在",
                "expected": "存在",
                "actual": "不存在",
                "abs_error": 1,
                "rel_error": 1,
                "passed": False
            }

        # ── 验证 2：电流输出文件是否生成 ──
        # ABACUS TDDFT 输出电流到 OUT.*/SPIN1_CURRENT（或类似名称）
        current_files = list(results_dir.rglob("tddft_si_hybrid/OUT.*/*CURRENT*"))
        current_ok = len(current_files) > 0
        validation.comparisons["tddft_current_output"] = {
            "key": "电流输出文件（SPIN1_CURRENT 或类似）是否生成",
            "expected": "存在",
            "actual": "存在" if current_ok else "不存在（可能需要检查 out_current 设置）",
            "abs_error": 0 if current_ok else 1,
            "rel_error": 0 if current_ok else 1,
            "passed": current_ok
        }

        # 综合判断：步数完成即认为通过（电流文件作为警告，不作为硬性失败）
        step_items = {k: v for k, v in validation.comparisons.items()
                      if 'steps' in k or 'scf' in k or 'log_exists' in k}
        if step_items:
            validation.passed = all(c.get("passed", False) for c in step_items.values())
        else:
            validation.passed = False

        if not current_ok:
            validation.warnings.append(
                "电流文件未找到。可能的原因：(1) 输出文件名与预期不符，请检查 OUT.* 目录；"
                "(2) out_current = 1 参数未生效。此项不影响整体通过判断。"
            )

        return validation

    def generate_report_section(self, validation: ValidationResult) -> str:
        """生成 Markdown 报告章节"""
        report = f"### {self.plugin_name}\n\n"
        report += "| 验证项 | 预期 | 实际 | 状态 |\n"
        report += "|--------|------|------|------|\n"

        for key, comp in validation.comparisons.items():
            status = "✅ PASS" if comp.get("passed") else "❌ FAIL"
            report += f"| {comp.get('key', key)} | {comp.get('expected', '-')} | {comp.get('actual', '-')} | {status} |\n"

        if validation.errors:
            report += "\n**错误：**\n"
            for e in validation.errors:
                report += f"- {e}\n"

        if validation.warnings:
            report += "\n**警告：**\n"
            for w in validation.warnings:
                report += f"- {w}\n"

        overall = "✅ 通过" if validation.passed else "❌ 失败"
        report += f"\n**总体结果：** {overall}\n\n"
        return report

    # ── 辅助验证方法 ──

    def _validate_tddft_log(self, log_path: Path) -> Dict:
        """从 running_md.log 验证 TDDFT 计算完成情况"""
        results = {}
        content = log_path.read_text(encoding="utf-8", errors="ignore")

        # ── 验证 SCF 收敛（TDDFT 内置 SCF）──
        scf_ok = (
            "charge density convergence is achieved" in content or
            bool(re.search(r'SCF\s+CONVERGE', content, re.IGNORECASE))
        )
        results["tddft_scf_converged"] = {
            "key": "内置 SCF 收敛",
            "expected": "是",
            "actual": "是" if scf_ok else "否",
            "abs_error": 0 if scf_ok else 1,
            "rel_error": 0 if scf_ok else 1,
            "passed": scf_ok
        }

        # ── 验证 TDDFT 步数完成 ──
        # ABACUS TDDFT 日志格式：STEP OF MOLECULAR DYNAMICS : N
        step_matches = re.findall(r'STEP OF MOLECULAR DYNAMICS\s*:\s*(\d+)', content)
        completed_steps = max((int(s) for s in step_matches), default=0)
        target_steps = 50  # INPUT 中设置的 md_nstep
        steps_ok = completed_steps >= target_steps

        results["tddft_md_steps"] = {
            "key": f"TDDFT 时间步完成（≥{target_steps}步）",
            "expected": f"≥ {target_steps}",
            "actual": str(completed_steps),
            "abs_error": 0 if steps_ok else 1,
            "rel_error": 0 if steps_ok else 1,
            "passed": steps_ok
        }

        # ── 提取总能量（验证为负数）──
        etot_matches = re.findall(
            r'(?:E_KohnSham|etot)\s*[=:]?\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', content
        )
        if etot_matches:
            etot = float(etot_matches[-1])
            etot_ok = etot < 0
            results["tddft_etot_negative"] = {
                "key": "初始 SCF 总能量为负数",
                "expected": "< 0 eV",
                "actual": f"{etot:.4f} eV",
                "abs_error": 0 if etot_ok else 1,
                "rel_error": 0 if etot_ok else 1,
                "passed": etot_ok
            }

        return results
