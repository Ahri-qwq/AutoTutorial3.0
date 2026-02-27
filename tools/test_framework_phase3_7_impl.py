"""
AutoTutorial 3.0 - 测试脚本框架
Phase 3-7: 完整实现
"""

import subprocess
import json
import time
import zipfile
import re
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orbital_db import ORBITAL_CORRECTIONS, query_github_candidates
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime, timedelta
import shutil


class OrbitalNotFoundError(Exception):
    """
    轨道文件下载失败且无已知映射时抛出。
    携带 GitHub 候选列表，供 Claude 在 testCLAUDE.md 流程中处理。

    Attributes:
        filename: 请求的错误文件名
        candidates: GitHub 返回的该元素候选文件名列表（可能为空）
    """
    def __init__(self, filename: str, candidates: list):
        self.filename = filename
        self.candidates = candidates
        super().__init__(
            f"轨道文件不存在：{filename}\n"
            f"候选文件：{candidates if candidates else '（GitHub 查询失败，请手动提供）'}"
        )


class BohriumJobManager:
    """Bohrium任务管理器"""

    def __init__(self, project_id: int = 205855):
        self.project_id = project_id
        self.default_machine = "c16_m32_cpu"
        self.default_image = "registry.dp.tech/dptech/abacus:LTSv3.10.1"

    def create_job_config(self, job_name: str, input_dir: Path, log_file: str = "OUT.*/running_*.log") -> Dict:
        """创建job.json配置"""
        # 设置环境变量，确保ABACUS能找到赝势和轨道文件
        # 这是为了解决Bohrium容器中环境变量未设置的问题
        command = "export ABACUS_PP_PATH=./ && export ABACUS_ORB_PATH=./ && OMP_NUM_THREADS=1 mpirun -np 8 abacus"

        return {
            "job_name": job_name,
            "command": command,
            "log_file": log_file,
            "backward_files": [],
            "project_id": self.project_id,
            "platform": "ali",
            "job_type": "container",
            "machine_type": self.default_machine,
            "image_address": self.default_image
        }

    def submit_job(self, job_config: Dict, input_dir: Path) -> Optional[str]:
        """提交单个任务"""
        # 写入job.json
        job_file = input_dir / "job.json"
        with open(job_file, 'w', encoding='utf-8') as f:
            json.dump(job_config, f, indent=2)

        # 提交任务
        result = subprocess.run(
            ['bohr', 'job', 'submit', '-i', 'job.json', '-p', './'],
            cwd=str(input_dir),
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            # 解析Job ID
            match = re.search(r'JobId:\s+(\d+)', result.stdout)
            if match:
                job_id = match.group(1)
                print(f"  [OK] Job ID: {job_id}")
                return job_id
            else:
                print(f"  [ERROR] 无法解析Job ID")
                print(f"  输出: {result.stdout}")
                return None
        else:
            print(f"  [ERROR] 提交失败")
            print(f"  错误: {result.stderr}")
            return None

    def get_job_status(self, job_id: str) -> Dict:
        """获取任务状态"""
        result = subprocess.run(
            ['bohr', 'job', 'describe', '-j', job_id, '--json'],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )

        if result.returncode == 0 and result.stdout:
            try:
                data = json.loads(result.stdout)
                # bohr job describe 返回的是列表，取第一个元素
                if isinstance(data, list) and len(data) > 0:
                    data = data[0]

                return {
                    'status': data.get('statusStr', 'Unknown'),
                    'spend_time': data.get('spendTime', 0),
                    'error_info': data.get('errorInfo', '')
                }
            except (json.JSONDecodeError, KeyError, IndexError) as e:
                print(f"  [WARNING] 无法解析 Job {job_id} 的状态: {e}")
                return {'status': 'Unknown', 'spend_time': 0, 'error_info': 'JSON parse error'}
        else:
            return {'status': 'Unknown', 'spend_time': 0, 'error_info': 'Failed to get status'}

    def check_jobs_running(self, job_ids: List[str]) -> Dict:
        """检查任务是否正常运行"""
        status_summary = {
            'Running': 0,
            'Finished': 0,
            'Failed': 0,
            'Other': 0
        }

        for job_id in job_ids:
            status_info = self.get_job_status(job_id)
            status = status_info['status']

            if status in status_summary:
                status_summary[status] += 1
            else:
                status_summary['Other'] += 1

        return status_summary

    def download_job_result(self, job_id: str, output_dir: Path):
        """下载任务结果"""
        result_dir = output_dir / f"job_{job_id}"
        result_dir.mkdir(exist_ok=True, parents=True)

        result = subprocess.run(
            ['bohr', 'job', 'download', '-j', job_id, '-o', str(result_dir)],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            # 解压out.zip
            for zip_file in result_dir.rglob("out.zip"):
                with zipfile.ZipFile(zip_file, 'r') as zf:
                    zf.extractall(zip_file.parent)
                print(f"  [OK] Job {job_id} 结果已下载并解压")
            return True
        else:
            print(f"  [ERROR] Job {job_id} 下载失败: {result.stderr}")
            return False


class PseudopotentialManager:
    """赝势和轨道文件管理器（自动缓存）"""

    # 轨道文件映射表 - 从共享数据库导入，不要在此处添加新映射
    # 新的映射请添加到 tools/orbital_db.py 的 ORBITAL_CORRECTIONS 中
    ORBITAL_MAPPING = ORBITAL_CORRECTIONS

    def __init__(self, cache_dir: Path = Path("tools/pseudopotentials")):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(exist_ok=True, parents=True)

        self.orbital_cache_dir = Path("tools/orbitals")
        self.orbital_cache_dir.mkdir(exist_ok=True, parents=True)

    def get_file(self, filename: str, file_type: str = "pseudopotential") -> Path:
        """获取文件（自动下载并缓存）"""
        # 如果是轨道文件，检查是否需要映射
        original_filename = filename
        if file_type == "orbital" and filename in self.ORBITAL_MAPPING:
            filename = self.ORBITAL_MAPPING[filename]
            print(f"  [映射] {original_filename} -> {filename}")

        cache_dir = self.cache_dir if file_type == "pseudopotential" else self.orbital_cache_dir
        cache_file = cache_dir / filename

        # 如果缓存存在，验证格式
        if cache_file.exists():
            # 验证不是HTML文件或404错误
            try:
                with open(cache_file, 'rb') as f:
                    header = f.read(100)
                if (b'<html' in header.lower() or
                    b'<!doctype' in header.lower() or
                    b'404: Not Found' in header or
                    b'404 Not Found' in header):
                    print(f"  [警告] {filename} 是无效文件，重新下载")
                    try:
                        cache_file.unlink()
                    except PermissionError:
                        print(f"  [警告] 无法删除 {filename}，跳过")
                        pass
                else:
                    print(f"  [缓存] {filename}")
                    return cache_file
            except Exception as e:
                print(f"  [警告] 验证文件失败: {e}")
                pass

        # 下载文件
        print(f"  [下载] {filename}...")

        # 尝试多个下载源
        urls = [
            f"https://raw.githubusercontent.com/deepmodeling/abacus-develop/develop/tests/PP_ORB/{filename}",
            f"https://github.com/deepmodeling/abacus-develop/raw/develop/tests/PP_ORB/{filename}",
        ]

        for url in urls:
            result = subprocess.run(
                ['curl', '-L', '-o', str(cache_file), url],
                capture_output=True,
                text=True
            )

            if result.returncode == 0 and cache_file.exists():
                # 验证下载的文件不是HTML或404错误
                try:
                    with open(cache_file, 'rb') as f:
                        header = f.read(100)
                    if (b'<html' not in header.lower() and
                        b'<!doctype' not in header.lower() and
                        b'404: Not Found' not in header and
                        b'404 Not Found' not in header):
                        print(f"  [OK] {filename} 下载完成")
                        return cache_file
                    else:
                        print(f"  [失败] {url} 返回错误页面")
                        try:
                            cache_file.unlink()
                        except PermissionError:
                            pass
                except Exception as e:
                    print(f"  [警告] 验证下载文件失败: {e}")

        # 如果是轨道文件，查询 GitHub 候选列表并抛出结构化异常
        if file_type == "orbital":
            element = filename.split("_")[0] if "_" in filename else ""
            candidates = query_github_candidates(element) if element else []
            raise OrbitalNotFoundError(filename, candidates)
        raise Exception(f"下载失败: {filename} - 所有下载源都失败")

    def copy_to_dir(self, filename: str, target_dir: Path, file_type: str = "pseudopotential"):
        """复制文件到目标目录"""
        source = self.get_file(filename, file_type)
        target = target_dir / filename
        shutil.copy(source, target)


class ResultComparator:
    """结果对比器"""

    def __init__(self, tolerance: float = 0.05):
        self.tolerance = tolerance

    def extract_elastic_results(self, result_dir: Path) -> Dict:
        """从abacustest后处理结果提取弹性常数"""
        metrics_file = result_dir / "metrics_elastic.json"

        if not metrics_file.exists():
            return {}

        with open(metrics_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # 提取结果
        results = {}
        if 'elastic_tensor' in data:
            tensor = data['elastic_tensor']
            results['C11'] = tensor.get('C11', 0)
            results['C12'] = tensor.get('C12', 0)
            results['C44'] = tensor.get('C44', 0)

        if 'elastic_moduli' in data:
            moduli = data['elastic_moduli']
            results['bulk_modulus'] = moduli.get('bulk_modulus', 0)
            results['shear_modulus'] = moduli.get('shear_modulus', 0)
            results['young_modulus'] = moduli.get('young_modulus', 0)
            results['poisson_ratio'] = moduli.get('poisson_ratio', 0)

        return results

    def compare(self, expected: Dict, actual: Dict) -> Dict:
        """对比预期和实际结果"""
        comparison = {}

        for key in expected:
            if key not in actual:
                comparison[key] = {
                    'expected': expected[key],
                    'actual': None,
                    'relative_error': None,
                    'passed': False,
                    'message': '缺少实际结果'
                }
                continue

            exp_val = float(expected[key])
            act_val = float(actual[key])

            if exp_val == 0:
                rel_error = abs(act_val - exp_val)
                passed = rel_error <= self.tolerance
            else:
                rel_error = abs(act_val - exp_val) / abs(exp_val)
                passed = rel_error <= self.tolerance

            comparison[key] = {
                'expected': exp_val,
                'actual': act_val,
                'relative_error': rel_error,
                'passed': passed,
                'message': 'PASS' if passed else 'FAIL'
            }

        return comparison

    def generate_comparison_table(self, comparison: Dict) -> str:
        """生成对比表格（Markdown格式）"""
        lines = []
        lines.append("| 参数 | 预期值 | 实际值 | 相对误差 | 状态 |")
        lines.append("|------|--------|--------|---------|------|")

        for key, data in comparison.items():
            if data['actual'] is None:
                lines.append(f"| {key} | {data['expected']:.2f} | N/A | N/A | ❌ FAIL |")
            else:
                status = "✅ PASS" if data['passed'] else "❌ FAIL"
                error_pct = f"{data['relative_error']*100:.2f}%"
                lines.append(
                    f"| {key} | {data['expected']:.2f} | {data['actual']:.2f} | {error_pct} | {status} |"
                )

        return "\n".join(lines)


class ReportGenerator:
    """测试报告生成器"""

    def __init__(self, test_dir: Path):
        self.test_dir = test_dir

    def generate(self, analysis, job_ids: List[str], comparison: Dict, start_time: datetime, end_time: datetime):
        """生成完整测试报告"""
        report_file = self.test_dir / "test_report.md"

        # 计算统计
        total_tests = len(comparison)
        passed_tests = sum(1 for v in comparison.values() if v['passed'])
        failed_tests = total_tests - passed_tests
        pass_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0

        # 生成报告
        comparator = ResultComparator()
        comparison_table = comparator.generate_comparison_table(comparison)

        duration = end_time - start_time
        duration_str = str(duration).split('.')[0]  # 去掉微秒

        overall_status = "✅ PASS" if failed_tests == 0 else "❌ FAIL"

        report_content = f"""# 测试报告

## 教程信息
- **标题**: {analysis.title}
- **主题**: {analysis.topic}
- **测试时间**: {start_time.strftime('%Y-%m-%d %H:%M:%S')}
- **测试状态**: {overall_status}

## 测试配置
- **工具依赖**: {'ASE, ' if analysis.needs_ase else ''}{'abacustest' if analysis.needs_abacustest else ''}
- **提交任务数**: {len(job_ids)}
- **总运行时间**: {duration_str}

## 测试结果统计
- **总测试项**: {total_tests}
- **通过**: {passed_tests}
- **失败**: {failed_tests}
- **通过率**: {pass_rate:.1f}%

## 结果对比

{comparison_table}

## 任务详情
- **Job IDs**: {', '.join(job_ids)}
- **结果目录**: {self.test_dir / 'results'}

## 结论
{'✅ 教程内容准确，计算结果与预期一致。' if failed_tests == 0 else '❌ 部分结果与预期不符，需要检查教程内容。'}

---
*报告生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
"""

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)

        print(f"\n[OK] 测试报告已生成: {report_file}")
        return report_file
