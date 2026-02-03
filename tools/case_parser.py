#!/usr/bin/env python3
"""
案例解析工具 - 命令行版本
用于解析用户提供的案例文件（docx/md），提取关键信息

使用方法：
    python tools/case_parser.py --input "path/to/case.docx"
"""

import os
import sys
import argparse
import re
from pathlib import Path

# 添加父目录到路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, ROOT_DIR)


class CaseParser:
    """案例解析器"""

    def __init__(self, file_path):
        """
        初始化解析器

        Args:
            file_path: 案例文件路径
        """
        self.file_path = file_path
        self.content = ""
        self.file_type = Path(file_path).suffix.lower()

    def read_file(self):
        """读取文件内容"""
        if not os.path.exists(self.file_path):
            print(f"[错误] 文件不存在: {self.file_path}", file=sys.stderr)
            sys.exit(1)

        try:
            if self.file_type == ".docx":
                self.content = self._read_docx()
            elif self.file_type in [".md", ".txt"]:
                self.content = self._read_text()
            else:
                print(f"[错误] 不支持的文件格式: {self.file_type}", file=sys.stderr)
                print("[提示] 支持的格式: .docx, .md, .txt", file=sys.stderr)
                sys.exit(1)
        except Exception as e:
            print(f"[错误] 读取文件失败: {e}", file=sys.stderr)
            sys.exit(1)

    def _read_docx(self):
        """读取docx文件"""
        try:
            import docx
        except ImportError:
            print("[错误] 缺少python-docx库", file=sys.stderr)
            print("[提示] 请安装: pip install python-docx", file=sys.stderr)
            sys.exit(1)

        doc = docx.Document(self.file_path)
        paragraphs = [para.text for para in doc.paragraphs]
        return "\n".join(paragraphs)

    def _read_text(self):
        """读取文本文件"""
        with open(self.file_path, 'r', encoding='utf-8', errors='ignore') as f:
            return f.read()

    def parse(self):
        """
        解析案例内容

        Returns:
            dict: 解析结果
        """
        result = {
            "file_structure": self._extract_file_structure(),
            "parameters": self._extract_parameters(),
            "workflow": self._extract_workflow(),
            "special_settings": self._extract_special_settings(),
        }
        return result

    def _extract_file_structure(self):
        """提取文件结构"""
        files = []

        # 常见的ABACUS文件
        common_files = ["INPUT", "STRU", "KPT", "KPOINTS", "CONTROL", "NUMERICAL_ORBITAL"]

        for file_name in common_files:
            # 查找文件名（不区分大小写）
            pattern = re.compile(rf'\b{file_name}\b', re.IGNORECASE)
            if pattern.search(self.content):
                files.append(file_name)

        # 查找其他可能的文件（如赝势文件、轨道文件）
        # 赝势文件：.upf, .UPF
        upf_pattern = re.compile(r'\b\w+\.upf\b', re.IGNORECASE)
        upf_files = upf_pattern.findall(self.content)
        if upf_files:
            files.append(f"赝势文件 ({len(set(upf_files))}个)")

        # 轨道文件：.orb
        orb_pattern = re.compile(r'\b\w+\.orb\b', re.IGNORECASE)
        orb_files = orb_pattern.findall(self.content)
        if orb_files:
            files.append(f"轨道文件 ({len(set(orb_files))}个)")

        return files if files else ["未识别到文件结构"]

    def _extract_parameters(self):
        """提取关键参数"""
        parameters = {}

        # 常见的ABACUS参数模式
        # 格式1: parameter = value
        # 格式2: parameter value
        # 格式3: parameter: value

        # 定义关键参数列表
        key_params = [
            "calculation", "basis_type", "ecutwfc", "ecutrho",
            "scf_thr", "scf_nmax", "cal_stress", "cal_force",
            "stress_thr", "force_thr", "relax_nmax",
            "smearing_method", "smearing_sigma", "mixing_type",
            "mixing_beta", "ks_solver", "nbands", "nspin",
            "ntype", "nelec", "symmetry", "kspacing"
        ]

        for param in key_params:
            # 尝试多种模式匹配
            patterns = [
                rf'{param}\s*=\s*([^\s\n]+)',  # parameter = value
                rf'{param}\s+([^\s\n]+)',       # parameter value
                rf'{param}\s*:\s*([^\s\n]+)',   # parameter: value
            ]

            for pattern in patterns:
                match = re.search(pattern, self.content, re.IGNORECASE)
                if match:
                    parameters[param] = match.group(1).strip()
                    break

        return parameters if parameters else {"提示": "未识别到参数，可能需要手动提取"}

    def _extract_workflow(self):
        """提取计算流程"""
        workflow = []

        # 查找常见的流程关键词
        workflow_keywords = {
            "结构优化": ["relax", "optimization", "优化", "结构优化"],
            "自洽计算": ["scf", "self-consistent", "自洽"],
            "非自洽计算": ["nscf", "non-scf", "非自洽"],
            "能带计算": ["band", "能带"],
            "态密度计算": ["dos", "density of states", "态密度"],
            "应力计算": ["stress", "应力"],
            "力计算": ["force", "受力"],
            "分子动力学": ["md", "molecular dynamics", "分子动力学"],
        }

        for step_name, keywords in workflow_keywords.items():
            for keyword in keywords:
                if re.search(rf'\b{keyword}\b', self.content, re.IGNORECASE):
                    workflow.append(step_name)
                    break

        # 尝试提取编号的步骤
        numbered_steps = re.findall(r'(?:步骤|Step)\s*[：:]\s*(.+?)(?:\n|$)', self.content, re.IGNORECASE)
        if numbered_steps:
            workflow.extend(numbered_steps[:5])  # 最多提取5个步骤

        return workflow if workflow else ["未识别到明确的计算流程"]

    def _extract_special_settings(self):
        """提取特殊设置"""
        special = []

        # 查找特殊设置的关键词
        special_keywords = [
            ("高精度", ["high accuracy", "高精度", "精确"]),
            ("自定义阈值", ["threshold", "阈值", "thr"]),
            ("特殊k点", ["k-point", "k点", "kspacing"]),
            ("自旋极化", ["spin", "nspin", "自旋"]),
            ("范德华修正", ["vdw", "van der waals", "范德华"]),
            ("DFT+U", ["dft+u", "hubbard", "u值"]),
        ]

        for setting_name, keywords in special_keywords:
            for keyword in keywords:
                if re.search(rf'\b{keyword}\b', self.content, re.IGNORECASE):
                    special.append(setting_name)
                    break

        # 查找注释或说明中的特殊设置
        comment_pattern = re.compile(r'(?:注意|注|说明|备注)[：:]\s*(.+?)(?:\n|$)', re.IGNORECASE)
        comments = comment_pattern.findall(self.content)
        if comments:
            special.extend([f"说明: {c.strip()}" for c in comments[:3]])  # 最多3条

        return special if special else ["未识别到特殊设置"]


def format_output(result, file_path):
    """格式化输出结果"""
    # 设置输出编码为UTF-8（Windows兼容）
    if sys.platform == 'win32':
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    print("=" * 80)
    print("案例解析结果")
    print("=" * 80)
    print()
    print(f"文件: {file_path}")
    print()

    print("## 文件结构")
    for item in result["file_structure"]:
        print(f"- {item}")
    print()

    print("## 关键参数")
    if isinstance(result["parameters"], dict):
        if "提示" in result["parameters"]:
            print(f"[提示] {result['parameters']['提示']}")
        else:
            for param, value in result["parameters"].items():
                print(f"{param} = {value}")
    print()

    print("## 计算流程")
    for i, step in enumerate(result["workflow"], 1):
        print(f"{i}. {step}")
    print()

    print("## 特殊设置")
    for item in result["special_settings"]:
        print(f"- {item}")
    print()


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="案例解析工具 - 解析ABACUS案例文件",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    python tools/case_parser.py --input "data/input/example.docx"
    python tools/case_parser.py --input "data/input/example.md"
        """
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="案例文件路径 (支持 .docx, .md, .txt)"
    )

    args = parser.parse_args()

    # 初始化解析器
    case_parser = CaseParser(args.input)

    # 读取文件
    case_parser.read_file()

    # 解析内容
    result = case_parser.parse()

    # 输出结果
    format_output(result, args.input)


if __name__ == "__main__":
    main()
