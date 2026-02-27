"""
tools/orbital_validator.py
扫描教程 Markdown 文件中的轨道文件名，验证并提供修正建议。

用法：
    # 只检查，不修改
    python tools/orbital_validator.py path/to/tutorial.md

    # 检查并自动修正（原文件被覆盖）
    python tools/orbital_validator.py path/to/tutorial.md --fix

    # 检查并输出修正后内容到新文件
    python tools/orbital_validator.py path/to/tutorial.md --fix --output path/to/fixed.md
"""

import re
import sys
import argparse
from pathlib import Path

# 导入共享数据库
sys.path.insert(0, str(Path(__file__).parent))
from orbital_db import validate_orbital


def find_orbital_filenames(content: str):
    """
    从 Markdown 内容中提取所有 .orb 文件名及其位置。
    Returns: list of (filename, line_number)
    """
    results = []
    for i, line in enumerate(content.splitlines(), start=1):
        matches = re.findall(r'\b(\w+_gga_\w+\.orb)\b', line)
        for m in matches:
            results.append((m, i))
    return results


def validate_tutorial(tutorial_path: str, fix: bool = False, output_path: str = None):
    """
    验证教程文件中的轨道文件名。

    Args:
        tutorial_path: 教程文件路径
        fix: 是否自动修正
        output_path: 修正后输出路径（None 表示覆盖原文件）

    Returns:
        issues: list of (filename, line, corrected, message)
        fixed_count: 自动修正的数量
    """
    path = Path(tutorial_path)
    content = path.read_text(encoding='utf-8')

    orb_files = find_orbital_filenames(content)
    if not orb_files:
        print(f"[OK] 未发现 .orb 文件名引用。")
        return [], 0

    issues = []
    fixed_count = 0
    fixed_content = content

    print(f"\n[扫描] 在 {path.name} 中发现 {len(orb_files)} 处 .orb 引用：\n")

    for filename, line_no in orb_files:
        is_valid, corrected, message = validate_orbital(filename)

        if is_valid:
            print(f"  第{line_no}行  OK  {filename}")
        else:
            issues.append((filename, line_no, corrected, message))
            if corrected:
                print(f"  第{line_no}行  NG  {filename}")
                print(f"             -> {corrected}  ({message})")
                if fix:
                    fixed_content = fixed_content.replace(filename, corrected)
                    fixed_count += 1
            else:
                print(f"  第{line_no}行  ??  {filename}")
                print(f"             -> {message}")

    if fix and fixed_count > 0:
        out_path = Path(output_path) if output_path else path
        out_path.write_text(fixed_content, encoding='utf-8')
        print(f"\n[修正] 已自动修正 {fixed_count} 处，保存到：{out_path}")
    elif issues:
        print(f"\n[摘要] 发现 {len(issues)} 处问题，{sum(1 for _, _, c, _ in issues if c)} 处可自动修正。")
        print(f"       使用 --fix 参数自动修正可修正的问题。")

    return issues, fixed_count


def main():
    parser = argparse.ArgumentParser(
        description="验证教程中的 ABACUS 轨道文件名"
    )
    parser.add_argument("tutorial", help="教程 Markdown 文件路径")
    parser.add_argument("--fix", action="store_true", help="自动修正已知错误文件名")
    parser.add_argument("--output", help="修正后输出到指定文件（默认覆盖原文件）")
    args = parser.parse_args()

    issues, fixed = validate_tutorial(args.tutorial, fix=args.fix, output_path=args.output)

    # 退出码：0=无问题或已全部修正，1=有未修正问题
    unfixed = sum(1 for _, _, c, _ in issues if not c or not args.fix)
    sys.exit(0 if unfixed == 0 else 1)


if __name__ == "__main__":
    main()
