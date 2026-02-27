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


def find_input_param_issues(content: str):
    """
    从 Markdown 内容中检测 INPUT 参数兼容性问题。
    Returns: list of (param_desc, line_number, fix_action)
      fix_action: 'remove_line' = 删除该行
    """
    results = []
    for i, line in enumerate(content.splitlines(), start=1):
        # nbands auto - ABACUS v3.10.x 不支持，删除后 ABACUS 自动计算
        if re.search(r'\bnbands\s+auto\b', line, re.IGNORECASE):
            results.append(('nbands auto', i, 'remove_line'))
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
    issues = []
    fixed_count = 0
    fixed_content = content

    if not orb_files:
        print(f"[OK] 未发现 .orb 文件名引用。")
    else:
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

    # ── INPUT 参数兼容性检查 ──────────────────────────────────────────
    param_issues = find_input_param_issues(content)

    if param_issues:
        print(f"\n[扫描] 发现 {len(param_issues)} 处 INPUT 参数兼容性问题：\n")
        if fix:
            lines = fixed_content.splitlines(keepends=True)
            lines_to_remove = set()
            for param_desc, line_no, action in param_issues:
                print(f"  第{line_no}行  NG  {param_desc}")
                print(f"             -> 删除此行（ABACUS 自动计算，无需指定）")
                if action == 'remove_line':
                    lines_to_remove.add(line_no)
            if lines_to_remove:
                fixed_content = ''.join(
                    line for i, line in enumerate(lines, start=1)
                    if i not in lines_to_remove
                )
                fixed_count += len(lines_to_remove)
        else:
            for param_desc, line_no, action in param_issues:
                print(f"  第{line_no}行  NG  {param_desc}")
                print(f"             -> 删除此行（ABACUS 自动计算，无需指定）")
            print(f"\n[摘要] 发现 {len(param_issues)} 处 INPUT 参数问题，均可自动修正。")
            print(f"       使用 --fix 参数自动修正。")

    # ── 最终保存 ─────────────────────────────────────────────────────
    if fix and fixed_count > 0:
        out_path = Path(output_path) if output_path else path
        out_path.write_text(fixed_content, encoding='utf-8')
        print(f"\n[修正] 已自动修正 {fixed_count} 处，保存到：{out_path}")
    elif issues or param_issues:
        fixable = sum(1 for _, _, c, _ in issues if c) + len(param_issues)
        print(f"\n[摘要] 发现问题：轨道文件 {len(issues)} 处，INPUT参数 {len(param_issues)} 处。")
        print(f"       共 {fixable} 处可自动修正，使用 --fix 参数修正。")

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
