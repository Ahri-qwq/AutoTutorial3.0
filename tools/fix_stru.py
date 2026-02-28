"""
STRU文件路径修复工具
将STRU文件中的赝势和轨道文件路径修改为相对路径
解决Bohrium容器中环境变量未设置的问题
"""

from pathlib import Path
import re
import sys


def fix_stru_paths(stru_file: Path, verbose: bool = True) -> bool:
    """
    修改STRU文件中的赝势和轨道路径为相对路径

    Args:
        stru_file: STRU文件路径
        verbose: 是否打印详细信息

    Returns:
        bool: 是否成功修改
    """
    if not stru_file.exists():
        if verbose:
            print(f"[ERROR] STRU文件不存在: {stru_file}")
        return False

    try:
        # 读取文件内容
        content = stru_file.read_text(encoding='utf-8')
        original_content = content

        # 分行处理
        lines = content.split('\n')
        fixed_lines = []
        modifications = []

        in_atomic_species = False
        in_numerical_orbital = False

        for i, line in enumerate(lines):
            original_line = line

            # 检测ATOMIC_SPECIES块
            if 'ATOMIC_SPECIES' in line:
                in_atomic_species = True
                fixed_lines.append(line)
                continue

            # 检测NUMERICAL_ORBITAL块
            if 'NUMERICAL_ORBITAL' in line:
                in_numerical_orbital = True
                fixed_lines.append(line)
                continue

            # 处理ATOMIC_SPECIES块中的赝势路径
            if in_atomic_species and ('.upf' in line.lower() or '.UPF' in line):
                parts = line.split()
                if len(parts) >= 3:
                    # 格式: Element Mass PseudopotentialPath
                    old_path = parts[2]
                    if '/' in old_path or '\\' in old_path:
                        # 提取文件名
                        filename = Path(old_path).name
                        new_path = f'./{filename}'
                        parts[2] = new_path
                        line = ' '.join(parts)
                        modifications.append(f"  行{i+1}: {old_path} -> {new_path}")

            # 处理NUMERICAL_ORBITAL块中的轨道路径
            if in_numerical_orbital and ('.orb' in line.lower() or '.ORB' in line):
                old_path = line.strip()
                if '/' in old_path or '\\' in old_path:
                    # 提取文件名
                    filename = Path(old_path).name
                    new_path = f'./{filename}'
                    # 保持原有的缩进
                    indent = len(line) - len(line.lstrip())
                    line = ' ' * indent + new_path
                    modifications.append(f"  行{i+1}: {old_path} -> {new_path}")

            # 检测块结束（遇到新的大写关键字）
            if line.strip() and not line.startswith(' ') and not line.startswith('\t'):
                if line.strip()[0].isupper() and 'ATOMIC_SPECIES' not in line and 'NUMERICAL_ORBITAL' not in line:
                    in_atomic_species = False
                    in_numerical_orbital = False

            fixed_lines.append(line)

        # 如果有修改，写回文件
        if modifications:
            new_content = '\n'.join(fixed_lines)
            stru_file.write_text(new_content, encoding='utf-8')

            if verbose:
                print(f"[OK] STRU文件路径已修复: {stru_file.name}")
                print(f"[修改] 共修改 {len(modifications)} 处路径:")
                for mod in modifications:
                    print(mod)

            return True
        else:
            if verbose:
                print(f"[INFO] STRU文件无需修改: {stru_file.name}")
            return True

    except Exception as e:
        if verbose:
            print(f"[ERROR] 修复STRU文件失败: {e}")
        return False


def fix_input_paths(input_file: Path, verbose: bool = True) -> bool:
    """
    修改INPUT文件中的赝势和轨道目录路径为相对路径

    Args:
        input_file: INPUT文件路径
        verbose: 是否打印详细信息

    Returns:
        bool: 是否成功修改
    """
    if not input_file.exists():
        if verbose:
            print(f"[ERROR] INPUT文件不存在: {input_file}")
        return False

    try:
        # 读取文件内容
        content = input_file.read_text(encoding='utf-8')
        lines = content.split('\n')
        fixed_lines = []
        modifications = []

        for i, line in enumerate(lines):
            original_line = line

            # 检查pseudo_dir和orbital_dir
            if 'pseudo_dir' in line.lower():
                parts = line.split()
                if len(parts) >= 2:
                    old_path = parts[1]
                    if old_path != './':
                        new_path = './'
                        parts[1] = new_path
                        line = ' '.join(parts)
                        modifications.append(f"  行{i+1}: pseudo_dir {old_path} -> {new_path}")

            if 'orbital_dir' in line.lower():
                parts = line.split()
                if len(parts) >= 2:
                    old_path = parts[1]
                    if old_path != './':
                        new_path = './'
                        parts[1] = new_path
                        line = ' '.join(parts)
                        modifications.append(f"  行{i+1}: orbital_dir {old_path} -> {new_path}")

            fixed_lines.append(line)

        # 如果有修改，写回文件
        if modifications:
            new_content = '\n'.join(fixed_lines)
            input_file.write_text(new_content, encoding='utf-8')

            if verbose:
                print(f"[OK] INPUT文件路径已修复: {input_file.name}")
                print(f"[修改] 共修改 {len(modifications)} 处路径:")
                for mod in modifications:
                    print(mod)

            return True
        else:
            if verbose:
                print(f"[INFO] INPUT文件无需修改: {input_file.name}")
            return True

    except Exception as e:
        if verbose:
            print(f"[ERROR] 修复INPUT文件失败: {e}")
        return False


def fix_directory(directory: Path, verbose: bool = True) -> bool:
    """
    修复目录中所有STRU和INPUT文件的路径

    Args:
        directory: 目录路径
        verbose: 是否打印详细信息

    Returns:
        bool: 是否全部成功
    """
    if not directory.exists():
        if verbose:
            print(f"[ERROR] 目录不存在: {directory}")
        return False

    success = True

    # 修复STRU文件
    stru_files = list(directory.glob("**/STRU"))
    if stru_files:
        if verbose:
            print(f"\n[处理] 找到 {len(stru_files)} 个STRU文件")
        for stru_file in stru_files:
            if not fix_stru_paths(stru_file, verbose):
                success = False

    # 修复INPUT文件
    input_files = list(directory.glob("**/INPUT"))
    if input_files:
        if verbose:
            print(f"\n[处理] 找到 {len(input_files)} 个INPUT文件")
        for input_file in input_files:
            if not fix_input_paths(input_file, verbose):
                success = False

    return success


def main():
    """命令行入口"""
    if len(sys.argv) < 2:
        print("用法:")
        print("  python fix_stru.py <STRU文件路径>")
        print("  python fix_stru.py <目录路径>")
        print("\n示例:")
        print("  python fix_stru.py test_workspace/Si_relax/STRU")
        print("  python fix_stru.py test_workspace/")
        sys.exit(1)

    path = Path(sys.argv[1])

    if path.is_file():
        # 单个文件
        if path.name == 'STRU':
            success = fix_stru_paths(path)
        elif path.name == 'INPUT':
            success = fix_input_paths(path)
        else:
            print(f"[ERROR] 不支持的文件类型: {path.name}")
            success = False
    elif path.is_dir():
        # 目录
        success = fix_directory(path)
    else:
        print(f"[ERROR] 路径不存在: {path}")
        success = False

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
