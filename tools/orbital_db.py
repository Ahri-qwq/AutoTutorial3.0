"""
tools/orbital_db.py
轨道文件可信数据库 - 单一可信数据源

KNOWN_ORBITALS: 已在 ABACUS GitHub tests/PP_ORB/ 验证存在的轨道文件
ORBITAL_CORRECTIONS: 错误文件名 → 正确文件名 的纠错映射

维护方式：
  当测试发现新的不存在文件时，在 ORBITAL_CORRECTIONS 中添加映射。
  当从 ABACUS GitHub 新增验证文件时，在 KNOWN_ORBITALS 中添加。
"""

# 已验证存在于 ABACUS GitHub tests/PP_ORB/ 的轨道文件
# 来源：tools/orbitals/ 缓存目录中的实际下载文件 + ABACUS 官方仓库
KNOWN_ORBITALS = {
    "Si": [
        "Si_gga_8au_100Ry_2s2p1d.orb",
        "Si_gga_8au_60Ry_2s2p1d.orb",
    ],
    "O": [
        "O_gga_7au_100Ry_2s2p1d.orb",
    ],
    "Ti": [
        "Ti_gga_8au_100Ry_4s2p2d1f.orb",
    ],
    "C": [
        "C_gga_8au_100Ry_2s2p1d.orb",
    ],
    "H": [
        "H_gga_8au_100Ry_2s1p.orb",
    ],
    "Mg": [
        "Mg_gga_8au_100Ry_4s2p1d.orb",
    ],
    "Fe": [
        "Fe_gga_9au_100Ry_4s2p2d1f.orb",
    ],
    "Al": [
        "Al_gga_7au_100Ry_4s4p1d.orb",
    ],
}

# 纠错映射表：错误文件名 → 正确文件名
# 来源：测试失败记录 + ABACUS 官方仓库验证
ORBITAL_CORRECTIONS = {
    # Si: RAG 知识库中的 7au 版本实际不存在，应使用 8au
    "Si_gga_7au_100Ry_2s2p1d.orb": "Si_gga_8au_100Ry_2s2p1d.orb",
    # Ti: 2s2p2d1f 版本不存在，应使用 4s2p2d1f
    "Ti_gga_7au_100Ry_2s2p2d1f.orb": "Ti_gga_8au_100Ry_4s2p2d1f.orb",
    "Ti_gga_8au_100Ry_2s2p2d1f.orb": "Ti_gga_8au_100Ry_4s2p2d1f.orb",
    # O: 6au 版本不存在，应使用 7au
    "O_gga_6au_100Ry_2s2p1d.orb": "O_gga_7au_100Ry_2s2p1d.orb",
    # C: 7au 版本不存在，应使用 8au
    "C_gga_7au_100Ry_2s2p1d.orb": "C_gga_8au_100Ry_2s2p1d.orb",
}

# 扁平化已知文件名集合（用于快速查询）
ALL_KNOWN_ORBITALS = {f for files in KNOWN_ORBITALS.values() for f in files}


def validate_orbital(filename: str):
    """
    验证轨道文件名是否正确。

    Returns:
        (is_valid: bool, corrected: str, message: str)
        - is_valid=True: 文件名已知且正确
        - is_valid=False, corrected!=None: 有纠错建议
        - is_valid=False, corrected=None: 未知文件名，无法自动修正
    """
    if filename in ALL_KNOWN_ORBITALS:
        return True, filename, "OK"

    if filename in ORBITAL_CORRECTIONS:
        corrected = ORBITAL_CORRECTIONS[filename]
        return False, corrected, f"已知错误文件名，应改为 {corrected}"

    # 尝试从文件名推断元素并给出建议
    element = filename.split("_")[0] if "_" in filename else None
    if element and element in KNOWN_ORBITALS:
        suggestions = KNOWN_ORBITALS[element]
        return False, None, f"未知文件名。{element} 的已知轨道文件：{', '.join(suggestions)}"

    return False, None, f"未知文件名，无法自动修正。请到 https://github.com/deepmodeling/abacus-develop/tree/develop/tests/PP_ORB 查找正确文件名"


def query_github_candidates(element: str) -> list:
    """
    查询 ABACUS GitHub，返回该元素的所有已知 .orb 文件名。

    Args:
        element: 元素符号，如 "Si"、"Ti"

    Returns:
        list of str: 文件名列表，网络失败时返回空列表
    """
    import urllib.request
    import json

    url = "https://api.github.com/repos/deepmodeling/abacus-develop/contents/tests/PP_ORB"
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "AutoTutorial3.0-OrbitalChecker"}
        )
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
        return [
            item["name"]
            for item in data
            if item["name"].startswith(element + "_") and item["name"].endswith(".orb")
        ]
    except Exception:
        # 任何失败（网络、rate limit、API 结构变更）均返回空列表
        # 空列表触发 testCLAUDE.md 2.3b 的「手动输入」分支，不会静默丢失错误
        return []


def add_correction(wrong: str, correct: str):
    """
    向 ORBITAL_CORRECTIONS 添加新映射，同时更新内存状态和源文件。

    Args:
        wrong: 错误的文件名（如 "Si_gga_7au_100Ry_2s2p1d.orb"）
        correct: 正确的文件名（如 "Si_gga_8au_100Ry_2s2p1d.orb"）
    """
    # 前置校验 1：wrong 与 correct 相同时无意义，拒绝写入
    if wrong == correct:
        raise ValueError(f"wrong 和 correct 相同：{wrong}，无需添加映射")
    # 前置校验 2：key 已存在时内存已是最新，跳过文件写入（幂等保证）
    if wrong in ORBITAL_CORRECTIONS:
        return

    from pathlib import Path as _Path

    # 1. 更新内存状态
    ORBITAL_CORRECTIONS[wrong] = correct
    ALL_KNOWN_ORBITALS.add(correct)

    # 2. 读取源文件
    db_path = _Path(__file__).resolve()
    content = db_path.read_text(encoding='utf-8')

    # 3. 在 ORBITAL_CORRECTIONS 结束的 } 前插入新条目
    #    文件中的标记：}\n\n# 扁平化已知文件名集合
    element = wrong.split("_")[0] if "_" in wrong else "unknown"
    new_entry = f'    # {element}: 运行时自动发现（测试阶段用户确认）\n    "{wrong}": "{correct}",\n'
    marker = '}\n\n# 扁平化已知文件名集合'
    if marker not in content:
        raise RuntimeError("orbital_db.py 格式已变化，无法自动插入新映射，请手动添加。")
    content = content.replace(marker, new_entry + marker, 1)

    # 4. 写回文件
    db_path.write_text(content, encoding='utf-8')
    print(f"[orbital_db] 已写入新映射：{wrong} -> {correct}")
