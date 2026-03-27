#!/usr/bin/env python3
"""
static_checker.py — judgeCLAUDE 静态指标检查脚本

检查 7 项可客观判定的教程质量指标：
  1.3 信息精准与一致性   (max 5)
  1.4 文件来源明确度     (max 3)
  1.5 时效与兼容性       (max 4)
  2.4 可导航性与步骤连贯性 (max 1)
  2.5 结构完整性         (max 5)
  3.2 环境条件明确度     (max 3)
  4.1 图文/代码邻近性    (max 3)

用法：
  python tools/static_checker.py --input <tutorial.md> --output <scores.json>
"""

import argparse
import json
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def load_markdown(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read()


def split_code_blocks(text: str):
    """返回 (prose_text, [(lang, code_content), ...])"""
    code_blocks = []
    prose_parts = []
    pos = 0
    for m in re.finditer(r"```(\w*)\n(.*?)```", text, re.DOTALL):
        prose_parts.append(text[pos:m.start()])
        code_blocks.append((m.group(1), m.group(2)))
        pos = m.end()
    prose_parts.append(text[pos:])
    prose = "\n".join(prose_parts)
    return prose, code_blocks


# ---------------------------------------------------------------------------
# 1.3 信息精准与一致性（满分 5）
# ---------------------------------------------------------------------------
# 策略：提取正文中明确的"参数 = 值"表述，再在代码块中找同名参数，
# 比对数值。不一致条数越多，扣分越多。

PARAM_PROSE_RE = re.compile(
    r"(?:设置?|建议?|推荐|设为|取|使用|取值为?)?\s*"
    r"`?(\w+)`?\s*[=＝为设]\s*([0-9]+(?:\.[0-9]+)?)",
    re.UNICODE,
)
PARAM_CODE_RE = re.compile(r"^\s*(\w+)\s*[=\s]\s*([0-9]+(?:\.[0-9]+)?)\s*(?:#.*)?$",
                            re.MULTILINE)


def check_1_3(prose: str, code_blocks) -> dict:
    prose_params: dict[str, set] = {}
    for m in PARAM_PROSE_RE.finditer(prose):
        k, v = m.group(1).lower(), m.group(2)
        prose_params.setdefault(k, set()).add(v)

    code_params: dict[str, set] = {}
    for _, code in code_blocks:
        for m in PARAM_CODE_RE.finditer(code):
            k, v = m.group(1).lower(), m.group(2)
            code_params.setdefault(k, set()).add(v)

    conflicts = []
    for k in prose_params:
        if k in code_params:
            pv = prose_params[k]
            cv = code_params[k]
            if not pv & cv:  # 无交集 → 不一致
                conflicts.append(f"{k}: 正文{pv} vs 代码{cv}")

    n_conflict = len(conflicts)
    if n_conflict == 0:
        score, note = 5, "无明显数值冲突"
    elif n_conflict == 1:
        score, note = 4, f"1处冲突：{conflicts[0]}"
    elif n_conflict == 2:
        score, note = 3, f"2处冲突：{'; '.join(conflicts)}"
    elif n_conflict <= 4:
        score, note = 1, f"{n_conflict}处冲突"
    else:
        score, note = 0, f"{n_conflict}处冲突（严重）"

    return {"score": score, "max": 5, "conflicts": conflicts, "note": note}


# ---------------------------------------------------------------------------
# 1.4 文件来源明确度（满分 3）
# ---------------------------------------------------------------------------
# 检查：赝势/轨道/结构文件的下载途径是否有交代。
# 关键词：dplibrary / SG15 / ONCV / PBE / 官网 / Github / deepmodeling
#         wget / 下载 / abacus-develop / PP_ORB

FILE_SOURCE_KEYWORDS = [
    r"dplibrary", r"SG15", r"ONCV", r"deepmodeling",
    r"abacus-develop", r"PP_ORB", r"ABACUS-orbitals",
    r"wget\s+http", r"下载", r"官网", r"github\.com",
    r"赝势.*来源", r"轨道.*来源", r"伪势.*下载",
    r"pseudopotential.*from", r"orbital.*from",
]

FILE_TYPE_KEYWORDS = {
    "pseudopotential": [r"\.upf\b", r"赝势", r"伪势", r"POTCAR", r"pseudopotential"],
    "orbital":         [r"\.orb\b", r"轨道文件", r"numerical.*orbital"],
    "structure":       [r"\.cif\b", r"STRU\b", r"结构文件", r"晶体结构"],
}


def check_1_4(text: str) -> dict:
    text_lower = text.lower()
    source_found = any(re.search(kw, text, re.IGNORECASE) for kw in FILE_SOURCE_KEYWORDS)

    mentioned_types = []
    for ftype, patterns in FILE_TYPE_KEYWORDS.items():
        if any(re.search(p, text, re.IGNORECASE) for p in patterns):
            mentioned_types.append(ftype)

    if source_found and len(mentioned_types) >= 2:
        score, note = 3, "赝势/轨道/结构来源均有交代"
    elif source_found and len(mentioned_types) == 1:
        score, note = 2, f"仅部分文件类型有来源说明：{mentioned_types}"
    elif source_found:
        score, note = 2, "有文件来源说明，但覆盖文件类型不全"
    elif len(mentioned_types) >= 1:
        score, note = 1, "提及了文件类型，但未说明获取途径"
    else:
        score, note = 0, "未交代任何文件来源"

    return {"score": score, "max": 3, "mentioned_types": mentioned_types,
            "source_found": source_found, "note": note}


# ---------------------------------------------------------------------------
# 1.5 时效与兼容性（满分 4）
# ---------------------------------------------------------------------------
# 检查废弃参数列表（ABACUS 已知废弃项）

DEPRECATED_PARAMS = [
    (r"\bnbands\s+auto\b",       "nbands auto（v3.10+ 不支持，应删除）"),
    (r"\bforce_thr\b",           "force_thr（旧版参数，建议用 force_thr_ev）"),
    (r"\bpseudopot_ntype\b",     "pseudopot_ntype（已废弃）"),
    (r"\bout_charge\s+[12]\b",   "out_charge 1/2（请改用 out_chg）"),
    (r"\bout_potential\s+[12]\b","out_potential（旧参数）"),
    (r"\bread_file_dir\b",       "read_file_dir（已废弃）"),
    (r"\bpara_k_num\b",          "para_k_num（旧版并行参数）"),
    (r"\bscaling_factor\b",      "scaling_factor（旧STRU格式参数）"),
]


def check_1_5(text: str) -> dict:
    found = []
    for pattern, desc in DEPRECATED_PARAMS:
        if re.search(pattern, text, re.IGNORECASE):
            found.append(desc)

    if len(found) == 0:
        score, note = 4, "未发现废弃参数"
    elif len(found) == 1:
        score, note = 2, f"发现1处废弃参数：{found[0]}"
    elif len(found) == 2:
        score, note = 1, f"发现{len(found)}处废弃参数"
    else:
        score, note = 0, f"发现{len(found)}处废弃参数（严重）"

    return {"score": score, "max": 4, "deprecated_found": found, "note": note}


# ---------------------------------------------------------------------------
# 2.4 可导航性与步骤连贯性（满分 1）
# ---------------------------------------------------------------------------
# 检查：
#   a) 是否有目录（TOC）或标题层级合理（不越级）
#   b) 步骤未过度分割（连续 H3/H4 标题超过 10 个视为过度分割）

HEADING_RE = re.compile(r"^(#{1,6})\s+(.+)", re.MULTILINE)


def check_2_4(text: str) -> dict:
    headings = HEADING_RE.findall(text)
    levels = [len(h) for h, _ in headings]

    # TOC：文中是否有 [xxx](#xxx) 锚点跳转或"目录"/"Table of Contents"字样
    has_toc = bool(re.search(r"(?:目录|Table of Contents|\[.*?\]\(#)", text, re.IGNORECASE))

    # 越级检查：从 H1 直接跳到 H3+（中间跳过一级）
    skip_level = False
    for i in range(1, len(levels)):
        if levels[i] - levels[i - 1] > 1:
            skip_level = True
            break

    # 过度分割：连续出现 >8 个 H3/H4 标题
    fine_headings = sum(1 for l in levels if l >= 3)
    over_split = fine_headings > 8

    issues = []
    if skip_level:
        issues.append("存在标题越级（如 H1→H3）")
    if over_split:
        issues.append(f"步骤过度分割（H3/H4 共{fine_headings}个）")

    if not issues:
        score, note = 1, "标题结构合理" + ("，有TOC" if has_toc else "")
    else:
        score, note = 0, "；".join(issues)

    return {"score": score, "max": 1, "has_toc": has_toc,
            "skip_level": skip_level, "over_split": over_split, "note": note}


# ---------------------------------------------------------------------------
# 2.5 结构完整性（满分 5）
# ---------------------------------------------------------------------------
# 必要模块：背景/引言、参数配置/输入文件、运行/执行、结果/输出、总结/附录

REQUIRED_SECTIONS = [
    ("background",  ["背景", "引言", "简介", "introduction", "overview", "前言",
                     "理论", "方法", "theory", "background"]),
    ("config",      ["配置", "输入", "参数", "INPUT", "STRU", "KPT",
                     "setup", "configuration", "设置", "准备"]),
    ("execution",   ["运行", "执行", "计算", "calculation", "run",
                     "任务", "提交", "job"]),
    ("results",     ["结果", "输出", "result", "output", "分析",
                     "数据", "验证", "收敛"]),
    ("summary",     ["总结", "小结", "附录", "进阶", "summary",
                     "conclusion", "reference", "参考"]),
]


def check_2_5(text: str) -> dict:
    found_sections = []
    missing_sections = []
    text_lower = text.lower()

    for sec_id, keywords in REQUIRED_SECTIONS:
        if any(kw.lower() in text_lower for kw in keywords):
            found_sections.append(sec_id)
        else:
            missing_sections.append(sec_id)

    score = len(found_sections)
    if missing_sections:
        note = f"缺失模块：{missing_sections}"
    else:
        note = "5个必要模块均齐全"

    return {"score": score, "max": 5, "found": found_sections,
            "missing": missing_sections, "note": note}


# ---------------------------------------------------------------------------
# 3.2 环境条件明确度（满分 3）
# ---------------------------------------------------------------------------
# 检查：版本号、并行/MPI 说明、资源需求（内存/时间/核数）

VERSION_RE = re.compile(
    r"(?:ABACUS|v|version|版本)\s*(?:v|V)?\s*\d+\.\d+(?:\.\d+)?",
    re.IGNORECASE,
)
MPI_RE = re.compile(
    r"(?:mpirun|mpiexec|MPI|OMP|omp_num_threads|并行|并发|核数|nprocs|-np\s+\d+)",
    re.IGNORECASE,
)
RESOURCE_RE = re.compile(
    r"(?:内存|memory|RAM|GB|MB|核|core|CPU|hour|小时|分钟|min|资源|resource)",
    re.IGNORECASE,
)


def check_3_2(text: str) -> dict:
    has_version  = bool(VERSION_RE.search(text))
    has_mpi      = bool(MPI_RE.search(text))
    has_resource = bool(RESOURCE_RE.search(text))

    score = sum([has_version, has_mpi, has_resource])
    missing = []
    if not has_version:  missing.append("版本号")
    if not has_mpi:      missing.append("MPI/并行说明")
    if not has_resource: missing.append("资源需求")

    note = "完整交代" if score == 3 else f"缺失：{missing}"
    return {"score": score, "max": 3,
            "has_version": has_version, "has_mpi": has_mpi,
            "has_resource": has_resource, "note": note}


# ---------------------------------------------------------------------------
# 4.1 图文/代码邻近性（满分 3）
# ---------------------------------------------------------------------------
# 策略：将文本按代码块分段，检查每个代码块前后是否有"孤立"段落
# （代码块与最近说明文字之间夹了超过 N 行空白行/其他段落）
# 简化实现：检查代码块前/后间隔行数是否过大（>15行视为远离）

MAX_GAP = 15  # 超过此行数视为代码与说明不邻近


def check_4_1(text: str) -> dict:
    lines = text.splitlines()
    code_block_line_ranges = []

    in_code = False
    start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("```"):
            if not in_code:
                in_code = True
                start = i
            else:
                in_code = False
                code_block_line_ranges.append((start, i))

    if not code_block_line_ranges:
        return {"score": 3, "max": 3, "note": "无代码块，不适用"}

    isolated = 0
    for cb_start, cb_end in code_block_line_ranges:
        # 检查代码块前的间距
        gap_before = 0
        for j in range(cb_start - 1, max(-1, cb_start - MAX_GAP - 5), -1):
            if lines[j].strip():
                break
            gap_before += 1

        # 检查代码块后的间距
        gap_after = 0
        for j in range(cb_end + 1, min(len(lines), cb_end + MAX_GAP + 5)):
            if lines[j].strip():
                break
            gap_after += 1

        if gap_before > MAX_GAP and gap_after > MAX_GAP:
            isolated += 1

    n_blocks = len(code_block_line_ranges)
    isolation_rate = isolated / n_blocks if n_blocks else 0

    if isolation_rate == 0:
        score, note = 3, f"所有{n_blocks}个代码块均与说明邻近"
    elif isolation_rate < 0.2:
        score, note = 2, f"{isolated}/{n_blocks}个代码块与说明距离较远"
    elif isolation_rate < 0.4:
        score, note = 1, f"{isolated}/{n_blocks}个代码块孤立"
    else:
        score, note = 0, f"{isolated}/{n_blocks}个代码块严重孤立"

    return {"score": score, "max": 3,
            "total_code_blocks": n_blocks,
            "isolated_blocks": isolated, "note": note}


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def run_checks(tutorial_path: str) -> dict:
    text = load_markdown(tutorial_path)
    prose, code_blocks = split_code_blocks(text)

    results = {
        "input_file": tutorial_path,
        "indicators": {
            "1.3": check_1_3(prose, code_blocks),
            "1.4": check_1_4(text),
            "1.5": check_1_5(text),
            "2.4": check_2_4(text),
            "2.5": check_2_5(text),
            "3.2": check_3_2(text),
            "4.1": check_4_1(text),
        }
    }

    total_score = sum(v["score"] for v in results["indicators"].values())
    total_max   = sum(v["max"]   for v in results["indicators"].values())
    results["total_score"] = total_score
    results["total_max"]   = total_max

    return results


def main():
    parser = argparse.ArgumentParser(description="judgeCLAUDE 静态指标检查")
    parser.add_argument("--input",  required=True, help="教程 Markdown 文件路径")
    parser.add_argument("--output", required=False, help="输出 JSON 路径（省略则打印到 stdout）")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"ERROR: 文件不存在：{args.input}", file=sys.stderr)
        sys.exit(1)

    results = run_checks(args.input)

    output_str = json.dumps(results, ensure_ascii=False, indent=2)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(output_str)
        print(f"[OK] 静态检查完成，结果已写入：{args.output}")
        print(f"   总分：{results['total_score']} / {results['total_max']}")
        for ind_id, v in results["indicators"].items():
            note = v['note']
            print(f"   [{ind_id}] {v['score']}/{v['max']} | {note}")
    else:
        print(output_str)


if __name__ == "__main__":
    main()
