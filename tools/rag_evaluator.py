#!/usr/bin/env python3
"""
rag_evaluator.py — RAG系统评估工具

评估RAG检索质量，专注于：
1. 工具归属错误（tool attribution error）检测
2. 精确率/召回率/忠实度/相关性指标
3. A/B测试前后对比
4. 回归测试套件

用法：
  # 运行完整评估套件
  python tools/rag_evaluator.py --mode eval

  # 只运行归属错误测试
  python tools/rag_evaluator.py --mode attribution

  # A/B对比（需要提前保存baseline结果）
  python tools/rag_evaluator.py --mode baseline --output results/baseline.json
  python tools/rag_evaluator.py --mode compare --baseline results/baseline.json

  # 回归测试
  python tools/rag_evaluator.py --mode regression --threshold 0.85
"""

import os
import sys
import json
import time
import hashlib
import argparse
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, field, asdict
from pathlib import Path

import chromadb
import dashscope
from dashscope import TextEmbedding
import yaml

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")

# ---------------------------------------------------------------------------
# API Key 初始化
# ---------------------------------------------------------------------------

def setup_api():
    try:
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
            api_key = (config.get("llm") or {}).get("api_key") or config.get("api_key")
            if api_key:
                dashscope.api_key = api_key
    except Exception:
        pass
    if not dashscope.api_key:
        dashscope.api_key = os.getenv("DASHSCOPE_API_KEY", "")


def get_embedding(text: str) -> Optional[List[float]]:
    try:
        resp = TextEmbedding.call(
            model=TextEmbedding.Models.text_embedding_v3,
            input=[text]
        )
        if resp.status_code == 200:
            return resp.output["embeddings"][0]["embedding"]
    except Exception as e:
        print(f"  [警告] embedding失败: {e}", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class TestCase:
    """单个测试用例"""
    id: str
    query: str
    # 期望结果
    expected_sources: List[str]        # 必须出现的source文件名（部分匹配）
    forbidden_sources: List[str]       # 不应出现的source文件名（归属错误检测）
    expected_keywords: List[str]       # 期望在检索内容中出现的关键词
    forbidden_keywords: List[str]      # 不应在高排名结果中出现的关键词
    top_k: int = 5
    category: str = "general"
    description: str = ""


@dataclass
class RetrievalResult:
    """单次检索结果"""
    query: str
    documents: List[str]
    sources: List[str]
    section_paths: List[str]
    distances: List[float]


@dataclass
class CaseResult:
    """单个测试用例的评估结果"""
    case_id: str
    query: str
    category: str
    # 精确率：期望source在前top_k中出现的比例
    precision: float
    # 召回率：期望source有多少被检索到
    recall: float
    # 归属错误：forbidden_source是否出现在前3名
    attribution_error: bool
    attribution_error_details: List[str]
    # 关键词命中
    keyword_hit_rate: float
    forbidden_keyword_in_top3: bool
    # 原始结果
    retrieved_sources: List[str]
    retrieved_section_paths: List[str]
    pass_: bool  # 综合是否通过


@dataclass
class EvalReport:
    """完整评估报告"""
    timestamp: str
    total_cases: int
    passed: int
    failed: int
    # 聚合指标
    mean_precision: float
    mean_recall: float
    mean_keyword_hit_rate: float
    attribution_error_rate: float  # 有归属错误的用例比例
    # 按类别分组
    by_category: Dict[str, Dict]
    # 明细
    case_results: List[Dict]
    # 整体通过率
    pass_rate: float


# ---------------------------------------------------------------------------
# 测试用例库
# ---------------------------------------------------------------------------

# 核心测试用例：覆盖已知失败模式和常见查询场景
GOLDEN_TEST_CASES: List[TestCase] = [

    # ============================================================
    # 归属错误测试组（Tool Attribution Error）
    # 这是最重要的一组测试
    # 失败模式：chunk来自"ABACUS+Atomkit 计算态密度和能带.md"，
    # 但正文内容讨论abacus-plot；agent错误地将abacus-plot的工作流
    # 归属到Atomkit
    # ============================================================

    TestCase(
        id="attr_001",
        query="使用abacus-plot绘制PDOS的步骤",
        expected_sources=["ABACUS 计算 PDOS.md"],
        forbidden_sources=["ABACUS+Atomkit 计算态密度和能带.md"],
        expected_keywords=["abacus-plot", "config.json", "plot-tools"],
        forbidden_keywords=["atomkit", "Atomkit"],
        top_k=5,
        category="attribution_error",
        description="abacus-plot PDOS绘制工作流不应归属到Atomkit文档"
    ),

    TestCase(
        id="attr_002",
        query="abacus-plot安装和配置方法",
        expected_sources=["ABACUS 计算 PDOS.md"],
        forbidden_sources=["ABACUS+Atomkit 计算态密度和能带.md"],
        expected_keywords=["abacus-plot", "setup.py", "plot-tools"],
        forbidden_keywords=["atomkit", "Atomkit"],
        top_k=5,
        category="attribution_error",
        description="abacus-plot安装不应归属到Atomkit"
    ),

    TestCase(
        id="attr_003",
        query="PDOS绘制工具 config.json pdosfile efermi species",
        expected_sources=["ABACUS 计算 PDOS.md"],
        forbidden_sources=["ABACUS+Atomkit 计算态密度和能带.md"],
        expected_keywords=["pdosfile", "efermi", "species"],
        forbidden_keywords=["Atomkit", "atomkit"],
        top_k=5,
        category="attribution_error",
        description="abacus-plot的config.json参数不应归属到Atomkit"
    ),

    TestCase(
        id="attr_004",
        query="能带态密度计算 Atomkit",
        # 正向测试：查Atomkit时确实应该找到Atomkit文档
        expected_sources=["ABACUS+Atomkit 计算态密度和能带.md"],
        forbidden_sources=[],
        expected_keywords=["Atomkit", "atomkit"],
        forbidden_keywords=[],
        top_k=5,
        category="attribution_error",
        description="正向：查询Atomkit时应找到Atomkit文档（基准对照）"
    ),

    TestCase(
        id="attr_005",
        query="abacus-plot -d -o 命令",
        expected_sources=["ABACUS 计算 PDOS.md"],
        forbidden_sources=["ABACUS+Atomkit 计算态密度和能带.md"],
        expected_keywords=["abacus-plot"],
        forbidden_keywords=[],
        top_k=3,
        category="attribution_error",
        description="abacus-plot命令行参数查询"
    ),

    # ============================================================
    # 精确率测试组：标准参数查询
    # ============================================================

    TestCase(
        id="prec_001",
        query="ABACUS INPUT文件 calculation scf relax",
        expected_sources=["ABACUS 使用教程｜电子自洽迭代.md",
                          "ABACUS 使用教程｜结构优化.md",
                          "abacus_user_guide.md"],
        forbidden_sources=[],
        expected_keywords=["calculation", "scf", "INPUT"],
        forbidden_keywords=[],
        top_k=5,
        category="precision",
        description="基本INPUT参数查询"
    ),

    TestCase(
        id="prec_002",
        query="DFT+U 强关联体系 NiO Hubbard U参数设置",
        expected_sources=["ABACUS 使用教程｜DFT+U 计算.md"],
        forbidden_sources=[],
        expected_keywords=["DFT+U", "dft_plus_u", "Hubbard"],
        forbidden_keywords=[],
        top_k=5,
        category="precision",
        description="DFT+U参数检索"
    ),

    TestCase(
        id="prec_003",
        query="弹性常数计算 elastic ABACUS",
        expected_sources=["ABACUS 使用教程｜结构优化.md"],
        forbidden_sources=[],
        expected_keywords=["弹性", "elastic", "stress"],
        forbidden_keywords=[],
        top_k=5,
        category="precision",
        description="弹性常数相关查询"
    ),

    TestCase(
        id="prec_004",
        query="ABACUS 磁性计算 nspin=2 初始磁矩设置",
        expected_sources=["ABACUS 使用教程｜磁性材料计算.md"],
        forbidden_sources=[],
        expected_keywords=["nspin", "磁矩", "磁性"],
        forbidden_keywords=[],
        top_k=5,
        category="precision",
        description="磁性计算参数查询"
    ),

    TestCase(
        id="prec_005",
        query="ABACUS 随机波函数DFT SDFT esolver_type",
        expected_sources=["ABACUS 随机波函数 DFT 方法使用教程.md"],
        forbidden_sources=[],
        expected_keywords=["sdft", "esolver_type", "随机波函数"],
        forbidden_keywords=[],
        top_k=5,
        category="precision",
        description="SDFT方法参数查询"
    ),

    # ============================================================
    # 忠实度/内容质量测试组
    # 确保检索到的内容包含实际可用的信息，而非仅标题匹配
    # ============================================================

    TestCase(
        id="faith_001",
        query="LCAO基组轨道文件格式 .orb文件命名规则",
        expected_sources=[],  # 任何含轨道文件说明的文档都可以
        forbidden_sources=[],
        expected_keywords=["orb", "au", "Ry"],  # 轨道文件名必须含这些单位
        forbidden_keywords=[],
        top_k=5,
        category="faithfulness",
        description="轨道文件命名规则查询，期望内容包含实际文件名格式"
    ),

    TestCase(
        id="faith_002",
        query="能带计算 高对称K点路径 KLINE_MODE",
        expected_sources=[],
        forbidden_sources=[],
        expected_keywords=["KLINE_MODE", "K点", "高对称"],
        forbidden_keywords=[],
        top_k=5,
        category="faithfulness",
        description="能带K点路径查询"
    ),

    TestCase(
        id="faith_003",
        query="ABACUS隐式溶剂模型 imp_sol epsilon参数",
        expected_sources=["ABACUS 隐式溶剂模型使用教程.md"],
        forbidden_sources=[],
        expected_keywords=["imp_sol", "epsilon", "溶剂"],
        forbidden_keywords=[],
        top_k=5,
        category="faithfulness",
        description="隐式溶剂参数查询"
    ),

    # ============================================================
    # 多源融合测试：一个查询应触发多个相关文档
    # ============================================================

    TestCase(
        id="multi_001",
        query="ABACUS SCF收敛 scf_thr mixing_type 自洽迭代",
        expected_sources=["ABACUS 使用教程｜电子自洽迭代.md",
                          "abacus_user_guide.md"],
        forbidden_sources=[],
        expected_keywords=["scf_thr", "mixing"],
        forbidden_keywords=[],
        top_k=8,
        category="multi_source",
        description="SCF参数应触发多个文档"
    ),

    TestCase(
        id="multi_002",
        query="态密度DOS能带band ABACUS计算流程",
        expected_sources=["ABACUS+Atomkit 计算态密度和能带.md",
                          "ABACUS 计算 PDOS.md"],
        forbidden_sources=[],
        expected_keywords=["DOS", "band", "nscf"],
        forbidden_keywords=[],
        top_k=8,
        category="multi_source",
        description="DOS+能带查询应融合多个文档"
    ),
]


# ---------------------------------------------------------------------------
# 检索执行
# ---------------------------------------------------------------------------

class RAGEvaluator:
    def __init__(self):
        setup_api()
        self.client = chromadb.PersistentClient(path=DB_PATH)
        self.collection = self.client.get_collection(name="abacus_knowledge")

    def retrieve(self, query: str, top_k: int = 5) -> RetrievalResult:
        """执行检索，返回结构化结果"""
        vec = get_embedding(query)
        if vec is None:
            return RetrievalResult(query=query, documents=[], sources=[],
                                   section_paths=[], distances=[])
        try:
            results = self.collection.query(
                query_embeddings=[vec],
                n_results=top_k,
                include=["documents", "metadatas", "distances"]
            )
        except Exception as e:
            print(f"  [错误] 检索失败: {e}", file=sys.stderr)
            return RetrievalResult(query=query, documents=[], sources=[],
                                   section_paths=[], distances=[])

        docs = results["documents"][0]
        metas = results["metadatas"][0]
        dists = results["distances"][0]

        sources = [m.get("source", "") for m in metas]
        section_paths = [m.get("section_path", "") for m in metas]

        return RetrievalResult(
            query=query,
            documents=docs,
            sources=sources,
            section_paths=section_paths,
            distances=dists
        )

    def evaluate_case(self, case: TestCase) -> CaseResult:
        """评估单个测试用例"""
        result = self.retrieve(case.query, case.top_k)

        retrieved_sources = result.sources
        retrieved_docs = result.documents

        # --- 精确率：前top_k中，expected_source至少部分命中的数量 ---
        def source_match(retrieved: str, expected: str) -> bool:
            return expected.lower() in retrieved.lower()

        matched_expected = set()
        for rs in retrieved_sources:
            for es in case.expected_sources:
                if source_match(rs, es):
                    matched_expected.add(es)

        if case.expected_sources:
            precision = len(matched_expected) / len(case.expected_sources)
            recall = len(matched_expected) / len(case.expected_sources)
        else:
            precision = 1.0  # 无期望source时不计算
            recall = 1.0

        # --- 归属错误：forbidden_source是否出现在前3名 ---
        top3_sources = retrieved_sources[:3]
        attribution_errors = []
        for fs in case.forbidden_sources:
            for rs in top3_sources:
                if source_match(rs, fs):
                    attribution_errors.append(
                        f"forbidden '{fs}' 出现在前3: '{rs}'"
                    )

        attribution_error = len(attribution_errors) > 0

        # --- 关键词命中率：期望关键词在检索内容中出现 ---
        all_content = " ".join(retrieved_docs).lower()
        if case.expected_keywords:
            hits = sum(1 for kw in case.expected_keywords if kw.lower() in all_content)
            keyword_hit_rate = hits / len(case.expected_keywords)
        else:
            keyword_hit_rate = 1.0

        # --- forbidden关键词检查（前3名内容）---
        top3_content = " ".join(retrieved_docs[:3]).lower()
        forbidden_keyword_in_top3 = any(
            kw.lower() in top3_content for kw in case.forbidden_keywords
        )

        # --- 综合pass判断 ---
        # 规则：
        #   1. 无归属错误（最严格）
        #   2. precision >= 0.5（至少一半期望source被找到）
        #   3. keyword_hit_rate >= 0.5
        #   4. forbidden_keyword_in_top3 不影响pass，但会记录
        pass_ = (
            not attribution_error
            and precision >= 0.5
            and keyword_hit_rate >= 0.5
        )

        return CaseResult(
            case_id=case.id,
            query=case.query,
            category=case.category,
            precision=precision,
            recall=recall,
            attribution_error=attribution_error,
            attribution_error_details=attribution_errors,
            keyword_hit_rate=keyword_hit_rate,
            forbidden_keyword_in_top3=forbidden_keyword_in_top3,
            retrieved_sources=retrieved_sources,
            retrieved_section_paths=result.section_paths,
            pass_=pass_
        )

    def run_eval(self, cases: List[TestCase] = None) -> EvalReport:
        """运行完整评估套件"""
        if cases is None:
            cases = GOLDEN_TEST_CASES

        print(f"\n{'='*60}")
        print(f"RAG 评估套件 — {len(cases)} 个测试用例")
        print(f"{'='*60}\n")

        case_results = []
        for case in cases:
            print(f"  [{case.id}] {case.description[:50]}...", end=" ", flush=True)
            cr = self.evaluate_case(case)
            case_results.append(cr)

            status = "✓ PASS" if cr.pass_ else "✗ FAIL"
            print(f"{status} | P={cr.precision:.2f} R={cr.recall:.2f} "
                  f"KW={cr.keyword_hit_rate:.2f}"
                  + (" [ATTR_ERR!]" if cr.attribution_error else ""))

            if cr.attribution_error_details:
                for d in cr.attribution_error_details:
                    print(f"      [归属错误] {d}")

            time.sleep(0.3)  # 避免API限速

        # 计算聚合指标
        total = len(case_results)
        passed = sum(1 for r in case_results if r.pass_)
        failed = total - passed

        mean_precision = sum(r.precision for r in case_results) / total
        mean_recall = sum(r.recall for r in case_results) / total
        mean_kw = sum(r.keyword_hit_rate for r in case_results) / total
        attr_err_rate = sum(1 for r in case_results if r.attribution_error) / total

        # 按类别分组
        by_category: Dict[str, Dict] = {}
        for cr in case_results:
            cat = cr.category
            if cat not in by_category:
                by_category[cat] = {"total": 0, "passed": 0,
                                     "attr_errors": 0, "precision": [],
                                     "recall": []}
            by_category[cat]["total"] += 1
            if cr.pass_:
                by_category[cat]["passed"] += 1
            if cr.attribution_error:
                by_category[cat]["attr_errors"] += 1
            by_category[cat]["precision"].append(cr.precision)
            by_category[cat]["recall"].append(cr.recall)

        # 计算每类别平均
        for cat, stats in by_category.items():
            stats["pass_rate"] = stats["passed"] / stats["total"]
            stats["mean_precision"] = sum(stats["precision"]) / len(stats["precision"])
            stats["mean_recall"] = sum(stats["recall"]) / len(stats["recall"])
            del stats["precision"]
            del stats["recall"]

        report = EvalReport(
            timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
            total_cases=total,
            passed=passed,
            failed=failed,
            mean_precision=mean_precision,
            mean_recall=mean_recall,
            mean_keyword_hit_rate=mean_kw,
            attribution_error_rate=attr_err_rate,
            by_category=by_category,
            case_results=[asdict(r) for r in case_results],
            pass_rate=passed / total
        )

        self._print_summary(report)
        return report

    def _print_summary(self, report: EvalReport):
        print(f"\n{'='*60}")
        print(f"评估摘要 ({report.timestamp})")
        print(f"{'='*60}")
        print(f"  总用例: {report.total_cases}  通过: {report.passed}  失败: {report.failed}")
        print(f"  通过率:           {report.pass_rate:.1%}")
        print(f"  平均精确率:       {report.mean_precision:.3f}")
        print(f"  平均召回率:       {report.mean_recall:.3f}")
        print(f"  平均关键词命中率: {report.mean_keyword_hit_rate:.3f}")
        print(f"  归属错误率:       {report.attribution_error_rate:.1%}  "
              f"({'⚠ 存在归属错误' if report.attribution_error_rate > 0 else '✓ 无归属错误'})")
        print()
        print("  按类别：")
        for cat, stats in report.by_category.items():
            print(f"    [{cat:20s}] pass={stats['pass_rate']:.1%}  "
                  f"P={stats['mean_precision']:.2f}  R={stats['mean_recall']:.2f}  "
                  f"attr_err={stats['attr_errors']}/{stats['total']}")
        print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# A/B 对比分析
# ---------------------------------------------------------------------------

def compare_reports(baseline_path: str, current: EvalReport) -> Dict:
    """对比基线和当前报告，计算改进幅度"""
    with open(baseline_path, "r", encoding="utf-8") as f:
        baseline = json.load(f)

    delta = {
        "pass_rate":              current.pass_rate - baseline["pass_rate"],
        "mean_precision":         current.mean_precision - baseline["mean_precision"],
        "mean_recall":            current.mean_recall - baseline["mean_recall"],
        "mean_keyword_hit_rate":  current.mean_keyword_hit_rate - baseline["mean_keyword_hit_rate"],
        "attribution_error_rate": current.attribution_error_rate - baseline["attribution_error_rate"],
    }

    # 逐用例对比（找到哪些用例变好/变差）
    baseline_by_id = {r["case_id"]: r for r in baseline["case_results"]}
    regressions = []
    improvements = []
    for cr in current.case_results:
        cid = cr["case_id"]
        if cid not in baseline_by_id:
            continue
        b = baseline_by_id[cid]
        if b["pass_"] and not cr["pass_"]:
            regressions.append(cid)
        elif not b["pass_"] and cr["pass_"]:
            improvements.append(cid)

    print("\n" + "="*60)
    print("A/B 对比报告")
    print("="*60)
    print(f"  基线时间: {baseline['timestamp']}")
    print(f"  当前时间: {current.timestamp}")
    print()
    print("  指标变化（当前 - 基线）:")
    for k, v in delta.items():
        arrow = "↑" if v > 0 else ("↓" if v < 0 else "→")
        # 归属错误率降低是好事
        good = (v > 0) if k != "attribution_error_rate" else (v < 0)
        mark = "✓" if good else ("✗" if (v < 0 and k != "attribution_error_rate") or
                                          (v > 0 and k == "attribution_error_rate") else "→")
        print(f"    {k:30s}: {arrow} {v:+.4f}  {mark}")

    if regressions:
        print(f"\n  ⚠ 退步用例 ({len(regressions)})：{regressions}")
    if improvements:
        print(f"  ✓ 改进用例 ({len(improvements)})：{improvements}")
    if not regressions and not improvements:
        print("\n  → 无用例级别变化")

    print("="*60 + "\n")

    return {
        "delta": delta,
        "regressions": regressions,
        "improvements": improvements,
        "baseline_timestamp": baseline["timestamp"],
        "current_timestamp": current.timestamp,
    }


# ---------------------------------------------------------------------------
# 回归测试（CI集成）
# ---------------------------------------------------------------------------

def run_regression(threshold: float = 0.85) -> bool:
    """
    回归测试：检查关键指标是否达到阈值
    返回True=通过，False=失败
    适合集成到CI/CD流程
    """
    evaluator = RAGEvaluator()

    # 只运行归属错误测试组（速度快，最关键）
    attribution_cases = [c for c in GOLDEN_TEST_CASES if c.category == "attribution_error"]
    report = evaluator.run_eval(attribution_cases)

    print(f"回归测试结果:")
    print(f"  通过率: {report.pass_rate:.1%}  (阈值: {threshold:.1%})")
    print(f"  归属错误率: {report.attribution_error_rate:.1%}  (目标: 0%)")

    ok = report.pass_rate >= threshold and report.attribution_error_rate == 0.0

    if ok:
        print("  [PASS] 回归测试通过")
    else:
        print("  [FAIL] 回归测试失败")
        if report.attribution_error_rate > 0:
            print("  ⚠ 存在工具归属错误，RAG改动可能引入新问题！")

    return ok


# ---------------------------------------------------------------------------
# 辅助：检查单个查询的归属情况（调试用）
# ---------------------------------------------------------------------------

def debug_attribution(query: str, top_k: int = 5):
    """调试单个查询的来源归属，显示详细的section_path信息"""
    evaluator = RAGEvaluator()
    result = evaluator.retrieve(query, top_k)

    print(f"\n查询: {query}")
    print(f"{'='*70}")
    for i, (doc, src, sp, dist) in enumerate(zip(
            result.documents, result.sources,
            result.section_paths, result.distances), 1):
        print(f"\n[{i}] 来源: {src}")
        print(f"    Section: {sp}")
        print(f"    距离: {dist:.4f}")
        print(f"    内容摘要: {doc[:150].replace(chr(10), ' ')}...")


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="RAG评估工具")
    parser.add_argument(
        "--mode",
        choices=["eval", "attribution", "baseline", "compare", "regression", "debug"],
        default="eval",
        help="运行模式"
    )
    parser.add_argument("--output", help="输出JSON路径")
    parser.add_argument("--baseline", help="基线JSON路径（用于compare模式）")
    parser.add_argument("--threshold", type=float, default=0.85,
                        help="回归测试通过率阈值（默认0.85）")
    parser.add_argument("--query", help="调试模式下的查询词")
    parser.add_argument("--top-k", type=int, default=5)

    args = parser.parse_args()

    if not dashscope.api_key:
        setup_api()
    if not dashscope.api_key:
        print("错误: 未配置API Key", file=sys.stderr)
        sys.exit(1)

    if args.mode == "debug":
        query = args.query or "使用abacus-plot绘制PDOS的步骤"
        debug_attribution(query, args.top_k)
        return

    if args.mode == "regression":
        ok = run_regression(args.threshold)
        sys.exit(0 if ok else 1)

    evaluator = RAGEvaluator()

    if args.mode == "attribution":
        cases = [c for c in GOLDEN_TEST_CASES if c.category == "attribution_error"]
        report = evaluator.run_eval(cases)
    elif args.mode in ("eval", "baseline"):
        report = evaluator.run_eval()
    else:
        report = None

    if args.mode == "compare":
        if not args.baseline:
            print("错误: compare模式需要 --baseline 参数", file=sys.stderr)
            sys.exit(1)
        report = evaluator.run_eval()
        compare_reports(args.baseline, report)

    if report and args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(asdict(report), f, ensure_ascii=False, indent=2)
        print(f"[OK] 评估结果已保存: {args.output}")
    elif report and args.mode == "baseline":
        out = args.output or "docs/dataBD/rag_eval_baseline.json"
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(asdict(report), f, ensure_ascii=False, indent=2)
        print(f"[OK] 基线已保存: {out}")


if __name__ == "__main__":
    main()
