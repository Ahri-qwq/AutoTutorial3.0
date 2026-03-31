import os
import glob
import uuid
import time
import shutil
import hashlib
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Iterable

import chromadb

import dashscope
from dashscope import TextEmbedding
import yaml

try:
    from docx import Document
except ImportError:
    Document = None

try:
    from pypdf import PdfReader
except ImportError:
    PdfReader = None

# ================= 路径配置 =================

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)

CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")
SOURCE_DIR = os.path.join(ROOT_DIR, "data", "knowledge_source")  # 重建来源&归档目标
ADD_DIR = os.path.join(ROOT_DIR, "data", "knowledge_add")        # 增量来源
DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
COLLECTION_NAME = "abacus_knowledge"

# 切分参数（字符级的二次切分上限；语义切分会优先保持段落/标题/代码块）
CHUNK_SIZE = 3000
OVERLAP = 200

# 支持的文件类型
SUPPORTED_EXTS = [".md", ".txt", ".docx", ".pdf"]

# 入库严格模式：embedding 失败时直接跳过该文件（不写 0 向量污染索引）
STRICT_EMBEDDING = os.getenv("AUTOTUTORIAL_STRICT_EMBEDDING", "1").strip() not in ("0", "false", "False")

# ================= API Key =================

def setup_api_key() -> None:
    api_key = ""
    try:
        if os.path.exists(CONFIG_PATH):
            with open(CONFIG_PATH, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f) or {}
            if config.get("llm"):
                api_key = (config["llm"] or {}).get("api_key", "")
            else:
                api_key = config.get("api_key", "")
    except Exception:
        pass

    if not api_key:
        api_key = os.getenv("DASHSCOPE_API_KEY", "")

    if api_key:
        dashscope.api_key = api_key
    else:
        print("❌ 警告: 未找到 API Key，后续调用将失败。")

# ================= Embedding =================

@dataclass
class EmbeddingResult:
    embeddings: List[List[float]]
    ok: bool
    error: str = ""


def embed_texts(texts: List[str], batch_size: int = 5, sleep_s: float = 0.2) -> EmbeddingResult:
    """显式调用 embedding API，便于在失败时选择跳过而不是写 0 向量。"""
    all_embeddings: List[List[float]] = []
    for i in range(0, len(texts), batch_size):
        batch = texts[i: i + batch_size]
        try:
            resp = TextEmbedding.call(
                model=TextEmbedding.Models.text_embedding_v3,
                input=batch
            )
            if resp.status_code == 200:
                embeddings = [item["embedding"] for item in resp.output["embeddings"]]
                all_embeddings.extend(embeddings)
            else:
                return EmbeddingResult([], ok=False, error=f"API Error: {resp.message}")
        except Exception as e:
            return EmbeddingResult([], ok=False, error=f"Network Exception: {e}")
        time.sleep(sleep_s)
    return EmbeddingResult(all_embeddings, ok=True)

# ================= 读取与切分 =================

def read_file(filepath: str) -> str:
    ext = os.path.splitext(filepath)[1].lower()
    content = ""

    try:
        if ext in [".md", ".txt"]:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

        elif ext == ".docx":
            if Document is None:
                print(f"⚠️ 跳过 DOCX（未安装 python-docx）: {filepath}")
                return ""
            doc = Document(filepath)
            content = "\n".join([p.text for p in doc.paragraphs])

        elif ext == ".pdf":
            if PdfReader is None:
                print(f"⚠️ 跳过 PDF（未安装 pypdf）: {filepath}")
                return ""
            reader = PdfReader(filepath)
            pages = []
            for page in reader.pages:
                pages.append(page.extract_text() or "")
            content = "\n".join(pages)

        else:
            print(f"⚠️ 不支持的文件类型 {ext}: {filepath}")
            return ""

    except Exception as e:
        print(f"⚠️ 无法读取 {filepath}: {e}")
        return ""

    return content


def _hash_text(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8", errors="ignore")).hexdigest()


def _stable_chunk_id(doc_relpath: str, section_path: str, chunk_index: int, chunk_text: str) -> str:
    """稳定可复用的 chunk id（用于证据回链与断言审计）。"""
    base = f"{doc_relpath}||{section_path}||{chunk_index}||{_hash_text(chunk_text)}"
    return hashlib.sha1(base.encode("utf-8", errors="ignore")).hexdigest()[:24]


def _guess_doc_type(doc_relpath: str) -> str:
    p = doc_relpath.lower()
    if any(k in p for k in ["manual", "docs", "document", "official", "abacus" + os.sep + "doc"]):
        return "official_doc"
    if any(k in p for k in ["example", "examples", "demo", "sample"]):
        return "example"
    if any(k in p for k in ["issue", "faq", "troubleshoot", "error"]):
        return "faq_issue"
    if any(k in p for k in ["paper", "article", "publication"]):
        return "paper"
    return "note"


@dataclass
class Chunk:
    text: str
    section_path: str
    span_start: int
    span_end: int


def _split_by_markdown_sections(text: str) -> List[Chunk]:
    """按 Markdown 标题切分，尽量保持语义完整；同时保留 code block 原子性。"""
    import re
    lines = text.splitlines(keepends=True)
    chunks: List[Chunk] = []

    in_code = False
    header_stack: List[str] = []  # 记录 #/##/### 的路径

    buf: List[str] = []
    buf_start = 0
    pos = 0

    def flush(end_pos: int):
        nonlocal buf, buf_start
        if not buf:
            return
        content = "".join(buf).strip()
        if content:
            section_path = " > ".join(header_stack) if header_stack else "(no_heading)"
            # 将完整标题路径写入 chunk 文本开头（参与 embedding，提供语境）
            # 使用原始文档标题（包含 H1），不做过滤
            # Agent 结合正文内容综合判断，不被元数据强制锚定
            if header_stack:
                header_prefix = " > ".join(header_stack)
                content = f"[{header_prefix}]\n{content}"
            chunks.append(Chunk(text=content, section_path=section_path, span_start=buf_start, span_end=end_pos))
        buf = []

    for ln in lines:
        # code fence toggle
        if ln.strip().startswith("```"):
            in_code = not in_code
            buf.append(ln)
            pos += len(ln)
            continue

        # heading line (only when not in code)
        m = None
        if not in_code:
            m = re.match(r"^(#{1,6})\s+(.*)\s*$", ln.strip())

        if m:
            flush(pos)
            buf_start = pos

            level = len(m.group(1))
            title = m.group(2).strip()
            if level <= len(header_stack):
                header_stack[:] = header_stack[: level - 1]
            header_stack.append(title)

            buf.append(ln)
            pos += len(ln)
            continue

        buf.append(ln)
        pos += len(ln)

    flush(pos)
    return chunks


def _split_long_chunk_preserve_overlap(chunk: Chunk, max_len: int = CHUNK_SIZE, overlap: int = OVERLAP) -> List[Chunk]:
    """对过长块进行字符级二次切分，但保持 section_path 与 span 映射。"""
    t = chunk.text
    if len(t) <= max_len:
        return [chunk]

    out: List[Chunk] = []
    start = 0
    while start < len(t):
        end = min(start + max_len, len(t))
        piece = t[start:end]
        out.append(Chunk(
            text=piece,
            section_path=chunk.section_path,
            span_start=chunk.span_start + start,
            span_end=chunk.span_start + end
        ))
        if end == len(t):
            break
        start += (max_len - overlap)
    return out


def split_document(filepath: str, text: str) -> List[Chunk]:
    ext = os.path.splitext(filepath)[1].lower()

    if not text or not text.strip():
        return []

    # 检查 Markdown frontmatter 中的 chunk_strategy
    if ext == ".md" and text.startswith("---"):
        lines = text.split('\n', 20)  # 只检查前20行
        for line in lines[1:]:
            if line.strip() == "---":
                break
            if "chunk_strategy:" in line and "whole" in line:
                # 整个文档作为一个 chunk
                return [Chunk(text=text, section_path="", span_start=0, span_end=len(text))]

    if ext == ".md":
        base = _split_by_markdown_sections(text)
        refined: List[Chunk] = []
        for c in base:
            refined.extend(_split_long_chunk_preserve_overlap(c))
        return refined

    if ext == ".txt":
        paras = []
        spans = []
        pos = 0
        buf = []
        buf_start = 0
        for ln in text.splitlines(keepends=True):
            if ln.strip() == "":
                if buf:
                    para = "".join(buf).strip()
                    if para:
                        paras.append(para)
                        spans.append((buf_start, pos))
                    buf = []
                buf_start = pos + len(ln)
            else:
                if not buf:
                    buf_start = pos
                buf.append(ln)
            pos += len(ln)
        if buf:
            para = "".join(buf).strip()
            if para:
                paras.append(para)
                spans.append((buf_start, pos))

        refined: List[Chunk] = []
        for i, para in enumerate(paras):
            c = Chunk(text=para, section_path=f"para_{i}", span_start=spans[i], span_end=spans[i])[1]
            refined.extend(_split_long_chunk_preserve_overlap(c))
        return refined

    c = Chunk(text=text.strip(), section_path="(fulltext)", span_start=0, span_end=len(text))
    return _split_long_chunk_preserve_overlap(c)

# ================= 扫描与入库通用函数 =================

def scan_files(root_dir: str) -> List[str]:
    if not os.path.exists(root_dir):
        return []
    patterns = [os.path.join(root_dir, f"**/*{ext}") for ext in SUPPORTED_EXTS]
    files: List[str] = []
    for p in patterns:
        files.extend(glob.glob(p, recursive=True))
    return files


def ingest_files(
    collection,
    files: List[str],
    ingest_root_dir: str,
    ingest_type: Optional[str] = None
) -> Tuple[int, List[str]]:
    """返回：(新增chunk数, 成功写入过至少一个chunk的文件路径列表)。"""

    total_chunks = 0
    ingested_files: List[str] = []

    for idx, filepath in enumerate(files):
        fname = os.path.basename(filepath)
        relpath = os.path.relpath(filepath, ingest_root_dir).replace("\\", "/")
        ext = os.path.splitext(filepath)[1].lower()


        print(f"[{idx+1}/{len(files)}] 处理: {relpath}...", end="", flush=True)

        content = read_file(filepath)
        if not content.strip():
            print(" [空]")
            continue

        doc_type = _guess_doc_type(relpath)
        chunks = split_document(filepath, content)
        if not chunks:
            print(" [无可切分内容]")
            continue

        docs: List[str] = []
        ids: List[str] = []
        metadatas: List[Dict] = []

        ingest_time = int(time.time())

        for i, c in enumerate(chunks):
            text = c.text.strip()
            if not text:
                continue

            chunk_text_hash = _hash_text(text)
            chunk_id = _stable_chunk_id(relpath, c.section_path, i, text)

            docs.append(text)
            ids.append(chunk_id)
            md: Dict = {
                "source": fname,
                "doc_path": relpath,
                "doc_type": doc_type,
                "ext": ext,
                "chunk_index": i,
                "section_path": c.section_path,
                "span_start": c.span_start,
                "span_end": c.span_end,
                "text_hash": chunk_text_hash,
                "ingest_time": ingest_time,
            }
            if ingest_type:
                md["ingest_type"] = ingest_type
            metadatas.append(md)

        if not docs:
            print(" [无有效chunk]")
            continue

        emb = embed_texts(docs)
        if not emb.ok:
            msg = f" ❌ embedding失败: {emb.error}"
            if STRICT_EMBEDDING:
                print(msg + "（严格模式：跳过该文件）")
                continue
            else:
                print(msg + "（非严格模式：跳过该文件）")
                continue

        try:
            collection.add(documents=docs, metadatas=metadatas, ids=ids, embeddings=emb.embeddings)
            print(f" ✅ {len(docs)} 片段")
            total_chunks += len(docs)
            ingested_files.append(filepath)
        except Exception as e:
            print(f" ❌ 写入失败: {e}")

    return total_chunks, ingested_files

# ================= 两种模式 =================

def rebuild_db() -> None:
    print("==========================================")
    print(" 🚀 AutoTutorial 知识库重建 ")
    print("==========================================")
    print(f"💾 数据库路径: {DB_PATH}")
    print(f"📂 数据来源: {SOURCE_DIR}")

    client = chromadb.PersistentClient(path=DB_PATH)

    try:
        client.delete_collection(COLLECTION_NAME)
        print("🧹 已清理旧集合")
    except Exception:
        pass

    print("🔌 连接阿里云 Embedding API...")
    collection = client.create_collection(name=COLLECTION_NAME)

    files = scan_files(SOURCE_DIR)
    print(f"📦 找到 {len(files)} 个文件，开始处理...")
    if not files:
        return

    total_chunks, _ = ingest_files(collection, files, ingest_root_dir=SOURCE_DIR, ingest_type="rebuild")
    print(f"\n🎉 重建完成！总片段: {total_chunks}")
    print(f"💾 数据库已保存至: {DB_PATH}")


def append_db() -> None:
    print("==========================================")
    print(" 🧩 AutoTutorial 知识库增量追加 ")
    print("==========================================")
    print(f"💾 数据库路径: {DB_PATH}")
    print(f"📂 增量目录: {ADD_DIR}")

    client = chromadb.PersistentClient(path=DB_PATH)
    collection = client.get_or_create_collection(name=COLLECTION_NAME)

    files = scan_files(ADD_DIR)
    print(f"📦 找到 {len(files)} 个增量文件，开始处理...")
    if not files:
        return

    total_chunks, ingested_files = ingest_files(collection, files, ingest_root_dir=ADD_DIR, ingest_type="append")
    print(f"\n🎉 增量追加完成！本次新增片段: {total_chunks}")
    print(f"💾 数据库路径: {DB_PATH}")

    os.makedirs(SOURCE_DIR, exist_ok=True)
    moved = 0

    for src_path in ingested_files:
        rel_path = os.path.relpath(src_path, ADD_DIR)
        dst_path = os.path.join(SOURCE_DIR, rel_path)
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)

        try:
            shutil.move(src_path, dst_path)
            moved += 1
        except Exception as e:
            print(f"⚠️ 归档失败 {src_path} -> {dst_path}: {e}")

    print(f"📂 已自动归档 {moved} 个文件到 knowledge_source。")

# ================= 主入口：运行时选择 =================

def main() -> None:
    setup_api_key()

    print("\n请选择操作模式：")
    print("1) 重建（清空并重建集合，读取 knowledge_source）")
    print("2) 添加（增量追加，读取 knowledge_add，成功后归档到 knowledge_source）")
    choice = input("请输入 1 或 2：").strip()

    if choice == "1":
        print("确认要进行重建？这将清空原有数据库，并重新读取 knowledge_source 中所有文件。")
        confirm = input("Y/N：").strip().lower()
        if confirm == "y":
            rebuild_db()
        else:
            print("已取消。")

    elif choice == "2":
        print("已选择增量追加，即将读取 knowledge_add……")
        append_db()

    else:
        print("❌ 无效输入，退出。")


if __name__ == "__main__":
    main()
