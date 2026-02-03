#!/usr/bin/env python3
"""
RAG检索工具 - 命令行版本
用于从知识库检索相关文档

使用方法：
    python tools/retriever.py --query "查询词" --top_k 10
"""

import os
import sys
import argparse
import chromadb
import dashscope
from dashscope import TextEmbedding
import yaml

# 添加父目录到路径，以便导入src模块
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, ROOT_DIR)

DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")

# 初始化 API Key
try:
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
            api_key = config.get("llm", {}).get("api_key") or config.get("api_key")
            if api_key:
                dashscope.api_key = api_key
except:
    pass

if not dashscope.api_key:
    dashscope.api_key = os.getenv("DASHSCOPE_API_KEY", "")


class CommandLineRetriever:
    """命令行版本的检索器"""

    def __init__(self):
        """初始化ChromaDB客户端"""
        try:
            self.client = chromadb.PersistentClient(path=DB_PATH)
            self.collection = self.client.get_collection(name="abacus_knowledge")
        except Exception as e:
            print(f"[错误] 无法连接到知识库: {e}", file=sys.stderr)
            print(f"[提示] 请确保知识库已初始化，路径: {DB_PATH}", file=sys.stderr)
            sys.exit(1)

    def _get_embedding(self, text: str):
        """计算文本的向量表示"""
        try:
            resp = TextEmbedding.call(
                model=TextEmbedding.Models.text_embedding_v3,
                input=[text]
            )
            if resp.status_code == 200:
                return resp.output["embeddings"][0]["embedding"]
            else:
                print(f"[错误] Embedding失败: {resp.message}", file=sys.stderr)
                return None
        except Exception as e:
            print(f"[错误] 网络错误: {e}", file=sys.stderr)
            return None

    def search(self, query: str, top_k: int = 5):
        """
        执行检索

        Args:
            query: 查询词
            top_k: 返回文档数量

        Returns:
            检索结果列表，每个元素是 (文档内容, 来源)
        """
        # 计算查询向量
        query_vec = self._get_embedding(query)
        if not query_vec:
            return []

        # 执行检索
        try:
            results = self.collection.query(
                query_embeddings=[query_vec],
                n_results=top_k
            )
        except Exception as e:
            print(f"[错误] 检索失败: {e}", file=sys.stderr)
            return []

        # 解析结果
        retrieved_docs = []
        if results and results['documents']:
            docs = results['documents'][0]
            metas = results['metadatas'][0]
            for i, doc in enumerate(docs):
                source = metas[i].get('source', 'Unknown')
                retrieved_docs.append((doc, source))

        return retrieved_docs


def format_output(docs, query):
    """格式化输出结果"""
    # 设置输出编码为UTF-8
    import sys
    if sys.platform == 'win32':
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    print("=" * 80)
    print(f"检索查询: {query}")
    print(f"返回文档数: {len(docs)}")
    print("=" * 80)
    print()

    if not docs:
        print("[提示] 未找到相关文档")
        return

    for i, (content, source) in enumerate(docs, 1):
        print(f"文档{i}:")
        print(f"来源: {source}")
        print(f"内容:")
        # 处理可能的特殊字符
        try:
            print(content)
        except UnicodeEncodeError:
            print(content.encode('utf-8', errors='replace').decode('utf-8'))
        print()
        print("-" * 80)
        print()


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="RAG检索工具 - 从知识库检索相关文档",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    python tools/retriever.py --query "弹性常数计算 物理原理" --top_k 10
    python tools/retriever.py --query "ABACUS INPUT参数" --top_k 5
        """
    )

    parser.add_argument(
        "--query",
        type=str,
        required=True,
        help="检索查询词"
    )

    parser.add_argument(
        "--top_k",
        type=int,
        default=5,
        help="返回文档数量 (默认: 5)"
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help="安静模式，只输出文档内容，不输出格式化信息"
    )

    args = parser.parse_args()

    # 检查API Key
    if not dashscope.api_key:
        print("[错误] 未找到DashScope API Key", file=sys.stderr)
        print("[提示] 请在config.yaml中配置api_key，或设置环境变量DASHSCOPE_API_KEY", file=sys.stderr)
        sys.exit(1)

    # 初始化检索器
    retriever = CommandLineRetriever()

    # 执行检索
    docs = retriever.search(args.query, args.top_k)

    # 输出结果
    if args.quiet:
        # 安静模式：只输出文档内容
        for content, source in docs:
            print(f"[来源: {source}]")
            print(content)
            print()
    else:
        # 正常模式：格式化输出
        format_output(docs, args.query)


if __name__ == "__main__":
    main()
