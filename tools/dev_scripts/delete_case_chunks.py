#!/usr/bin/env python3
"""删除案例文档的旧 chunks"""

import os
import sys
import chromadb

# 添加父目录到路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, ROOT_DIR)

DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
COLLECTION_NAME = "abacus_knowledge"

def delete_case_chunks():
    """删除案例文档的所有 chunks"""
    client = chromadb.PersistentClient(path=DB_PATH)
    collection = client.get_collection(name=COLLECTION_NAME)

    # 获取所有数据
    all_data = collection.get()

    # 案例文件名
    case_files = [
        "band_Si_diamond.md",
        "dos_MgO.md",
        "elastic_Al_FCC_org.md",
        "elastic_Al_FCC_relax.md"
    ]

    # 找出需要删除的 IDs
    ids_to_delete = []
    for i, metadata in enumerate(all_data['metadatas']):
        source = metadata.get('source', '')
        for case_file in case_files:
            if case_file in source:
                ids_to_delete.append(all_data['ids'][i])
                break

    if ids_to_delete:
        print(f"找到 {len(ids_to_delete)} 个旧案例 chunks")
        collection.delete(ids=ids_to_delete)
        print(f"删除完成")
    else:
        print("未找到旧案例 chunks")


if __name__ == "__main__":
    delete_case_chunks()
