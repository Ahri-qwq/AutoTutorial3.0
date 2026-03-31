#!/usr/bin/env python3
"""检查数据库中的文档数量"""

import os
import sys
import chromadb

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, ROOT_DIR)

DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
COLLECTION_NAME = "abacus_knowledge"

def check_db():
    client = chromadb.PersistentClient(path=DB_PATH)
    collection = client.get_collection(name=COLLECTION_NAME)

    all_data = collection.get()
    total = len(all_data['ids'])

    print(f"数据库总片段数: {total}")

    # 统计案例文档
    case_files = ["band_Si_diamond.md", "dos_MgO.md", "elastic_Al_FCC_org.md", "elastic_Al_FCC_relax.md"]
    case_count = 0

    for metadata in all_data['metadatas']:
        source = metadata.get('source', '')
        for case_file in case_files:
            if case_file in source:
                case_count += 1
                break

    print(f"案例文档片段数: {case_count}")
    print(f"其他文档片段数: {total - case_count}")

if __name__ == "__main__":
    check_db()
