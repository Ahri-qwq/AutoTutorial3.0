import re

with open('03_draft_full.md', 'r', encoding='utf-8') as f:
    lines = f.readlines()

print('=== 超过50个中文字的正文行（非代码块、非表格）===')
in_code = False
for i, line in enumerate(lines, 1):
    stripped = line.strip()
    if stripped.startswith('```'):
        in_code = not in_code
    if in_code:
        continue
    if stripped.startswith('|') or stripped.startswith('#') or stripped.startswith('$'):
        continue
    cn_chars = len(re.findall(r'[\u4e00-\u9fff]', stripped))
    if cn_chars > 50:
        print(f'Line {i:3d} ({cn_chars}chars): {stripped[:100]}')
