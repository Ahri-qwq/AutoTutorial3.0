import sys
import re
sys.stdout.reconfigure(encoding='utf-8')

with open('_workspace/20260302_隐式溶剂模型/process/03_draft_full.md', encoding='utf-8') as f:
    lines = f.readlines()

print("=== 超长句检查（>80字的正文行）===")
found = False
for i, line in enumerate(lines, 1):
    ls = line.strip()
    if not ls:
        continue
    if ls[0] in '#|->*':
        continue
    if ls.startswith('```') or ls.startswith('`'):
        continue
    if len(ls) > 80:
        found = True
        print('行%d(%d字): %s...' % (i, len(ls), ls[:60]))
if not found:
    print("无超长句")

print("\n=== AI腔检查 ===")
ai_patterns = ['在当今', '随着.*的发展', '值得注意的是', '综上所述', '不得不说']
found2 = False
for i, line in enumerate(lines, 1):
    for pat in ai_patterns:
        if re.search(pat, line):
            found2 = True
            print('行%d: %s' % (i, line.strip()))
if not found2:
    print("未发现AI腔表达")
