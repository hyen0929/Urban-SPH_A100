#!/usr/bin/env python3
"""runs/timeline_dt01/ 재생성 — 세 분기 런의 유효 구간만 심볼릭 링크 + CSV 병합. README.md 참조."""
import os, re, csv, glob
R = os.path.dirname(os.path.abspath(__file__)) + '/..'
T = os.path.join(R, 'timeline_dt01')
SEG = [('2026-09-14_leipzig-dt01-full', 0, 11000),
       ('2026-09-15_leipzig-dt01-sponge', 11001, 12000),
       ('2026-09-15_leipzig-dt01-sponge-uy', 12001, 30000)]
KIND = {'slab': r'fluid_condition_(\d+)stp\.vtk', 'vort': r'vorticity_(\d+)stp\.vtk',
        'full': r'fluid_(\d+)stp\.vtk', 'LDM': r'LDM_(\d+)stp\.vtk'}
for k in KIND: os.makedirs(os.path.join(T, k), exist_ok=True)
for d, lo, hi in SEG:
    for f in sorted(glob.glob(os.path.join(R, d, 'out', '*.vtk'))):
        b = os.path.basename(f)
        for k, pat in KIND.items():
            m = re.fullmatch(pat, b)
            if not m: continue
            s = int(m.group(1))
            if not (lo <= s <= hi): continue
            link = os.path.join(T, k, b)
            if not os.path.islink(link): os.symlink(os.path.relpath(f, os.path.dirname(link)), link)
rows = []
for d, lo, hi in SEG:
    for r in csv.DictReader(open(os.path.join(R, d, 'out', 'PPE_diag_explicit.csv'))):
        try: s = int(r['step'])
        except ValueError: continue
        if lo <= s <= hi: rows.append((s, d, r))
rows.sort()
cols = ['step','time','bi_max','resid_max','div_un1_max','div_un1_avg','div_un1_plain_max','div_un1_plain_avg']
with open(os.path.join(T, 'PPE_diag_timeline.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(cols + ['run'])
    for s, d, r in rows: w.writerow([r[c] for c in cols] + [d])
SUMK = re.compile(rb"^(\d+)\s+tu_max=([\d.eE+-]+)\s+.*?uz_max=([\d.eE+-]+)\s+ftotal_max=([\d.eE+-]+)(?:.*\n)*?"
                  rb".*?\xce\xa3k_turb=([\d.eE+-]+)\s+\xce\xa3\xc2\xbd\|u\|\xc2\xb2=([\d.eE+-]+)", re.M)
er = []
for d, lo, hi in SEG:
    for m in SUMK.finditer(open(os.path.join(R, d, 'run.log'), 'rb').read()):
        s = int(m.group(1))
        if lo <= s <= hi: er.append((s, *[float(m.group(i)) for i in (2,3,4,5,6)], d))
er.sort()
with open(os.path.join(T, 'domain_timeline.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['step','time','tu_max','uz_max','ftotal_max','sum_k_turb','sum_half_u2','run'])
    for s, tu, uz, ft, k, ke, d in er: w.writerow([s, round(s*0.01, 2), tu, uz, ft, k, ke, d])
print('timeline 재생성 완료')
