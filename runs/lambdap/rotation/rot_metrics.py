#!/usr/bin/env python3
"""풍향 회전 비교 — 0/90/180(/270)° 의 소스별 침적 지표와 λp 상관.

지표 정의·침적 판정은 deposition/dep_metrics.py 를 그대로 재사용한다.
도시를 회전시켰으므로 소스 좌표도 (657.5, 0) 기준으로 같이 회전시킨다.
λp 는 소스와 주변 건물이 함께 돌아 회전 불변이므로 원본 값을 쓴다.
"""
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEP  = os.path.normpath(os.path.join(HERE, '..', 'deposition'))
MORPH= os.path.normpath(os.path.join(HERE, '..'))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
ns = {'__file__': os.path.join(DEP, 'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
load, footprint, spearman = ns['load'], ns['footprint'], ns['spearman']

CX, CY = 657.5, 0.0
COURT = {3, 6, 10}
ANGLES = [(0,   '2026-09-15_leipzig-dt01-sponge-uy'),
          (90,  '2026-09-23_leipzig-rot90'),
          (180, '2026-09-23_leipzig-rot180'),
          (270, '2026-09-26_leipzig-rot270')]

S0 = {int(r['id']): (float(r['x']), float(r['y'])) for r in
      csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}

def rot(xy, deg):
    t = np.radians(deg); c, s = np.cos(t), np.sin(t)
    dx, dy = xy[0]-CX, xy[1]-CY
    return (CX + dx*c - dy*s, CY + dx*s + dy*c)

IDS = sorted(S0); keep = np.array([i not in COURT for i in IDS])
KEYS = ['f_dep', 'area90_m2', 'f_r150', 'shift_m', 'elong']
rows, corr = [], []
for deg, d in ANGLES:
    rd = os.path.join(RUNS, d)
    if not os.path.exists(os.path.join(rd, 'out', 'LDM_30000stp.vtk')):
        print(f'{deg:4d}°  없음 ({d})'); continue
    p, src, dep = load(rd)
    M = {}
    for i in IDS:
        m = src == i
        f = footprint(p[m & dep][:, :2], np.array(rot(S0[i], deg)))
        f['f_dep'] = float((m & dep).sum()/max(m.sum(), 1))
        M[i] = f
        rows.append(dict(deg=deg, src=i, lambda_p150=lam[150][i], **{k: f[k] for k in KEYS}))
    c = {'deg': deg, 'dep_total_pct': 100*dep.mean()}
    for k in KEYS:
        y = np.array([M[i][k] for i in IDS])
        c[k] = spearman(np.array([lam[150][i] for i in IDS])[keep], y[keep])
    corr.append(c)
    print(f'{deg:4d}°  침적 {c["dep_total_pct"]:5.2f}%   ' +
          '  '.join(f'{k} {c[k]:+.2f}' for k in KEYS))

with open(os.path.join(HERE, 'rot_metrics.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(os.path.join(HERE, 'rot_corr.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(corr[0].keys())); w.writeheader(); w.writerows(corr)
print(f'\n저장: rot_metrics.csv ({len(rows)}행), rot_corr.csv')
