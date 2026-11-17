#!/usr/bin/env python3
"""방출고 민감도 — 평균장 런 3종(z = 2.5 / 5 / 20 m)을 비교한다.

세 런은 유동장·난수·침적 모델(DD_WALL_DEP=1, d_w=2 m)이 모두 같고 LDM 입력의 z 열만 다르다.
동기: 침적량은 Σ(Δt/dz) 로 결정되고(REPORT §5 검산 ρ=+0.96), λp 가 그 안으로 들어오는 통로는
'침적층 안에서의 평균 높이'(+0.71)다. 방출고를 침적층 안/위로 옮기면 이 통로가 어떻게 되는가?

출력: release_height.csv, 화면 요약
"""
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
COURT = {3, 6, 10}
R_MAIN = 150
NREL = 99988

CASES = [(2.5, '2026-09-17_leipzig-z25-meanflow'),
         (5.0, '2026-09-17_leipzig-dt01-wall-meanflow'),      # 기존 기준 런
         (20.0, '2026-09-17_leipzig-z20-meanflow')]

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read = ns['read_vtk']

S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
IDS = sorted(S)
KEEP = [i for i in IDS if i not in COURT]


def sp(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


def area90(px, py, cell=10.0):
    """침적 입자의 90 %를 담는 지면 면적 [m²] — dep_metrics.py 와 같은 정의."""
    if len(px) < 10: return 0.0
    ix = np.floor(px/cell).astype(int); iy = np.floor(py/cell).astype(int)
    _, c = np.unique(ix.astype(np.int64)*100000 + iy, return_counts=True)
    c = np.sort(c)[::-1]
    return float(np.searchsorted(np.cumsum(c), 0.9*c.sum()) + 1)*cell*cell


rows = []
for z0, d in CASES:
    f = os.path.join(RUNS, d, 'out', 'LDM_30000stp.vtk')
    if not os.path.exists(f):
        print(f'!! 없음: {f}'); continue
    p, F = read(f)
    src = F['p_type'].astype(int) - 100
    buf = F['buffer_type'].astype(int)
    for i in IDS:
        m = src == i
        g = m & (buf == 1); w = m & (buf == 2); air = m & (buf == 0)
        sx, sy = S[i]
        r = np.hypot(p[g, 0]-sx, p[g, 1]-sy)
        rows.append(dict(z_rel=z0, run=d, src=i, lam=lam[R_MAIN][i],
                         f_ground=g.sum()/NREL, f_wall=w.sum()/NREL,
                         f_dep=(g.sum()+w.sum())/NREL,
                         area90_m2=area90(p[g, 0], p[g, 1]),
                         f_r150=float((r <= 150).mean()) if g.sum() else np.nan,
                         z_air_med=float(np.median(p[air, 2])) if air.sum() else np.nan,
                         f_air_below2=float((p[air, 2] <= 2).mean()) if air.sum() else np.nan))

with open(os.path.join(HERE, 'release_height.csv'), 'w', newline='') as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
    wtr.writeheader(); wtr.writerows(rows)

print(f'{"z_rel":>6} {"침적률":>8} {"지면":>8} {"벽면":>8} {"벽/지면":>8} '
      f'{"ρ(f_dep,λp)":>12} {"ρ(지면,λp)":>11} {"ρ(A90,λp)":>10} {"공중z중앙":>9}')
for z0, d in CASES:
    R = [r for r in rows if r['z_rel'] == z0]
    if not R: continue
    K = [r for r in R if r['src'] not in COURT]
    L = [r['lam'] for r in K]
    fd = np.mean([r['f_dep'] for r in R]); fg = np.mean([r['f_ground'] for r in R])
    fw = np.mean([r['f_wall'] for r in R])
    print(f'{z0:6.1f} {fd:8.3f} {fg:8.3f} {fw:8.3f} {fw/max(fg,1e-9):8.3f} '
          f'{sp(L, [r["f_dep"] for r in K]):12.2f} {sp(L, [r["f_ground"] for r in K]):11.2f} '
          f'{sp(L, [r["area90_m2"] for r in K]):10.2f} '
          f'{np.nanmean([r["z_air_med"] for r in R]):9.1f}')

print('\n소스별 침적률')
print('  src   λp   ' + ''.join(f'{z0:>8.1f} m' for z0, _ in CASES) + '   z2.5/z5')
for i in IDS:
    v = [next((r['f_dep'] for r in rows if r['z_rel'] == z0 and r['src'] == i), np.nan)
         for z0, _ in CASES]
    mk = '*' if i in COURT else ' '
    ratio = v[0]/v[1] if len(v) > 1 and v[1] else np.nan
    print(f' {mk}S{i:<3}{lam[R_MAIN][i]:.3f} ' + ''.join(f'{x:10.3f}' for x in v) + f'{ratio:9.2f}')
