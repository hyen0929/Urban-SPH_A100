#!/usr/bin/env python3
"""평균장이 전이 유동을 어느 규모까지 재현하는가 — 셀 크기에 따른 침적 지도 상관.

발표자료 8번(5×5 이상화 배열)과 13번(Leipzig 실도시)이 서로 다른 결론처럼 보이던 문제를
같은 지표로 맞춰 푼다. 핵심: 선형 카운트 Pearson 은 소스 근처 고농도 셀이 지배해 항상 0.95+ 가
나온다. 순위상관(또는 log)으로 보면 두 자료 모두 **셀이 작아질수록 급격히 무너진다** — 즉
평균장은 덩어리는 맞히고 블록 규모 패턴은 잃는다는 같은 결론이다.

출력: meanflow_scale.csv
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
ARCH = os.path.join(RUNS, '..', 'archive', 'Ph.D. Defence (CaseCS)')

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def load(')[0], ns)
read = ns['read_vtk']


def deposited(path):
    p, F = read(path)
    m = (F['buffer_type'].astype(int) != 0) if 'buffer_type' in F \
        else (np.linalg.norm(F['velocity'], axis=1) == 0)
    return p[m][:, 0], p[m][:, 1]


def grid(px, py, lim, cell):
    nx = int((lim[1]-lim[0])/cell)+1; ny = int((lim[3]-lim[2])/cell)+1
    ix = np.clip(((px-lim[0])/cell).astype(int), 0, nx-1)
    iy = np.clip(((py-lim[2])/cell).astype(int), 0, ny-1)
    G = np.zeros((nx, ny)); np.add.at(G, (ix, iy), 1.0)
    return G.ravel()


def spear(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


CASES = [
    # (이름, 평균장 파일, 전이 파일, 영역, 셀 목록, 건물 크기[m])
    ('5x5 array', os.path.join(ARCH, '*CaseCS1(LDM averaged)', 'LDM_20000stp.vtk'),
     os.path.join(ARCH, '*CaseCS1(LES-LDM transient)', 'LDM_30000stp.vtk'),
     [40., 255., -60., 70.], [2., 3., 5., 7.5, 10., 15., 20.], 10.0),
    ('Leipzig', os.path.join(RUNS, '2026-10-23_leipzig-s13-most-meanflow', 'out', 'LDM_30000stp.vtk'),
     os.path.join(RUNS, '2026-10-20_leipzig-s13-most-wall', 'out', 'LDM_30000stp.vtk'),
     [0., 1320., -660., 660.], [2.5, 5., 10., 20., 40., 80.], 30.0),
]

rows = []
for name, fm, ft, lim, cells, bsize in CASES:
    if not (os.path.exists(fm) and os.path.exists(ft)):
        print(f'!! 없음: {name}'); continue
    mx, my = deposited(fm); tx, ty = deposited(ft)
    print(f'{name}: 평균장 {len(mx):,} / 전이 {len(tx):,} 침적 입자')
    for c in cells:
        g1, g2 = grid(mx, my, lim, c), grid(tx, ty, lim, c)
        ok = (g1+g2) > 0
        rows.append(dict(case=name, cell_m=c, cell_per_building=c/bsize, n_cells=int(ok.sum()),
                         pearson_linear=float(np.corrcoef(g1[ok], g2[ok])[0, 1]),
                         pearson_log=float(np.corrcoef(np.log10(g1[ok]+1), np.log10(g2[ok]+1))[0, 1]),
                         spearman=spear(g1[ok], g2[ok])))
        print(f'   셀 {c:5.1f} m  n={ok.sum():6d}  선형 {rows[-1]["pearson_linear"]:.3f}  '
              f'log {rows[-1]["pearson_log"]:.3f}  순위 {rows[-1]["spearman"]:.3f}')

with open(os.path.join(HERE, 'meanflow_scale.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('\nmeanflow_scale.csv 저장')
