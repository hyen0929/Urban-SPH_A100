#!/usr/bin/env python3
"""소스별 지면 침적 격자 — 노출 격자(ped_grid_fine.npz)와 같은 규격(5 m, 같은 lim)으로 만든다.

기존 grid_dep.npz 는 13개 소스를 합친 10 m 격자라 '소스별로 색을 달리한 지도'를 못 그린다.
여기서는 전이 런의 마지막 LDM 스냅샷에서 침적 입자를 소스별로 나눠 담는다.

출력: grid_dep_src.npz  (S1..S13, cell, lim, nrel)
"""
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
RUN = sys.argv[1] if len(sys.argv) > 1 else '2026-09-15_leipzig-dt01-sponge-uy'
STEP = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
CELL = 5.0
LIM = [0.0, 1320.0, -660.0, 660.0]

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def load(')[0], ns)
read = ns['read_vtk']
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}

p, F = read(os.path.join(RUNS, RUN, 'out', f'LDM_{STEP}stp.vtk'))
src = F['p_type'].astype(int) - 100
dep = (F['buffer_type'].astype(int) == 1) if 'buffer_type' in F \
    else (np.linalg.norm(F['velocity'], axis=1) == 0)      # 지면 침적만
nx = int((LIM[1]-LIM[0])/CELL)+1; ny = int((LIM[3]-LIM[2])/CELL)+1
ix = np.clip(((p[:, 0]-LIM[0])/CELL).astype(int), 0, nx-1)
iy = np.clip(((p[:, 1]-LIM[2])/CELL).astype(int), 0, ny-1)

out, nrel = {}, []
for i in sorted(S):
    m = dep & (src == i)
    G = np.zeros((nx, ny), np.float32)
    if m.any(): np.add.at(G, (ix[m], iy[m]), 1.0)
    out[f'S{i}'] = G
    nrel.append(int((src == i).sum()))
    print(f'  S{i:<3} 침적 {int(m.sum()):6d}')
out['cell'] = np.array(CELL); out['lim'] = np.array(LIM); out['nrel'] = np.array(nrel)
np.savez_compressed(os.path.join(HERE, 'grid_dep_src.npz'), **out)
print(f'{CELL:g} m 격자 {nx}×{ny}, 총 침적 {int(dep.sum()):,} → grid_dep_src.npz')
