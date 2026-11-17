#!/usr/bin/env python3
"""보행자층(z<2 m) 노출 지도만 고해상도로 다시 적분한다.

`ped_exposure.py` 는 고도대 6개 × 소스 13개를 10 m 격자로 뽑는다. 지도 그림에서는
10 m 셀의 계단이 보이므로, 여기서는 **기본 밴드 하나만** 세밀한 격자로 다시 쌓는다.
지표(면적·최대·총량)는 셀 크기에 의존하므로 건드리지 않고 그림용으로만 쓴다.

  python3 ped_grid_fine.py [run] [cell_m]     기본 5 m → ped_grid_fine.npz
"""
import csv, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
RUN = os.path.join(RUNS, sys.argv[1] if len(sys.argv) > 1 else '2026-10-20_leipzig-s13-most-wall')
CELL = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
EVERY, DT, PED_Z = 100, 1.0, 2.0
LIM = [0.0, 1320.0, -660.0, 660.0]

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read = ns['read_vtk']
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}

steps = sorted(int(re.search(r'_(\d+)stp', f).group(1))
               for f in os.listdir(os.path.join(RUN, 'out')) if f.startswith('LDM_'))
steps = [s for s in steps if s >= 15000 and s % EVERY == 0]
nx = int((LIM[1]-LIM[0])/CELL)+1; ny = int((LIM[3]-LIM[2])/CELL)+1
acc = {i: np.zeros((nx, ny), dtype=np.float32) for i in S}
nrel = {i: 0 for i in S}
for n, s in enumerate(steps):
    p, F = read(os.path.join(RUN, 'out', f'LDM_{s}stp.vtk'))
    src = F['p_type'].astype(int) - 100
    air = (F['buffer_type'].astype(int) == 0) if 'buffer_type' in F \
        else (np.linalg.norm(F['velocity'], axis=1) != 0)
    ped = air & (p[:, 2] > 0) & (p[:, 2] <= PED_Z)
    ix = np.clip(((p[:, 0]-LIM[0])/CELL).astype(int), 0, nx-1)
    iy = np.clip(((p[:, 1]-LIM[2])/CELL).astype(int), 0, ny-1)
    for i in S:
        m = ped & (src == i)
        if m.any(): np.add.at(acc[i], (ix[m], iy[m]), DT)
        nrel[i] = max(nrel[i], int((src == i).sum()))
    if n % 30 == 0: print(f'  {n+1}/{len(steps)}', flush=True)

out = {f'S{i}': acc[i] for i in S}
out['cell'] = np.array(CELL); out['lim'] = np.array(LIM)
out['nrel'] = np.array([nrel[i] for i in sorted(S)])
np.savez_compressed(os.path.join(HERE, os.environ.get('PGF_OUT', 'ped_grid_fine.npz')), **out)
tot = sum(acc[i].sum() for i in S)
print(f'{CELL:g} m 격자 {nx}×{ny}, 체류 {tot:.0f} 입자·s, 파일 ped_grid_fine.npz')
