#!/usr/bin/env python3
"""보행자 높이 노출 — 침적량보다 도시형태와 직접 얽히는 양.

정의: 공중에 있는 LDM 입자가 보행자층(z <= PED_Z)에 머무는 시간을 수평 격자에 적분한다.
  E(x,y) = Σ_t N(x,y,t)·Δt / (셀부피)       [입자·s/m³]
  소스별 정규화 노출 D = Σ_cells E·A_cell·PED_Z / N_released = 방출 입자 1개당 보행자층 체류시간 [s]
침적한 입자(buffer_type!=0 또는 속도 0)는 공기 중에 없으므로 제외한다.

출력: ped_exposure.csv (소스별), ped_grid.npz (지도), 화면 요약
"""
import csv, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
RUN = os.path.join(RUNS, sys.argv[1] if len(sys.argv) > 1 else '2026-09-15_leipzig-dt01-sponge-uy')
CELL, EVERY, DT = 10.0, 100, 1.0                  # 격자 [m], step 간격(1 s)
# 고도대 민감도: (하한, 상한). 첫 번째가 기본 보행자층.
BANDS = [(0.0, 2.0), (0.0, 1.5), (0.0, 4.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0)]
COURT = {3, 6, 10}

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read = ns['read_vtk']
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}


def spearman(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


steps = sorted(int(re.search(r'_(\d+)stp', f).group(1))
               for f in os.listdir(os.path.join(RUN, 'out')) if f.startswith('LDM_'))
steps = [s for s in steps if s >= 15000 and s % EVERY == 0]
lim = None
acc = {(i, b): None for i in S for b in range(len(BANDS))}
nrel = {i: 0 for i in S}
for s in steps:
    p, F = read(os.path.join(RUN, 'out', f'LDM_{s}stp.vtk'))
    src = F['p_type'].astype(int) - 100
    air = (F['buffer_type'].astype(int) == 0) if 'buffer_type' in F else (np.linalg.norm(F['velocity'], axis=1) != 0)
    if lim is None:
        lim = [0.0, 1320.0, -660.0, 660.0]
        nx = int((lim[1]-lim[0])/CELL)+1; ny = int((lim[3]-lim[2])/CELL)+1
        for k in acc: acc[k] = np.zeros((nx, ny), dtype=np.float32)
    ix = np.clip(((p[:, 0]-lim[0])/CELL).astype(int), 0, nx-1)
    iy = np.clip(((p[:, 1]-lim[2])/CELL).astype(int), 0, ny-1)
    for b, (z0b, z1b) in enumerate(BANDS):
        ped = air & (p[:, 2] > z0b) & (p[:, 2] <= z1b)
        for i in S:
            m = ped & (src == i)
            if m.any(): np.add.at(acc[(i, b)], (ix[m], iy[m]), DT)
    for i in S:
        nrel[i] = max(nrel[i], int((src == i).sum()))


def metrics(A, i, thick):
    E = A/(CELL*CELL*thick)
    tot = float(A.sum())
    gx = lim[0] + (np.arange(A.shape[0])+.5)*CELL; gy = lim[2] + (np.arange(A.shape[1])+.5)*CELL
    GX, GY = np.meshgrid(gx, gy, indexing='ij')
    r = np.hypot(GX-S[i][0], GY-S[i][1])
    return dict(exposure_s_per_particle=tot/max(nrel[i], 1), peak_s_per_m3=float(E.max()),
                area_1pct_m2=float((E > 0.01*E.max()).sum()*CELL*CELL),
                frac_within_150m=float(A[r <= 150].sum()/max(tot, 1e-9)),
                r_mean_m=float((A*r).sum()/max(tot, 1e-9)))

ids = sorted(S)
keep = np.array([i not in COURT for i in ids])
x75 = np.array([lam[75][i] for i in ids]); x150 = np.array([lam[150][i] for i in ids])
rows = []
print(f"{'고도대 [m]':>12} {'노출량':>8} {'최대농도':>9} {'노출면적':>9} {'150m 갇힘':>10} {'평균거리':>9}   (λp(75) Spearman, 중정 제외 10점)")
for b, (z0b, z1b) in enumerate(BANDS):
    M = [metrics(acc[(i, b)], i, z1b-z0b) for i in ids]
    for i, m in zip(ids, M):
        rows.append(dict(src=i, band=f'{z0b:g}-{z1b:g}', lambda_p75=lam[75][i], lambda_p150=lam[150][i], **m))
    def c(key):
        y = np.array([m[key] for m in M]); return spearman(x75[keep], y[keep])
    print(f"{z0b:5.1f}-{z1b:<5.1f} {c('exposure_s_per_particle'):8.2f} {c('peak_s_per_m3'):9.2f} "
          f"{c('area_1pct_m2'):9.2f} {c('frac_within_150m'):10.2f} {c('r_mean_m'):9.2f}")
with open(os.path.join(HERE, os.environ.get('PED_BANDS', 'ped_exposure_bands.csv')), 'w', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w_.writeheader(); w_.writerows(rows)
tot_all = {b: sum(float(acc[(i, b)].sum()) for i in ids) for b in range(len(BANDS))}
print("\n고도대별 총 체류시간 [입자·s] 과 단위두께당:")
for b, (z0b, z1b) in enumerate(BANDS):
    print(f"  {z0b:4.1f}-{z1b:<4.1f} m : {tot_all[b]:12,.0f}  ({tot_all[b]/(z1b-z0b):11,.0f} /m)")
